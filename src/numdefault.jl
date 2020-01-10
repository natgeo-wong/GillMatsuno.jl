"""
The Gill Matsuno circulation is the steady-state solution to a constant forcing of damped
shallow-water equations on a beta-plane.  When the shallow-water equations are subject to forcing and dissipative effects, the equations are:

    u_t + αu + ϕ_x - βyv   = 0
    v_t + αu + ϕ_y + βyu   = 0
    ϕ_t + αϕ + gH(u_x+v_y) = -Q

where α is the damping coefficient and Q is the heat flux term.

This script decomposes the above equations into their numerical form.  We take as default
initial conditions u = v = ϕ = 0.
"""

function ufield!(
    un::Array, u::Array, ϕ::Array, v::Array, y::Array, δt::Real, δx::Real;
    parameters::Array=[0.1,0.5,1,1], nsize::Array
)

    α,β,g,H = parameters; nx,ny = nsize;

    for jj = 1 : ny, ii = 1 : nx

        if ii != 1 && ii != nx; ϕx = (ϕ[ii+1,jj] - ϕ[ii-1,jj]) / (2*δx)
        elseif ii == 1;         ϕx = -3*ϕ[ii,jj] + 4*ϕ[ii+1,jj] - ϕ[ii+2,jj];
        else;                   ϕx =  3*ϕ[ii,jj] - 4*ϕ[ii-1,jj] + ϕ[ii-2,jj];
        end

        un[ii,jj] = u[ii,jj] - δt * (α*u[ii,jj] - β*y[ii,jj]*v[ii,jj] + ϕx);

    end

    return

end

function vfield!(
    vn::Array, v::Array, ϕ::Array, u::Array, y::Array, δt::Real, δy::Real;
    parameters::Array=[0.1,0.5,1,1], nsize::Array
)

    α,β,g,H = parameters; nx,ny = nsize;

    for jj = 1 : ny, ii = 1 : nx

        if jj != 1 && jj != ny; ϕy = (ϕ[ii,jj+1] - ϕ[ii,jj-1]) / (2*δy)
        elseif jj == 1;         ϕy = -3*ϕ[ii,jj] + 4*ϕ[ii,jj+1] - ϕ[ii,jj+2];
        else;                   ϕy =  3*ϕ[ii,jj] - 4*ϕ[ii,jj-1] + ϕ[ii,jj-2];
        end

        vn[ii,jj] = v[ii,jj] - δt * (α*v[ii,jj] + β*y[ii,jj]*u[ii,jj] + ϕy);

    end

    return

end

function ϕfield!(
    ϕn::Array, ϕ::Array, u::Array, v::Array, Q::Array, δt::Real, δx::Real, δy::Real;
    parameters::Array=[0.1,0.5,1,1], nsize::Array
)

    α,β,g,H = parameters; nx,ny = nsize;

    for jj = 1 : ny, ii = 1 : nx

        if jj != 1 && jj != ny; vy = (v[ii,jj+1] - v[ii,jj-1]) / (2*δy)
        elseif jj == 1;         vy = -3*v[ii,jj] + 4*v[ii,jj+1] - v[ii,jj+2];
        else;                   vy =  3*v[ii,jj] - 4*v[ii,jj-1] + v[ii,jj-2];
        end

        if ii != 1 && ii != nx; ux = (u[ii+1,jj] - u[ii-1,jj]) / (2*δx)
        elseif ii == 1;         ux = -3*u[ii,jj] + 4*u[ii+1,jj] - u[ii+2,jj];
        else;                   ux =  3*u[ii,jj] - 4*u[ii-1,jj] + u[ii-2,jj];
        end

        ϕn[ii,jj] = ϕ[ii,jj] - δt * (α*ϕ[ii,jj] + g*H*(ux+vy) + Q[ii,jj]);

    end

    return

end

function Qfield!(Q::Array, A::Real, L::Real, xvec::Array, yvec::Array, nsize::Array)

    nx,ny = nsize;

    for jj = 1 : ny
        for ii = 1 : nx
            if abs(xvec[ii]) < L
                  Q[ii,jj] = A * cos.(pi*xvec[ii]*0.5/L) .* exp.(-(yvec[jj].^2)/4);
            else; Q[ii,jj] = 0;
            end
        end
    end

    return

end

function GMcalc(;xmin::Real=-25.0, xmax::Real=50.0, δx::Real=0.1,
                ymin::Real=-10.0, ymax::Real=10.0, δy::Real=0.1,
                nt::Integer=5000, δt::Real=0.001,
                A::Real=1, L::Real=2, α::Real=0.1, β::Real=0.5, g::Real=1.0, H::Real=1.0)

    xvec = convert(Array,xmin:δx:xmax); nx = length(xvec);
    yvec = convert(Array,ymin:δy:ymax); ny = length(yvec); nsize = [nx,ny];
    x = repeat(xvec,1,ny); y = convert(Array,transpose(repeat(yvec,1,nx)));
    Q = zeros(nx,ny); Qfield!(Q,A,L,xvec,yvec,nsize);
    ϕ = zeros(nx,ny); u = zeros(nx,ny); v = zeros(nx,ny);
    ϕn = zeros(nx,ny); un = zeros(nx,ny); vn = zeros(nx,ny);

    for tt = 1 : nt
        ϕfield!(ϕn,ϕ,u,v,Q,δt,δx,δy,parameters=[α,β,g,H],nsize=nsize)
        ufield!(un,u,ϕ,v,y,δt,δx,parameters=[α,β,g,H],nsize=nsize)
        vfield!(vn,v,ϕ,u,y,δt,δy,parameters=[α,β,g,H],nsize=nsize)
        ϕ .= ϕn; u .= un; v .= vn;
    end

    return ϕ,u,v

end
