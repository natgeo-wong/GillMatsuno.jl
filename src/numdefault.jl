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

function ufield(u::Array,ϕ::Array,v::Array,y::Array,
                δt::Real,δx::Real;parameters::Array=[0.1,0.5,1,1])

    α,β,g,H = parameters;

    ϕ_x = (circshift(ϕ,-1) - circshift(ϕ,1)) / (2*δx);
    ϕ_x[1,:] .= 0; ϕ_x[end,:] .= 0;
    # ϕ_x[1,:] .= 3*(ϕ_x[3,:]-ϕ_x[2,:]) - ϕ_x[4,:];
    # ϕ_x[end,:] .= 3*(ϕ_x[end-3,:]-ϕ_x[end-2,:]) - ϕ_x[end-4,:];

    return u - δt * (α*u - β.*y.*v + ϕ_x);

end

function vfield(v::Array,ϕ::Array,u::Array,y::Array,
                δt::Real,δy::Real;parameters::Array=[0.1,0.5,1,1])

    α,β,g,H = parameters;

    ϕ_y = (circshift(ϕ,(0,-1)) - circshift(ϕ,(0,1))) / (2*δy);
    ϕ_y[:,1] .= 0; ϕ_y[:,end] .= 0
    # ϕ_y[:,1] .= 3*(ϕ_y[:,3]-ϕ_y[:,2]) - ϕ_y[:,4];
    # ϕ_y[:,end] .= 3*(ϕ_y[:,end-3]-ϕ_y[:,end-2]) - ϕ_y[:,end-4];

    return v - δt * (α*v + β.*y.*u + ϕ_y);

end

function ϕfield(ϕ::Array,u::Array,v::Array,Q::Array,
                δt::Real,δx::Real,δy::Real;parameters::Array=[0.1,0.5,1,1])

    α,β,g,H = parameters;

    v_y = (circshift(v,(0,-1)) - circshift(v,(0,1))) / (2*δy);
    v_y[:,1] .= 0; v_y[:,end] .= 0
    # v_y[:,1] .= 3*(v_y[:,3]-v_y[:,2]) - v_y[:,4];
    # v_y[:,end] .= 3*(v_y[:,end-3]-v_y[:,end-2]) - v_y[:,end-4];

    u_x = (circshift(u,-1) - circshift(u,1)) / (2*δx);
    u_x[1,:] .= 0; u_x[end,:] .= 0;
    # u_x[1,:] .= 3*(u_x[3,:]-u_x[2,:]) - u_x[4,:];
    # u_x[end,:] .= 3*(u_x[end-3,:]-u_x[end-2,:]) - u_x[end-4,:];

    return ϕ - δt * (α*ϕ + g*H*(u_x+v_y) + Q);

end

function Qfield(Q::Array,A::Real,L::Real,xvec::Array,yvec::Array,nx::Integer,ny::Integer)

    for jj = 1 : ny
        for ii = 1 : nx
            if abs(xvec[ii]) < L
                  Q[ii,jj] = A * cos.(pi*xvec[ii]*0.5/L) .* exp.(-(yvec[jj].^2)/4);
            else; Q[ii,jj] = 0;
            end
        end
    end

    return Q

end

function GMcalc(;xmin::Real=-25.0,xmax::Real=50.0,δx::Real=0.1,
                ymin::Real=-10.0,ymax::Real=10.0,δy::Real=0.1,
                nt::Integer=2000,δt::Real=0.1,
                A::Real=1,L::Real=2,α::Real=0.1,β::Real=0.5,g::Real=1.0,H::Real=1.0)

    xvec = convert(Array,xmin:δx:xmax); nx = length(xvec);
    yvec = convert(Array,ymin:δy:ymax); ny = length(yvec);
    x = repeat(xvec,1,ny);
    y = convert(Array,transpose(repeat(yvec,1,nx)));
    Q = zeros(nx,ny); Q = Qfield(Q,A,L,xvec,yvec,nx,ny);
    ϕ = zeros(nx,ny); u = zeros(nx,ny); v = zeros(nx,ny);

    while δt > 0.5 * δx || δt > 0.5 * δy
        δt = δt / 2;
    end

    for ii = 1 : nt
        ϕn = ϕfield(ϕ,u,v,Q,δt,δx,δy,parameters=[α,β,g,H])
        un = ufield(u,ϕ,v,y,δt,δx,parameters=[α,β,g,H])
        vn = vfield(v,ϕ,u,y,δt,δy,parameters=[α,β,g,H])
        ϕ = deepcopy(ϕn); u = deepcopy(un); v = deepcopy(vn);
    end

    return ϕ,u,v

end
