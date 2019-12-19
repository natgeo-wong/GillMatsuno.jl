function calcQ0(x,A,L)

    if abs(x) < L
        return A * cos(pi*x/(2*L))
    else
        return 0
    end

end

function calcq0(x,A,α,k,L)

    if x <= -L
        return 0
    elseif abs(x) < L
        return -A / (α^2+k^2) * (α*cos(k*x) + k*(sin(k*x) + exp(-α*(x+L))))
    else
        return -A*k / (α^2+k^2) * (1+exp(-2*α*L)) * exp(α*(L-x))
    end

end

function calcq2(x,A,α,k,L)

    if x >= L
        return 0
    elseif abs(x) < L
        return A / ((3*α)^2+k^2) * (-3*α*cos(k*x) + k*(sin(k*x) - exp(3*α*(x-L))))
    else
        return -A*k / ((3*α)^2+k^2) * (1+exp(-6*α*L)) * exp(3*α*(L+x))
    end

end

function calcu(x,y,c,β,A,α,k,L)

    return c/2 * (calcq0(x,A,α,k,L) + calcq2(x,A,α,k,L) * (2*β*y^2/c - 3)) * exp(-β*y^2/(2*c))

end

function calcv(x,y,c,β,A,α,k,L)

    return c*y * (calcQ0(x,A,L) + calcq2(x,A,α,k,L) * (4*α/c)) * exp(-β*y^2/(2*c))

end

function calcϕ(x,y,c,β,A,α,k,L)

    return c^2/2 * (calcq0(x,A,α,k,L) + calcq2(x,A,α,k,L) * (2*β*y^2/c + 1)) * exp(-β*y^2/(2*c))

end

function calcw(x,y,c,β,A,α,k,L,H)

    return H*sqrt(c*β/2) *
            (2*calcQ0(x,A,L) + α*calcq0(x,A,α,k,L) + α*calcq2(x,A,α,k,L)*(2*β*y^2/c + 1)) *
            exp(-β*y^2/(2*c))

end

function GManalytic(xvec,yvec,g,β,A,α,L,H)

    c = sqrt(g*H); k = pi / (2*L)
    nx = length(xvec); ny = length(yvec);
    u = zeros(ny,nx); v = zeros(ny,nx); ϕ = zeros(ny,nx); w = zeros(ny,nx);

    for ii = 1 : nx
        for jj = 1 : ny
            u[jj,ii] = calcu(xvec[ii],yvec[jj],c,β,A,α,k,L)
            v[jj,ii] = calcv(xvec[ii],yvec[jj],c,β,A,α,k,L)
            ϕ[jj,ii] = calcϕ(xvec[ii],yvec[jj],c,β,A,α,k,L)
            w[jj,ii] = calcw(xvec[ii],yvec[jj],c,β,A,α,k,L,H)
        end
    end

    return u,v,ϕ

end
