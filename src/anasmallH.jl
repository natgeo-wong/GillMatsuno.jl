function calcuH(x,y,β,A,α,k)

    return - A / (2*α*β) * cos(k*x) * exp(-(y^2)/4)

end

function calcvH(x,y,β,A,α,k)

    if y != 0
        return A / (β*y) * (k/α * sin(k*x) - 0.5/β*cos(k*x)) * exp(-(y^2)/4)
    else
        return 0
    end

end

function calcϕH(x,y,A,α,k)

    return - A / α * cos(k*x) * exp(-(y^2)/4)

end

function GMsmallH(xvec,yvec,β,A,α,L)

    k = pi / (2*L); nx = length(xvec); ny = length(yvec);
    u = zeros(ny,nx); v = zeros(ny,nx); ϕ = zeros(ny,nx);

    for ii = 1 : nx
        for jj = 1 : ny
            if abs(xvec[ii]) <= L
                u[jj,ii] = calcuH(xvec[ii],yvec[jj],β,A,α,k)
                v[jj,ii] = calcvH(xvec[ii],yvec[jj],β,A,α,k)
                ϕ[jj,ii] = calcϕH(xvec[ii],yvec[jj],A,α,k)
            else
                u[jj,ii] = 0
                v[jj,ii] = 0
                ϕ[jj,ii] = 0
            end
        end
    end

    return u,v,ϕ

end
