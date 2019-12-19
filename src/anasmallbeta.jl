function calcuβ(x,y,A,α,k,K,β)

    return K/α * exp(-abs(x)-β*(y^2)*0.5/α) * sign(x)

end

function calcvβ(x,y,A,α,k,K,β,c,L)

    v = K * sqrt(pi*α*0.5*β) * (1/α - α/(c^2)) * exp(-x) * erf(y*sqrt(β*0.5*α)) + α*K*y/(c^2)
    if abs(x) > L
        v =+ A*sqrt(pi)/(c^2) * cos(k*x) * erf(y/2)
    end
    return v

end

function calcϕβ(x,y,A,α,k,K,β)

    return K * exp(-abs(x)-β*(y^2)*0.5/α)

end

function GMsmallβ(xvec,yvec,β,A,α,g,H,L)

    c = sqrt(g*H); k = pi / (2*L); K = -1;
    nx = length(xvec); ny = length(yvec);
    u = zeros(ny,nx); v = zeros(ny,nx); ϕ = zeros(ny,nx);

    for ii = 1 : nx
        for jj = 1 : ny
            u[jj,ii] = calcuβ(xvec[ii],yvec[jj],A,α,k,K,β)
            v[jj,ii] = calcvβ(xvec[ii],yvec[jj],A,α,k,K,β,c,L)
            ϕ[jj,ii] = calcϕβ(xvec[ii],yvec[jj],A,α,k,K,β)
        end
    end

    return u,v,ϕ

end
