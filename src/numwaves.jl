"""
The Linear Shallow Water Equations in 1 Dimension are given by:
u_t + g * η_x = 0
η_t + H * u_x = 0

which can then be reduced to:
η_tt - c^2 * η_xx = 0, where c^2 = g * H

which is simply the wave equation.  The discretized for is therefore given by
η_j^(t+1) = 2(1-r^2)η_j^t - η_j^(t-1) + r^2(η_(j+1)^t + η_(j+1)^t),
where r = cΔt/Δx
"""

function shallowsetup1D(height,width,domain,depth,gravity=9.81)
end

function shallowwave1D(ηn,ηm1,r,boundarycondition::AbstractString)

    ηp1 = zeros(size(ηn))

    if boundarycondition == "periodic"
        ηp1 = 2 * (1-r^2) * ηn - ηm1 + r^2 * (circshift(ηn,1) + circshift(ηn,-1))
    elseif boundarycondition == "reflective"
        ηp1[2:end-1] = 2 * (1-r^2) * ηn[2:end-1] - ηm1[2:end-1] + r^2 * (ηn[1:end-2] + ηn[3:end]);
        ηp1[1] = 2 * (1-r^2) * ηn[1] - ηm1[1] + r^2 * (ηn[1] + ηn[2]);
        ηp1[end] = 2 * (1-r^2) * ηn[end] - ηm1[end] + r^2 * (ηn[end] + ηn[end-1]);
    end

    return ηp1

end

function shallowwave2D(ηn,ηm1,r,boundarycondition::AbstractString)

    ηp1 = zeros(size(ηn));

    if boundarycondition == "periodic"
        ηp1 = 2 * (1-2*r^2) * ηn - ηm1 + r^2 * (circshift(ηn,(1,0)) + circshift(ηn,(-1,0)) +
        + circshift(ηn,(0,1)) + circshift(ηn,(0,-1)));
    elseif boundarycondition == "reflective"
        ηp1[2:end-1,2:end-1]  = 2 * (1-2*r^2) * ηn[2:end-1,2:end-1] - ηm1[2:end-1,2:end-1];
        ηp1[2:end-1,2:end-1] += r^2 * (ηn[1:end-2,2:end-1] + ηn[3:end,2:end-1] +
        + ηn[2:end-1,3:end] + ηn[2:end-1,1:end-2]);

        ηp1[1,2:end-1] = 2 * (1-2*r^2) * ηn[1,2:end-1] - ηm1[1,2:end-1] + r^2 * (ηn[1,2:end-1] +
        + ηn[2,2:end-1] + ηn[1,3:end] + ηn[1,1:end-2]);
        ηp1[end,2:end-1] = 2 * (1-2*r^2) * ηn[end,2:end-1] - ηm1[end,2:end-1] +
        + r^2 * (ηn[end,2:end-1] + ηn[end-1,2:end-1] + ηn[end,3:end] + ηn[end,1:end-2]);
        ηp1[2:end-1,1] = 2 * (1-2*r^2) * ηn[2:end-1,1] - ηm1[2:end-1,1] +
        + r^2 * (ηn[2:end-1,1] + ηn[2:end-1,2] + ηn[3:end,1] + ηn[1:end-2,1]);
        ηp1[2:end-1,end] = 2 * (1-2*r^2) * ηn[2:end-1,end] - ηm1[2:end-1,end] +
        + r^2 * (ηn[2:end-1,end] + ηn[2:end-1,end-1] + ηn[3:end,end] + ηn[1:end-2,end]);

        np1[1,1] = 2 * (1-2*r^2) * ηn[1,1] - ηm1[1,1] + r^2 * (2*ηn[1,1] + ηn[1,2] + ηn[2,1]);
        np1[1,end] = 2 * (1-2*r^2) * ηn[1,end] - ηm1[1,end] +
        + r^2 * (2*ηn[1,end] + ηn[1,end-1] + ηn[2,end]);
        np1[end,1] = 2 * (1-2*r^2) * ηn[end,1] - ηm1[end,1] +
        + r^2 * (2*ηn[end,1] + ηn[end,2] + ηn[end-1,1]);
        np1[end,end] = 2 * (1-2*r^2) * ηn[end,end] - ηm1[end,end] +
        + r^2 * (2*ηn[end,end] + ηn[end,end-1] + ηn[end-1,end]);
    end

    return ηp1

end
