# GillMatsuno

This package contains codes for both the analytic and numerical time-stepping solutions to
the Gill-Matsuno equations, as described by Vallis [2018] (Atmospheric and Oceanic Fluid
Dynamics, 2nd Edition):

$$
u\_t + \\alpha u + \\phi_x - βyv   &= 0
v\_t + \\alpha u + \\phi_y + βyu   &= 0 \\
ϕ\_t + \\alpha\\phi + gH(u\_x+v\_y) &= -Q
$$
