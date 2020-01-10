# GillMatsuno

This package contains codes for both the analytic and numerical time-stepping solutions to
the Gill-Matsuno equations, as described by Vallis [2018] (Atmospheric and Oceanic Fluid
Dynamics, 2nd Edition):

```
∂u/∂t + αu + ∂ϕ/∂x - βyv = 0
∂v/∂t + αu + ∂ϕ/∂y + βyu = 0
∂ϕ/∂t + αϕ + gH(∂u/∂x + ∂v/∂y) = -Q
```
