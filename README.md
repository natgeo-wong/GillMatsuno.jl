# **<div align="center">GillMatsuno.jl</div>**

**Created By:** Nathanael Wong (nathanaelwong@fas.harvard.edu)

## **Introduction**

`GillMatsuno.jl` is a Julia package that:
* numerically solves the Shallow-Water Equations on a $\beta$-plane
* allows the user to define custom heat-forcing $Q$

`GillMatsuno.jl` can be installed via
```
] add GillMatsuno
```

## Equations
The Shallow-Water Equations are as described by *Vallis [2018]* (Atmospheric and Oceanic Fluid Dynamics, 2nd Edition):
```
∂u/∂t + αu + ∂ϕ/∂x - βyv = 0
∂v/∂t + αu + ∂ϕ/∂y + βyu = 0
∂ϕ/∂t + αϕ + gH(∂u/∂x + ∂v/∂y) = -Q
```

We solve the equations numerically using the Finite-Difference Method.  The backend solvers can be found in `src/numdefault.jl`, but this README will not go through the details.

## Running the Gill-Matsuno Numerical Model
The numerical solution to the Gill-Matusno Model is called using `GMcalc()`.  This function accepts the following keyword arguments as modifications to the initial model:
* `xmin` and `xmax` define the left and right boundaries of the domain
* `ymin` and `ymax` define the lower and upper boundaries of the domain
* `δx` and `δy` define the grid-spacing in the x- and y- directions respectively
* `nt` and `δt` define the number of timesteps, and the time between each step, respectively
* `A` and `L` represent the amplitude and length of the heat-forcing `Q`
* `α` represents the damping coefficient on the winds induced by the heat-forcing (default: `α = 0.1`)
* `β` is the Coriolis Factor (nondimensionalized to `β = 0.5` as the default)
* `g` and `H` represent gravity and the height of the domain (both nondimensionalized to 1 as default)

An example is given below, with the default values as follows:
```
GMcalc(xmin=-25.0,xmax=50.0,δx=0.1,ymin=-10.0,ymax=10.0,δy=0.1,nt=5000,δt=0.001,
	   A=1,L=2,α=0.1,β=0.5,g=1.0,H=1.0)
```
