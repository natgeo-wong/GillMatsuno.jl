# **<div align="center">GillMatsuno.jl</div>**

<p align="center">
  <a href="https://www.repostatus.org/#active">
    <img alt="Repo Status" src="https://www.repostatus.org/badges/latest/active.svg?style=flat-square" />
  </a>
  <a href="https://travis-ci.com/github/natgeo-wong/GillMatsuno.jl">
    <img alt="Travis CI" src="https://travis-ci.com/natgeo-wong/GillMatsuno.jl.svg?branch=master&style=flat-square">
  </a>
  <a href="https://github.com/natgeo-wong/GillMatsuno.jl/actions?query=workflow%3ADocumentation">
    <img alt="Documentation Build" src="https://github.com/natgeo-wong/GillMatsuno.jl/workflows/Documentation/badge.svg">
  </a>
  <br>
  <a href="https://mit-license.org">
    <img alt="MIT License" src="https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square">
  </a>
  <img alt="Latest Release" src="https://img.shields.io/github/v/release/natgeo-wong/GillMatsuno.jl">
  <a href="https://natgeo-wong.github.io/GillMatsuno.jl/stable/">
    <img alt="Latest Documentation" src="https://img.shields.io/badge/docs-stable-blue.svg?style=flat-square">
  </a>
  <a href="https://natgeo-wong.github.io/GillMatsuno.jl/dev/">
    <img alt="Latest Documentation" src="https://img.shields.io/badge/docs-latest-blue.svg?style=flat-square">
  </a>
</p>

**Created By:** Nathanael Wong (nathanaelwong@fas.harvard.edu)

## Introduction

`GillMatsuno.jl` is a Julia package that:
* numerically solves the Shallow-Water Equations on a $\beta$-plane using finite-difference methods
* allows the user to define custom heat-forcing $Q$

`GillMatsuno.jl` can be installed via
```
] add GillMatsuno
```

Due to the recent improvements in memory allocations in Julia, `GillMatsuno.jl` `v2` works best in Julia `v1.5` and above, but can work from `v1.3` onwards.

## Using `GillMatsuno.jl`
There are four components to running a model in `GillMatsuno.jl`.  They are:
1. Grid `G`
2. Domain parameters `D`
3. Heat Source `Q`
4. Simulation Setup `S`

### 1. Defining the Grid `G`

`GillMatsuno.jl` uses a staggered Arakawa C-Grid.  The Grid `G` is generated via the function `GenerateGrid`, as follows
```
G = GenerateGrid(size = (nx,ny), x = (xmin,xmax), y = (ymin,ymax))
```

It is to be noted that the shallow-water equations to be solved have been nondimensionalized.  Typical values of `xmin` and `xmax` are *O*(25) (negative and positive respectively), and *O*(10) for `ymin` and `ymax`.

### 2. Defining the Domain Parameters `D`

We define the domain parameters using the `DomainProperties()` function.  The default values are:
* `α` represents the damping coefficient on the winds induced by the heat-forcing (default: `α = 0.1`)
* `β` is the Coriolis Factor (nondimensionalized to `β = 0.5` as the default)
* `g` and `H` represent gravity and the height of the domain (both nondimensionalized to 1 as default)
```
D = DomainParameters(α=0.2,β=0.5,g=1.0,H=1.0)
```

### 3. Defining the Heat Source *Q*

The heat source *Q* is analogous to a mass source/sink.  As of now, *Q* can only be defined as a gaussian peak (or the cumulative sums of gaussian peaks), though we aim to extend this to equatorial bands.

*Q* can be defined via the function `QfieldProperties`
```
Q = QfieldProperties(A=1.0,Lx=2.0,Ly=2.0,Qx=0,Qy=0)
```
Where we have that
* `A` is the amplitude of the source
* `Lx` and `Ly` are the non-dimensionalized widths of the `Q` in the `x`- and `y`-directions respectively
* `Qx` and `Qy` denote the location of the center of `Q`

### 4. Setting up the Simulation `S`

The simulation structure `S` is defined as follows:
```
S = CreateSimulation(δt=5e-4,tt=50,ft=0.5,fnc="test.nc")
```
Where we have that
* `δt` is the model timestep
* `tt` is the total model runtime
* `ft` is the output frequency in model runtime

So, using the parameters above, we see that the model is ran for 10^6 timesteps, with the fields output every 10000 steps to the netCDF file `test.nc`.

### 5. Running the Simulation

With the fields we have defined above, we put them into the function `runGillMatsuno(S,G,[Q],D)`, and then we can extract the fields and do plotting/analysis, as you wish!
