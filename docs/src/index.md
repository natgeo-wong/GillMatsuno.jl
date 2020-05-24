# GillMatsuno.jl
*Numerical Solutions to the Gill-Matsuno Equations*

`GillMatsuno.jl` is a Julia package that:
* numerically solves the Shallow-Water Equations on a $\beta$-plane
* allows the user to define custom heat-forcing $Q$

## Installation
`GillMatsuno.jl` can be installed using Julia's built-in package manager as follows:

```
julia> ]
(@v1.4) pkg> add GillMatsuno
```

You can update `GillMatsuno.jl` to the latest version using
```
(@v1.4) pkg> update GillMatsuno
```

And if you want to get the latest release without waiting for me to update the Julia Registry (although this generally isn't necessary since I make a point to release patch versions as soon as I find bugs or add new working features), you may fix the version to the `master` branch of the GitHub repository:
```
(@v1.4) pkg> add GillMatsuno#master
```

## Documentation

The documentation for `GillMatsuno.jl` is covers:
1. Theory - the equations behind the Gill-Matsuno Equations
2. Examples - the classical Gill-Matsuno equation with default values
3. API - comprehensive summary of all exported functionalities

## Getting help
If you are interested in using `GillMatsuno.jl` or are trying to figure out how to use it, please feel free to ask me questions and get in touch!  Please feel free to [open an issue](https://github.com/natgeo-wong/GillMatsuno.jl/issues/new) if you have any questions, comments, suggestions, etc!
