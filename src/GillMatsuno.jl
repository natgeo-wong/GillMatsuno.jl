module GillMatsuno

# Just a small module I made to write various forms of the shallow water equations that are relevant
# to the atmosphere

## Modules Used
using Seaborn, NetCDF

export
        GManalytic, GMcalc

## Including other files in the module

include("anadefault.jl")
include("anasmallbeta.jl")
include("anasmallH.jl")

include("numdefault.jl")
include("numwaves.jl")


end # module
