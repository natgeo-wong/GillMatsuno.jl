module GillMatsuno

# Just a small module I made to write various forms of the shallow water equations that are relevant
# to the atmosphere

## Modules Used
using Dates
using NCDatasets
using Printf
using TimerOutputs
import Base: eltype, show

export
        GenerateGrid,
        DomainProperties,
        QfieldProperties,
        CreateSimulation,
        runGillMatsuno,
        createQ

## Including other files in the module

include("Grid.jl")
include("Domain.jl")
include("QForcing.jl")
include("Simulation.jl")
include("Run.jl")

end # module
