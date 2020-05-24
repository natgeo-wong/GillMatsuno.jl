using Documenter
using GillMatsuno

makedocs(
    modules  = [GillMatsuno],
    doctest  = false,
    format   = Documenter.HTML(
        collapselevel = 1,
        prettyurls    = false
    ),
    authors  = "Nathanael Wong",
    sitename = "GillMatsuno.jl",
    pages    = [
        "Home"     => "index.md",
        # "Theory"   => "theory.md",
        # "Examples" => "examples.md"
    ]
)

deploydocs(
    repo = "github.com/natgeo-wong/GillMatsuno.jl.git",
)
