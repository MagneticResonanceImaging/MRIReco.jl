using Documenter, MRIReco

makedocs(
    modules = [MRIReco],
    format = :html,
    sitename = "MRIReco.jl",
    authors = "Tobias Knopp, Mirco Grosser",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "gettingStarted.md",
    ],
)

deploydocs(
    repo = "github.com/MagneticResonanceImaging/MRIReco.jl.git",
    target = "build",
    julia = "1.0",
    deps = nothing,
    make = nothing,
)
