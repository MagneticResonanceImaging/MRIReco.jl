using Documenter, MRIReco

makedocs(
    modules = [MRIReco],
    format = :html,
    sitename = "MRIReco.jl",
    pages = Any["index.md"]
)

deploydocs(
    repo = "github.com/tknopp/MRIReco.jl.git",
    target = "build",
    julia = "1.0",
    deps = nothing,
    make = nothing,
)
