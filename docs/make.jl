using Documenter, MRIReco

makedocs(
    modules = [MRIReco],
    format = :html,
    sitename = "Julia MRI Package",
    authors = "Tobias Knopp, Mirco Grosser",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "gettingStarted.md",
        "Trajectory" => "trajectories.md",
        "Acquisition Data" => "acquisitionData.md",
        "Operators" => "operators.md",
        "Simulation" => "simulation.md",
        "Reconstruction" => "reconstruction.md",
    ],
)

deploydocs(
    repo = "github.com/MagneticResonanceImaging/MRIReco.jl.git",
    target = "build",
    julia = "1.0",
    deps = nothing,
    make = nothing,
)
