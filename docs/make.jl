using Documenter, MRIReco

makedocs(
    modules = [MRIReco],
    format = :html,
    sitename = "Julia MRI Package",
    authors = "Tobias Knopp, Mirco Grosser",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "gettingStarted.md",
        "Acquisition Data" => "acquisitionData.md",
        "Images" => "image.md",
        "Offresonance" => "offresonance.md",
        "Parallel Imaging" => "SENSE.md",
        "Simulation" => "simulation.md",
        "Reconstruction" => "reconstruction.md",
        "Trajectory" => "trajectories.md",
        "Operators" => "operators.md",
    ],
    html_prettyurls = false, #!("local" in ARGS),
)

deploydocs(
    repo = "github.com/MagneticResonanceImaging/MRIReco.jl.git",
    target = "build",
)
