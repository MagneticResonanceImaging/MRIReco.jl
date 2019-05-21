import ImageMagick
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
        "File Handling" => "filehandling.md",
        "Images" => "image.md",
        "Offresonance" => "offresonance.md",
        "Parallel Imaging" => "SENSE.md",
        "Trajectory" => "trajectories.md",
        "Imaging Operators" => "operators.md",
        "Simulation" => "simulation.md",
        "Reconstruction" => "reconstruction.md",
        "API" => "API.md",
    ],
    html_prettyurls = false, #!("local" in ARGS),
)

deploydocs(
    repo = "github.com/MagneticResonanceImaging/MRIReco.jl.git",
    target = "build",
)
