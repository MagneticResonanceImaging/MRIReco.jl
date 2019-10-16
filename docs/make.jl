using Documenter, MRIReco

makedocs(
    modules = [MRIReco],
    format = Documenter.HTML(prettyurls = false),
    sitename = "Julia MRI Package",
    authors = "Tobias Knopp, Mirco Grosser",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "gettingStarted.md",
        "Overview" => "overview.md",
        "Acquisition Data" => "acquisitionData.md",
        "File Handling" => "filehandling.md",
        "Images" => "image.md",
        "Trajectory" => "trajectories.md",
        "Imaging Operators" => "operators.md",
        "Offresonance" => "offresonance.md",
        "Parallel Imaging" => "SENSE.md",
        "Compressed Sensing" => "compressedSensing.md",
        "Customize" => "custom.md",
        "API" => "API.md",
    ],
)

deploydocs(
    repo = "github.com/MagneticResonanceImaging/MRIReco.jl.git",
    target = "build",
)
