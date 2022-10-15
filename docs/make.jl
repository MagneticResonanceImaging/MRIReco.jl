using Documenter, MRIBase, MRIFiles, MRICoilSensitivities, 
                  MRIOperators, MRISampling, MRISimulation, MRIReco

makedocs(
    modules = [MRIReco],
    sitename = "Julia MRI Package",
    authors = "Tobias Knopp, Mirco Grosser, and Co-Authors",
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
        "API" => "API.md"
    ],
)

deploydocs(
    repo = "github.com/MagneticResonanceImaging/MRIReco.jl.git"
)
