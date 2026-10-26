using Documenter, Literate
using MRIBase, MRIFiles, MRICoilSensitivities, MRIOperators, MRISampling, MRISimulation, MRIReco

# Generates examples
include("generate.jl")

makedocs(
    modules = [MRIReco, MRIBase, MRIOperators, MRIFiles, MRISampling, MRICoilSensitivities, MRISimulation],
    sitename = "Julia MRI Package",
    authors = "Tobias Knopp, Mirco Grosser, and Co-Authors",
    repo = Remotes.GitHub("MagneticResonanceImaging", "MRIReco.jl"),
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://magneticresonanceimaging.github.io/MRIReco.jl/stable/",
        assets = String[],
        collapselevel = 1,
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "gettingStarted.md",
        "Tutorials" => [
            "Subspace Reconstruction" => "generated/tutorials/subspace.md",
        ],
        "How-to Guides" => [
            "Raw Data and Files" => [
                "Bin Profiles into Contrasts" => "generated/howto/binning.md",
                "Extract Cartesian k-Space" => "generated/howto/extract_kspace.md",
            ],
        ],
        "Explanations" => [
            "Package Architecture" => "explanations/architecture.md",
            "Data Model" => "generated/explanations/data_model.md",
        ],
        "Reference" => [
            "API" => [
                "MRIReco" => "API/MRIReco.md",
                "MRIBase" => "API/MRIBase.md",
                "MRIOperators" => "API/MRIOperators.md",
                "MRIFiles" => "API/MRIFiles.md",
                "MRISampling" => "API/MRISampling.md",
                "MRICoilSensitivities" => "API/MRICoilSensitivities.md",
                "MRISimulation" => "API/MRISimulation.md",
            ],
        ],
        "Legacy Pages" => [
            "Acquisition Data" => "acquisitionData.md",
            "File Handling" => "filehandling.md",
            "Trajectory" => "trajectories.md",
            "Imaging Operators" => "operators.md",
            "Offresonance" => "offresonance.md",
            "Parallel Imaging" => "SENSE.md",
            "Compressed Sensing" => "compressedSensing.md",
            "Customize" => "custom.md",
        ],
    ],
    pagesonly = true,
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/MagneticResonanceImaging/MRIReco.jl.git",
    devbranch = "master",
)
