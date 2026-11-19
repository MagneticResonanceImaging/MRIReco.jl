using Documenter, DocumenterVitepress, Literate
using MRIBase, MRIFiles, MRICoilSensitivities, MRIOperators, MRISampling, MRISimulation, MRIReco

# Generates examples
include("generate.jl")

# The contributing guide is maintained in CONTRIBUTING.md. Copy it into the docs and point
# its links to other files of the repository to GitHub.
contributing = replace(read(joinpath(@__DIR__, "..", "CONTRIBUTING.md"), String), "\r\n" => "\n")
contributing = replace(contributing, r"\]\((?!https?://|#)([^)]+)\)" =>
    s"](https://github.com/MagneticResonanceImaging/MRIReco.jl/blob/master/\1)")
write(joinpath(@__DIR__, "src", "contributing.md"),
      "```@meta\nEditURL = \"../../CONTRIBUTING.md\"\n```\n\n" * contributing)

makedocs(
    modules = [MRIReco, MRIBase, MRIOperators, MRIFiles, MRISampling, MRICoilSensitivities, MRISimulation],
    sitename = "MRIReco.jl",
    authors = "Tobias Knopp, Mirco Grosser, and Co-Authors",
    repo = Remotes.GitHub("MagneticResonanceImaging", "MRIReco.jl"),
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/MagneticResonanceImaging/MRIReco.jl",
        devbranch = "master",
        devurl = "dev",
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "gettingStarted.md",
        "Tutorials" => [
            "Direct Reconstruction" => "generated/tutorials/direct.md",
            "Iterative Reconstruction" => "generated/tutorials/iterative.md",
            "Multi-Coil Reconstruction" => "generated/tutorials/multicoil.md",
            "Multi-Echo Reconstruction" => "generated/tutorials/multiecho.md",
            "Multi-Coil Multi-Echo and Subspace Reconstruction" => "generated/tutorials/subspace.md",
        ],
        "How-to Guides" => [
            "Raw Data and Files" => [
                "Load and Save Raw Data" => "generated/howto/load_save.md",
                "Select and Prepare Profiles" => "generated/howto/prepare_profiles.md",
                "Bin Profiles into Contrasts" => "generated/howto/binning.md",
                "Extract Cartesian k-Space" => "generated/howto/extract_kspace.md",
                "Create a Trajectory" => "generated/howto/trajectory.md",
                "Create AcquisitionData from Arrays" => "generated/howto/from_arrays.md",
            ],
            "Sampling" => [
                "Undersample k-Space" => "generated/howto/undersample.md",
            ],
            "Coil Sensitivities" => [
                "Estimate Coil Sensitivities with ESPIRiT" => "generated/howto/espirit.md",
                "Compress Multi-Coil Data" => "generated/howto/coil_compression.md",
            ],
            "Simulation" => [
                "Simulate k-Space Data" => "generated/howto/simulation.md",
            ],
            "Encoding Operators" => [
                "Use Encoding Operators Directly" => "generated/howto/encoding_operators.md",
            ],
            "Reconstruction" => [
                "Configure a Reconstruction with Plans" => "generated/howto/plans.md",
                "Correct Off-Resonance" => "generated/howto/offresonance.md",
                "Run a Reconstruction on the GPU" => "generated/howto/gpu.md",
                "Use a Custom Sparsifying Transform" => "generated/howto/custom_transform.md",
                "Use a Plug-and-Play Prior" => "generated/howto/plug_and_play.md",
            ],
            "Contribute to MRIReco.jl" => "contributing.md",
        ],
        "Explanations" => [
            "Package Architecture" => "explanations/architecture.md",
            "Data Model" => "generated/explanations/data_model.md",
            "Signal Model and Encoding Operators" => "explanations/signal_model.md",
            "Reconstruction as an Inverse Problem" => "explanations/inverse_problem.md",
            "GPU Acceleration" => "explanations/gpu.md",
        ],
        "Reference" => [
            "Reconstruction Parameters" => "reference/parameters.md",
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
    ],
    pagesonly = true,
    checkdocs = :exports,
)

# VitePress copies only the pages and the files they reference into the site, so copy the notebooks that
# generate.jl writes next to the Literate pages into each built site
generated = joinpath(@__DIR__, "src", "generated")
for site in readdir(joinpath(@__DIR__, "build"); join = true)
    (isdir(site) && all(isdigit, basename(site))) || continue
    for (root, _, files) in walkdir(generated), file in filter(endswith(".ipynb"), files)
        target = joinpath(site, "generated", relpath(root, generated))
        mkpath(target)
        cp(joinpath(root, file), joinpath(target, file); force = true)
    end
end

# Pull requests are not deployed. For those from this repository,
# DocumenterVitepress.deploydocs would still post a status that links to a preview.
if get(ENV, "GITHUB_EVENT_NAME", "") != "pull_request"
    DocumenterVitepress.deploydocs(
        repo = "github.com/MagneticResonanceImaging/MRIReco.jl.git",
        target = joinpath(@__DIR__, "build"),
        devbranch = "master",
    )
end
