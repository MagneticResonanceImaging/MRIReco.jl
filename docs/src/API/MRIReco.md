# MRIReco

MRIReco provides the image reconstruction of the framework. Of the packages in this repository, it depends on MRIBase and MRIOperators. It re-exports both, together with [RegularizedLeastSquares.jl](https://juliaimagerecon.github.io/RegularizedLeastSquares.jl/stable/), so the data types, encoding operators, solvers and regularization terms are available after `using MRIReco`. MRIFiles, MRISampling, MRICoilSensitivities and MRISimulation are independent of MRIReco and have to be loaded separately.

## Reconstruction

```@docs
MRIReco.reconstruction
```

## Reconstruction Plans

`reconstruction` translates the parameter dictionary into a reconstruction plan of
[AbstractImageReconstruction.jl](https://github.com/JuliaImageRecon/AbstractImageReconstruction.jl) and runs it.
A plan for each reconstruction type is stored in `reco_plans/` and can be loaded with `MRIRecoPlan`, e.g.
`MRIRecoPlan("multiCoil")`, or run directly with `reconstruct("multiCoil", acqData; kwargs...)`.

```@autodocs
Modules = [MRIReco]
Pages = ["Algorithms/Algorithms.jl", "Algorithms/Storage.jl", "Algorithms/MRIRecoContext.jl", "Algorithms/Direct.jl"]
Private = false
```

### Iterative Reconstructions

```@autodocs
Modules = [MRIReco]
Pages = ["Iterative/Iterative.jl", "Iterative/IterativeContextParameters.jl", "Iterative/Standard.jl",
         "Iterative/MultiEcho.jl", "Iterative/MultiCoil.jl", "Iterative/MultiCoilMultiEcho.jl",
         "Iterative/MultiCoilMultiEchoSubspace.jl"]
Private = false
```

### Parameters

```@autodocs
Modules = [MRIReco]
Pages = ["Parameters/EncodingParameters.jl", "Parameters/WeightingParameters.jl", "Parameters/CoilParameters.jl",
         "Parameters/SolverParameters.jl"]
Private = false
```

## Deprecated Reconstruction Functions

Before reconstruction plans, `reconstruction` called one of the following functions. They are deprecated.

```@docs
MRIReco.reconstruction_simple
MRIReco.reconstruction_multiEcho
MRIReco.reconstruction_multiCoil
MRIReco.reconstruction_multiCoilMultiEcho
```

## Tools

```@autodocs
Modules = [MRIReco]
Pages = ["Tools/PreProcess.jl", "Tools/Regridding.jl", "Tools/Subspace.jl", "Tools/ErrorMeasures.jl"]
Private = false
```
