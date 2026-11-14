# MRIReco

MRIReco provides the image reconstruction of the framework. Of the packages in this repository, it depends on MRIBase and MRIOperators. It re-exports both, together with [RegularizedLeastSquares.jl](https://juliaimagerecon.github.io/RegularizedLeastSquares.jl/stable/), so the data types, encoding operators, solvers and regularization terms are available after `using MRIReco`. MRIFiles, MRISampling, MRICoilSensitivities and MRISimulation are independent of MRIReco and have to be loaded separately.

## Reconstruction

```@docs
MRIReco.reconstruction
```

## Lower-Level Reconstruction Functions

For the reconstruction types `"standard"`, `"multiEcho"`, `"multiCoil"` and `"multiCoilMultiEcho"`, `reconstruction` calls one of the following functions.

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
