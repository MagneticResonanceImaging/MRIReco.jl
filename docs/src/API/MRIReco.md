# MRIReco

MRIReco provides the image reconstruction of the framework. `using MRIReco` also loads MRIBase, MRIOperators and [RegularizedLeastSquares.jl](https://juliaimagerecon.github.io/RegularizedLeastSquares.jl/stable/) and re-exports their names, so the data types, encoding operators, solvers and regularization terms are available without further imports.

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
