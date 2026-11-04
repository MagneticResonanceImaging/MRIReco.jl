# MRIOperators

MRIOperators implements the encoding operators used in MRI reconstruction: Fourier encoding with and without field inhomogeneity correction, coil sensitivities and sparsifying transforms. It builds on [LinearOperatorCollection.jl](https://juliaimagerecon.github.io/LinearOperatorCollection.jl/stable/) and re-exports its operators, such as `FFTOp`, `NFFTOp`, `WaveletOp`, `DiagOp` and `normalOperator`, which are documented there. MRIOperators depends on MRIBase. MRIReco re-exports MRIOperators, so `using MRIReco` makes all names on this page available.

## Encoding Operators

```@autodocs
Modules = [MRIOperators]
Pages = ["src/EncodingOp.jl"]
Private = false
```

## Operator Types

```@autodocs
Modules = [MRIOperators]
Pages = ["src/ExplicitOp.jl", "src/FieldmapNFFTOp.jl", "src/SensitivityOp.jl", "src/SubspaceOp.jl", "src/MapSliceOp.jl"]
Private = false
```

## Sparsifying Transforms

```@autodocs
Modules = [MRIOperators]
Pages = ["src/SparseOp.jl"]
Private = false
```
