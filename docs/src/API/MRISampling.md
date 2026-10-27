# MRISampling

MRISampling creates undersampling patterns and applies them to k-space data. MRISampling depends on MRIBase. It has to be loaded with `using MRISampling`, because MRIReco does not depend on it.

## Undersampling

```@autodocs
Modules = [MRISampling]
Pages = ["src/MRISampling.jl"]
Private = false
```

## Sampling Patterns

```@autodocs
Modules = [MRISampling]
Pages = ["src/Regular.jl", "src/Random.jl", "src/Lines.jl", "src/PoissonDisk.jl", "src/VDPoissonDisk.jl"]
Private = false
```

## Incoherence

```@autodocs
Modules = [MRISampling]
Pages = ["src/Incoherence.jl"]
Private = false
```
