# MRISimulation

MRISimulation simulates MRI data from images, including coil sensitivities, field inhomogeneities and relaxation. MRIReco does not re-export MRISimulation, so it has to be loaded with `using MRISimulation`.

## Simulation

```@autodocs
Modules = [MRISimulation]
Pages = ["src/MRISimulation.jl"]
Private = false
```

## Coil Sensitivities, Field Maps and Relaxation

```@autodocs
Modules = [MRISimulation]
Pages = ["src/CoilSensitivities.jl", "src/Fieldmap.jl", "src/RelaxationMap.jl"]
Private = false
```

## Sequences

```@autodocs
Modules = [MRISimulation]
Pages = ["Sequences/MultiEcho.jl", "Sequences/BlochSimulation.jl"]
Private = false
```
