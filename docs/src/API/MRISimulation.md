# MRISimulation

MRISimulation simulates MRI data from images, including coil sensitivities, field inhomogeneities and relaxation. MRIReco does not re-export MRISimulation, so it has to be loaded with `using MRISimulation`.

## Simulation

```@docs
MRISimulation.simulation
MRISimulation.addNoise
MRISimulation.simulation_fast
MRISimulation.simulation_explicit
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
