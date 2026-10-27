# MRISimulation

MRISimulation simulates MRI data from images, including coil sensitivities, field inhomogeneities and relaxation. MRISimulation depends on MRIBase and MRIOperators. It has to be loaded with `using MRISimulation`, because MRIReco does not depend on it.

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
