# MRICoilSensitivities

MRICoilSensitivities estimates coil sensitivity maps, for example with ESPIRiT, and compresses multi-channel data to fewer virtual coils. MRICoilSensitivities depends on MRIBase. It has to be loaded with `using MRICoilSensitivities`, because MRIReco does not depend on it.

## Coil Sensitivities

```@autodocs
Modules = [MRICoilSensitivities]
Pages = ["src/espirit.jl", "src/simple.jl"]
Private = false
```

## Coil Compression

```@autodocs
Modules = [MRICoilSensitivities]
Pages = ["src/CoilCompression.jl"]
Private = false
```
