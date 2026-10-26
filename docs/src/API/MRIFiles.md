# MRIFiles

MRIFiles handles raw data files, for example ISMRMRD files (`ISMRMRDFile`) and Bruker datasets (`BrukerFile`). MRIFiles depends on MRIBase. MRIReco does not depend on MRIFiles, so it has to be loaded with `using MRIFiles`.

## ISMRMRD Files

```@autodocs
Modules = [MRIFiles]
Pages = ["ISMRMRD/ISMRMRD.jl"]
```
