# MRIFiles

MRIFiles handles raw data files, for example ISMRMRD files (`ISMRMRDFile`) and Bruker datasets (`BrukerFile`). MRIFiles depends on MRIBase. It has to be loaded with `using MRIFiles`, because MRIReco does not depend on it. The guide [Load and Save Raw Data](@ref load-save) shows how to use these file types.

## ISMRMRD Files

```@docs
ISMRMRDFile
MRIBase.RawAcquisitionData(::ISMRMRDFile, ::Any)
MRIBase.AcquisitionData(::ISMRMRDFile, ::Any)
MRIFiles.save(::ISMRMRDFile, ::RawAcquisitionData)
```

## Bruker Datasets

```@docs
BrukerFile
MRIBase.RawAcquisitionData(::BrukerFile)
RawAcquisitionDataFid_360
RawAcquisitionData_3DUTE
recoData
JcampdxFile
```

## Other Formats

```@docs
CFLFile
DFFile
readRawSiemens
```
