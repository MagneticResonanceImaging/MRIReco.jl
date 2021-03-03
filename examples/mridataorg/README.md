# Example of mridata.org reconstruction

This folder contains example code for downloading data from the website mridata.org 
and performing image reconstruction of the undersampled FSE data. The data is stored
in the ISMRMRD file format

## Installation

In order to use this code one first has to download Julia (version 1.5) and install
`MRIReco` by executing

```julia
using Pkg
Pkg.add("MRIReco")
```

Load the package by entering
```julia
using MRIReco
```

You can then switch to the directory of this example by entering
```julia
dir = joinpath(dirname(pathof(MRIReco)), "..","examples","mridataorg")
cd(dir)
```

## Execution
After installation and switching to the example directory, the example code can be
executed by entering

```julia
include("example.jl")
```

This will first download the ISMRMRD data (about 167 MB) and then perform a reconstruction.
Parameters of the reconstruction are documented in the Julia script and can be
changed. After the reconstruction is done, the script opens a PyPlot window
and show the reconstruction result. 
