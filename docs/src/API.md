# API

## Operators
Operators are implemented as subtypes of `AbstractLinearOperator`, which is defined in the package `LinearOperators.jl`. Such operators must provide a function implementing the product and a function implementing the product with the adjoint. Furthermore, the number of rows and columns of the operator must be specified.


```@docs
MRIReco.encodingOps2d_simple
MRIReco.encodingOps3d_simple
MRIReco.encodingOps2d_parallel
MRIReco.encodingOps3d_parallel
MRIReco.encodingOp2d_multiEcho
MRIReco.encodingOp3d_multiEcho
MRIReco.encodingOp2d_multiEcho_parallel
MRIReco.encodingOp3d_multiEcho_parallel
MRIReco.fourierEncodingOp2d
MRIReco.fourierEncodingOp3d
MRIReco.ExplicitOp(shape::NTuple{D,Int64}, tr::Trajectory, correctionmap::Array{ComplexF64,D}; MRIReco.echoImage::Bool=false, kargs...) where D
RegularizedLeastSquares.FFTOp(T::Type, shape::Tuple, shift=true)
MRIReco.NFFTOp(shape::Tuple, tr::Trajectory; nodes=nothing, kargs...)
MRIReco.FieldmapNFFTOp(shape::NTuple{D,Int64}, tr::Trajectory,
                        correctionmap::Array{ComplexF64,D};
                        method::String="nfft",
                        echoImage::Bool=true,
                        alpha::Float64=1.75,
                        m::Float64=3.0,
                        K=20,
                        kargs...) where D
MRIReco.SamplingOp
MRIReco.SensitivityOp
MRIReco.SparseOp
MRIReco.RegularizedLeastSquares.WeightingOp
```

## Datatypes
```@docs
MRIReco.AcquisitionData
MRIReco.AcquisitionData(tr::T,kdata::Array{Matrix{ComplexF64},3}
                        ; seqInfo=Dict{Symbol,Any}()
                        , idx=nothing
                        , encodingSize=Int64[0,0,0]
                        , fov=Float64[0,0,0]
                        , kargs...) where T <: Union{Trajectory,Vector{Trajectory}}
MRIReco.trajectory(acqData::AcquisitionData,i::Int64=1)
MRIReco.numContrasts(acqData::AcquisitionData)
MRIReco.numChannels
MRIReco.numSlices
MRIReco.numRepititions
MRIReco.kData
MRIReco.multiEchoData
MRIReco.multiCoilData
MRIReco.multiCoilMultiEchoData
MRIReco.profileData
MRIReco.samplingDensity
MRIReco.changeEncodingSize2D
MRIReco.convert3dTo2d
MRIReco.RawAcquisitionData
MRIReco.trajectory(f::RawAcquisitionData; slice::Int=1, contrast::Int=1)
MRIReco.rawdata(f::RawAcquisitionData)
MRIReco.AcquisitionData(f::RawAcquisitionData)
MRIReco.RawAcquisitionData(f::ISMRMRDFile, dataset="dataset")
MRIReco.AcquisitionData(f::ISMRMRDFile, dataset="dataset")
```

## Trajectories
```@docs
MRIReco.Trajectory
MRIReco.trajectory(trajName::AbstractString, numProfiles::Int, numSamplingPerProfile::Int; MRIReco.numSlices::Int64=1, TE::Float64=0.0, AQ::Float64=1.e-3, kargs...)
MRIReco.string(tr::Trajectory)
MRIReco.echoTime(tr::Trajectory)
MRIReco.acqTimePerProfile(tr::Trajectory)
MRIReco.numProfiles(tr::Trajectory)
MRIReco.numSamplingPerProfile(tr::Trajectory)
MRIReco.numSlices(tr::Trajectory)
MRIReco.isCircular(tr::Trajectory)
MRIReco.isCartesian(tr::Trajectory)
MRIReco.dims(tr::Trajectory)
MRIReco.kspaceNodes(tr::Trajectory)
MRIReco.readoutTimes(tr::Trajectory)
MRIReco.CartesianTrajectory
MRIReco.EPITrajectory
MRIReco.OneLine2dTrajectory
MRIReco.RadialTrajectory
MRIReco.SpiralTrajectory
MRIReco.SpiralTrajectoryVarDens
MRIReco.CartesianTrajectory3D
MRIReco.KooshballTrajectory
MRIReco.StackOfStarsTrajectory
```

## Sequences
```@docs
MRIReco.MESequence
MRIReco.numContrasts(seq::MESequence)
MRIReco.echoTimes(seq::MESequence)
MRIReco.flipAngles(seq::MESequence)
MRIReco.echoAmplitudes(seq::MESequence, R1::Float64, R2::Float64, numStates=nothing)
MRIReco.epgAmplitudes(seq::MESequence, R1::Real, R2::Real, numStates=nothing)
MRIReco.epgRotation
MRIReco.epgRelaxation
MRIReco.epgDephasing
MRIReco.rfRotation
```

## Sampling
```@docs
MRIReco.sample
MRIReco.sample_kspace(data::AbstractArray,redFac::Float64,patFunc::AbstractString;kargs...)
MRIReco.sample_kspace(acqData::AcquisitionData,redFac::Float64,
                       patFunc::AbstractString; rand=true, profiles=true,
                       seed = 1234, kargs...)
MRIReco.sample_regular(shape::Tuple, redFac::Float64; kargs...)
MRIReco.sample_random(shape::Tuple{Int64,Int64},redFac::Float64;calsize::Int64=0,kargs...)
MRIReco.sample_poissondisk(shape::Tuple{Int64,Int64},redFac::Float64;calsize::Int64=0, seed::Int64=1234,kargs...)
MRIReco.sample_vdpoisson(shape::Tuple{Int64,Int64},redFac::Float64; seed::Int64=1234,kargs...)
MRIReco.sample_lines(shape::Tuple{Int64,Int64},redFac::Float64;sampleFunc="random",kargs...)
MRIReco.calculateIncoherence(acqData::AcquisitionData, recoParams::Dict, slice=1)
```

## Simulation
```@docs
MRIReco.simulation(image::Array{T,3}, simParams::Dict) where T<:Union{ComplexF64,Float64}
MRIReco.simulation(image::Array{T,3}, simParams::Dict, filename::String;
                        force=false) where T<:Union{ComplexF64,Float64}
MRIReco.simulation(image::Array{T,2}, simParams::Dict) where T<:Union{ComplexF64,Float64}
MRIReco.simulation(tr::Trajectory
                    , image::Array{ComplexF64}
                    , correctionMap = []
                    ; opName="fast"
                    , senseMaps=[]
                    , verbose=true
                    , kargs...)
MRIReco.simulation(seq::AbstractSequence, tr::Vector{Trajectory}
                    , image::Array{ComplexF64,3}
                    ; opName="fast"
                    , r1map=[]
                    , r2map=[]
                    , fmap=[]
                    , senseMaps=[]
                    , verbose=true
                    , kargs...)
MRIReco.addNoise(x::Vector, snr::Float64, complex= true)
MRIReco.addNoise(acqData::AcquisitionData, snr::Float64)
MRIReco.addNoise!(acqData::AcquisitionData, snr::Float64)
MRIReco.birdcageSensitivity
MRIReco.quadraticFieldmap
```

## Reconstruction
```@docs
MRIReco.reconstruction(acqData::AcquisitionData, recoParams::Dict)
MRIReco.reconstruction(acqData::AcquisitionData, recoParams::Dict, filename::String;force=false)
MRIReco.setupIterativeReco
MRIReco.reconstruction_direct_2d
MRIReco.reconstruction_direct_3d
MRIReco.reconstruction_simple
MRIReco.reconstruction_multiEcho
MRIReco.reconstruction_multiCoil
MRIReco.reconstruction_multiCoilMultiEcho
MRIReco.espirit
MRIReco.nrmsd
```
