export simulation,simulation_fast,simulation_explicit,simulateTempSubspace,addNoise

include("Fieldmap.jl")
include("RelaxationMap.jl")
include("CoilSensitivity.jl")
include("ExpApproximation.jl")

"""
  simulation(image::Array{Float64}, simParams::Dict)

Simulate MRI raw data from given `image` data. All simulation parameters
are passed to the function in the form of a dictionary.
"""
function simulation(image::Array{T,3}, simParams::Dict) where T<:Union{ComplexF64,Float64}
  haskey(simParams, :correctionMap) ? cmap = reshape(simParams[:correctionMap],size(image)) : cmap = zeros(ComplexF64,size(image))

  seqName = get(simParams,:seqName,"SE")
  trajName = get(simParams,:trajName,"Cartesian")
  seq = sequence(seqName, trajName, simParams[:numProfiles], simParams[:numSamplingPerProfile]; simParams...)

  opName = get(simParams,:simulation,"fast")
  if opName!="fast" && opName!="explicit"
    error("simulation $(simParams[:simulation]) is not known...")
  end

  return simulation(seq, ComplexF64.(image); opName=opName, r2map=real.(cmap), fmap=imag.(cmap), simParams...)
end

# This version stores the simulation data into a file
function simulation(image::Array{T,3}, simParams::Dict, filename::String;
                        force=false) where T<:Union{ComplexF64,Float64}
  if !force && isfile(filename)
    return loadIBIFile(filename)
  else
    acq = simulation(image, simParams)
    saveasIBIFile(filename, acq)
    return acq
  end
end

function simulation(image::Array{T,2}, simParams::Dict) where T<:Union{ComplexF64,Float64}
  nx,ny = size(image)
  return simulation(reshape(image,nx,ny,1),simParams)
end

#####################################################################
# low level simulation routines
#####################################################################

"""
    Transforms a given image to k-space Domain using an exact
    evaluation of the discrete integral.
    Returns the demodulated signal.

...
# Arguments
* `tr::Abstract2DTrajectory` : Two-dimensional Trajectory
* `image::Matrix` : Image to be transformed
* `correctionMap=[]` : Correctionmap is summed up of the field offresonance and relaxation
* `verbose=true` : Prints the progress if true
...
# Examples
```jldoctest
julia> N=64; I = shepp_logan(N); tr = CartesianTrajectory(N,N);
julia> kdata = simulation(tr,I);

```

"""
function simulation(tr::Abstract2DTrajectory, image::Array{ComplexF64,3}, correctionMap=[]
              ; opName="fast", senseMaps=[], verbose=true, kargs...) #where T<:Union{ComplexF64,Float64}

  nx,ny,nz = size(image)

  if isempty(correctionMap)
    disturbanceTerm = zeros(ComplexF64,nx,ny,nz)
  else
    if size(correctionMap) != size(image)
      error("correctionMap and image should have the same size!")
    end
    disturbanceTerm = ComplexF64.(correctionMap)
  end

  if isempty(senseMaps)
    sensFac = ones(nx,ny,nz,1)
    nc=1
  else
    if size(senseMaps)[1:3] != (nx,ny,nz)
      error("senseMaps and image should have the same size!")
    end
    sensFac = senseMaps
    nc = size(sensFac,4)
  end

  nodes = kspaceNodes(tr)
  kdata = zeros(ComplexF64, size(nodes,2),nc,nz)
  if verbose==true
    p = Progress(nz*nc, 1, "Simulating data...")
  end
  for z = 1:nz
    E = fourierEncodingOp2d((nx,ny),tr,opName;slice=z,correctionMap=disturbanceTerm,echoImage=false,symmetrize=false)
    for c = 1:nc
      kdata[:,c,z] = E*vec(sensFac[:,:,z,c].*image[:,:,z])
      verbose && next!(p)
    end
  end

  return AcquisitionData(tr, vec(kdata),numCoils=nc,numSlices=nz)
end

"""
    Transforms a given image to k-space Domain using an exact
    evaluation of the discrete integral.
    Returns the demodulated signal.

...
# Arguments
* `tr::Abstract3DTrajectory` : Three-dimensional Trajectory
* `image::Matrix` : Image to be transformed
* `correctionMap=[]` : Correctionmap is summed up of the field offresonance and relaxation
* `verbose=true` : Prints the progress if true
...

"""
function simulation(tr::Abstract3DTrajectory, image::Array{ComplexF64,3}, correctionMap=[];
              opName="fast", senseMaps=[], verbose=true, kargs...) # where T<:Union{ComplexF64,Float64}

  nx,ny,nz = size(image)

  if isempty(correctionMap)
    disturbanceTerm = zeros(ComplexF64,nx,ny,nz)
  else
    if size(correctionMap) != size(image)
      error("correctionMap and image should have the same size!")
    end
    disturbanceTerm = ComplexF64.(correctionMap)
  end

  if isempty(senseMaps)
    sensFac = ones(nx,ny,nz,1)
    nc=1
  else
    if size(senseMaps)[1:3] != (nx,ny,nz)
      error("senseMaps and image should have the same size!")
    end
    sensFac = senseMaps
    nc = size(sensFac,4)
  end

  nodes = kspaceNodes(tr)
  kdata = zeros(ComplexF64, size(nodes,2),nc)
  if verbose==true
    p = Progress(nc, 1, "Simulating data...")
  end

  E = fourierEncodingOp3d((nx,ny,nz),tr,opName;correctionMap=disturbanceTerm,echoImage=false,symmetrize=false)
  for c = 1:nc
    kdata[:,c] = E*vec(sensFac[:,:,:,c].*image)
    verbose && next!(p)
  end

  return AcquisitionData(tr, vec(kdata) ,numCoils=nc,numSlices=nz)
end

function simulation(seq::AbstractSequence
                    , image::Array{ComplexF64,3}
                    ; opName="fast"
                    , r1map=[]
                    , r2map=[]
                    , fmap=[]
                    , senseMaps=[]
                    , verbose=true
                    , kargs...)
  nx,ny,nz = size(image)
  ne = numEchoes(seq)

  correctionMap = zeros(ComplexF64,nx,ny,nz)
  if isempty(r1map)
    r1map = zeros(nx,ny,nz)
  end
  if isempty(r2map)
    r2map = zeros(nx,ny,nz)
  else
    correctionMap = zeros(nx,ny,nz) # ComplexF64.(r2map)
  end
  if !isempty(fmap)
    correctionMap = correctionMap .+ 1im*fmap
  end

  # compute echo amplitudes
  ampl = zeros(ComplexF64, nx,ny,nz,ne )
  if verbose
    p = Progress(nx*ny*nz,1,"Compute echo amplitudes ")
  end
  for k = 1:nz
    for j=1:ny
      for i=1:nx
        ampl[i,j,k,:]= echoAmplitudes(seq, r1map[i,j,k], r2map[i,j,k] )
      end
    end
  end

  # this assumes the same number of readout points per echo
  size_kNodes = size(kspaceNodes( trajectory(seq) ), 2)
  isempty(senseMaps) ? nc=1 : nc=size(senseMaps,4)
  out = zeros( ComplexF64, size_kNodes, ne, nc )

  if verbose
    p = Progress(numEchoes(seq), 1, "Simulating data...")
  end
  for i = 1:ne
    tr = trajectory(seq,i)
    te = echoTime(tr)
    nodes = kspaceNodes(tr)
    t = readoutTimes(tr)

    # include correction map and compensate for relaxation before TE,
    # which is taken into account by the EPG-Simulation
    img_scaled = ampl[:,:,:,i].*image.*exp.(r2map*te)
    if length(findall(x->isfinite(x),img_scaled)) != length(img_scaled)
      error("scaled image contains NaN")
    end

    # out[:,i,:] = simulation(tr, ampl[:,:,:,i].*image.*exp.(r2map*te), correctionMap; senseMaps=senseMaps, verbose=true, kargs...).kdata
    out[:,i,:] = simulation(tr, ampl[:,:,:,i].*image, correctionMap; senseMaps=senseMaps, verbose=true, kargs...).kdata

    if length(findall(x->isfinite(x),out[:,i,:])) != length(out[:,i,:])
      error("out contains NaN")
    end
  end

  return AcquisitionData(seq, vec(out), numEchoes=ne, numCoils=nc, numSlices=nz)
end

#
# simulate temporal signals for a given parameter map and
# construct temporal subspace using SVD.
#
function simulateTempSubspace(seq::AbstractSequence
                              , r1sample::Array{Float64}
                              , r2sample::Array{Float64}
                              ; subsize::Int64=4)

  # simulate training data
  Nr1 = length(r1sample)
  Nr2 = length(r2sample)
  signal = Array(Float64,numEchoes(seq),Nr1*Nr2)
  param = Array(Float64,Nr1*Nr2,2)                # simulated parameters
  for i=1:Nr1, j=1:Nr2
    signal[:,(i-1)*Nr2+j] = abs( echoAmplitudes( seq, r1sample[i], r2sample[j] ) )
    param[(i-1)*Nr2+j,1] =  r1sample[i]
    param[(i-1)*Nr2+j,2] =  r2sample[j]
  end

  # calculate temporal subspace
  U,S,V = svd(signal)

  return U[:,1:subsize]#, signal, param
end

"""
    Getting (simulated) k-data with the help of NFFTs with lower computational effort

...
# Arguments
* `tr::Abstract2DTrajectory` : Two-dimensional Trajectory
* `image::Matrix` : Image to be transformed
* `correctionMap=[]` : Correctionmap is summed up of the field offresonance and relaxation
* `alpha::Float64 = 1.75` : Oversampling factor
* `m::Float64=4.0` : Kenrel size
* `K::Int64=28` : Rank how much coefficients will be used to approx. correctionterm
* `method="nfft"` : Which Method used to get the approx. coefficients of correctionterm
...

"""
function simulation_fast(tr::Abstract2DTrajectory
                        , image::Matrix
                        , correctionMap = []
                        ; alpha::Float64 = 1.75
                        , m::Float64=4.0
                        , K::Int64=16
                        , method="nfft"
                        , senseMaps = []
                        , kargs...)
  image = reshape(image,size(image)[1],size(image)[2],1)
  if !isempty(correctionMap)
    correctionMap = reshape(correctionMap,size(image))
  end

  return simulation(tr,ComplexF64.(image),correctionMap;opName="fast", alpha=alpha,m=m,K=K,method=method,senseMaps=senseMaps,kargs...)
end

function simulation_explicit(tr::Abstract2DTrajectory
                        , image::Matrix
                        , correctionMap = []
                        ; alpha::Float64 = 1.75
                        , m::Float64=4.0
                        , K::Int64=16
                        , method="nfft"
                        , senseMaps = []
                        , kargs...)
  image = reshape(image,size(image)[1],size(image)[2],1)
  if !isempty(correctionMap)
    correctionMap = reshape(correctionMap,size(image))
  end
  return simulation(tr,ComplexF64.(image),correctionMap;opName="explicit", alpha=alpha,m=m,K=K,method=method,senseMaps=senseMaps,kargs...)
end

"""
  Adds average white gaussian noise to the signal

"""
function addNoise(x::Vector, snr::Float64, complex= true)
  signalAmpl = sum(abs.(x))/length(x)

  # target noise amplitude
  noiseAmpl = signalAmpl/snr
  if complex
    noise = noiseAmpl/sqrt(2.)*( randn(size(x))+ 1im*randn(size(x)) )
  else
    noise = noiseAmpl*randn(size(x))
  end

  return x+noise
end

function addNoise(acqData::AcquisitionData, snr::Float64)
  noisyData = addNoise(acqData.kdata, snr, true)

  return AcquisitionData(acqData.seq, noisyData, acqData.numEchoes,
                         acqData.numCoils, acqData.numSlices, acqData.samplePointer,
                         acqData.subsampleIndices, acqData.encodingSize, acqData.fov)
end
