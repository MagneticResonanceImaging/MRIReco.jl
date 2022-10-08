export simulation,simulation_fast,simulation_explicit,simulateTempSubspace,addNoise

include("Fieldmap.jl")
include("RelaxationMap.jl")
include("CoilSensitivity.jl")
include("LoydMaxAlg.jl")
include("ExpApproximation.jl")

"""
    simulation(image::Array{T,3}, simParams::Dict) where T<:AbstractFloat

Simulate MRI raw data from given `image` data. All simulation parameters
are passed to the function in the form of a dictionary.
"""
function simulation(image::Array{T,3}, simParams::Dict) where T<:AbstractFloat
  return simulation(Complex.(image), simParams)
end

"""
    simulation(image::Array{T,3}, simParams::Dict) where T<:Complex{<:AbstractFloat}

Simulate MRI raw data from given `image` data. All simulation parameters
are passed to the function in the form of a dictionary.
"""
function simulation(image::Array{Complex{T},3}, simParams::Dict) where T<:AbstractFloat
  haskey(simParams, :correctionMap) ? cmap = reshape(simParams[:correctionMap],size(image)) : cmap = zeros(T,size(image))

  seqName = get(simParams,:seqName,"ME")
  trajName = get(simParams,:trajName,"Cartesian")
  seq = sequence(seqName, simParams[:numProfiles], simParams[:numSamplingPerProfile]; simParams...)
  numSamp, numProf = get(simParams, :numSamplingPerProfile,1), get(simParams, :numProfiles,1)
  tr = [trajectory(T, trajName, numProf,numSamp; simParams...) for i=1:numContrasts(seq)]

  opName = get(simParams,:simulation,"fast")
  if opName!="fast" && opName!="explicit"
    error("simulation $(simParams[:simulation]) is not known...")
  end

  return simulation(seq, tr, image; opName=opName, r2map=real.(cmap), fmap=imag.(cmap), simParams...)
end

"""
    simulation(image::Array{T,3}, simParams::Dict, filename::String) where T<:Union{Complex{<:AbstractFloat},AbstractFloat}

Performs the same simulation as `simulation(image, simParams)` and saves the result in
a file with name `filename`
"""
function simulation(image::Array{T,3}, simParams::Dict, filename::String;
                        force=false) where T<:Union{Complex{<:AbstractFloat},AbstractFloat}
  if !force && isfile(filename)
    f = ISMRMRDFile(filename)
    return AcquisitionData(f)
  else
    acq = simulation(image, simParams)
    f = ISMRMRDFile(filename)
    save(f, RawAcquisitionData(acq))
    return acq
  end
end

"""
    simulation(image::Array{T,2}, simParams::Dict) where T<:Union{Complex{<:AbstractFloat},AbstractFloat}

Simulate MRI raw data from given `image` data. All simulation parameters
are passed to the function in the form of a dictionary.
"""
function simulation(image::Array{T,2}, simParams::Dict) where T<:Union{Complex{<:AbstractFloat},AbstractFloat}
  nx,ny = size(image)
  return simulation(reshape(image,nx,ny,1),simParams)
end

#####################################################################
# low level simulation routines
#####################################################################

"""
simulation2d(tr::Trajectory, image::Array{T,3}, correctionMap=[]
              ; opName="fast", senseMaps=[], verbose=true, kargs...) where T<:Complex{<:AbstractFloat}

Transforms a given image to k-space Domain for a 2d Acquisition.
The Fourier integrals can be evaluated exactly or using NFFT
Returns the demodulated signal.

...
# Arguments
* `tr::Trajectory`             - two-dimensional Trajectory
* `image::Array{T,3}` - image to be transformed
* (`correctionMap=[]`)         - sum of the field offresonance (imaginary) map and relaxation map (real)
* (`opName="fast")             - name of operator to use ("explicit" or "fast")
* (`sensmaps=[]`)              - array of coil sensitivities
* (`verbose=true`)             - prints the progress if true
* `kargs...`                   - addional keyword arguments
...

"""
function simulation2d(tr::Trajectory{T}, image::Array{Complex{T},3}, correctionMap::Array{Complex{T}}=Complex{T}[]
              ; opName="fast", senseMaps=[], verbose=true, kargs...) where T<:AbstractFloat

  nx,ny,nz = size(image)

  if isempty(correctionMap)
    disturbanceTerm = zeros(Complex{T},nx,ny,nz)
  else
    if size(correctionMap) != size(image)
      error("correctionMap and image should have the same size!")
    end
    disturbanceTerm = Complex{T}.(correctionMap)
  end

  if isempty(senseMaps)
    sensFac = ones(nx,ny,nz,1)
    nc=1
  else
    if size(senseMaps)[1:3] != (nx,ny,nz)
      error("senseMaps and image should have the same size $(size(senseMaps)[1:3]) != $((nx,ny,nz))!")
    end
    sensFac = senseMaps
    nc = size(sensFac,4)
  end

  nodes = kspaceNodes(tr)
  # kdata = zeros(Complex{T}, size(nodes,2),nc,nz)
  kdata = [zeros(Complex{T},size(nodes,2),nc) for echo=1:1, slice=1:nz, rep=1:1]
  if verbose==true
    p = Progress(nz*nc, 1, "Simulating data...")
  end

  @sync for z = 1:nz
    Threads.@spawn begin
      E = fourierEncodingOp((nx,ny),tr,opName;slice=z,correctionMap=disturbanceTerm,echoImage=false)
      for c = 1:nc
        kdata[1,z,1][:,c] .= E*vec(sensFac[:,:,z,c].*image[:,:,z])
        verbose && next!(p)
      end
    end
  end

  return AcquisitionData(tr, kdata,encodingSize=[nx,ny,nz])
end

"""
    simulation3d(tr::Trajectory, image::Array{T,3}, correctionMap=[];
              opName="fast", senseMaps=[], verbose=true, kargs...) where T<:Complex{<:AbstractFloat}

Transforms a given image to k-space Domain 3d Acquisition.
The Fourier integrals can be evaluated exactly or using NFFT
Returns the demodulated signal.

...
# Arguments
* `tr::Trajectory`             - three-dimensional Trajectory
* `image::Array{T,3}` - image to be transformed
* (`correctionMap=[]`)         - sum of the field offresonance (imaginary) map and relaxation map (real)
* (`opName="fast")             - name of operator to use ("explicit" or "fast")
* (`sensmaps=[]`)              - array of coil sensitivities
* (`verbose=true`)             - prints the progress if true
* `kargs...`                   - addional keyword arguments
...

"""
function simulation3d(tr::Trajectory{T}, image::Array{Complex{T},3}, correctionMap::Array{Complex{T}}=Complex{T}[];
              opName="fast", senseMaps=[], verbose=true, kargs...) where T<:AbstractFloat

  nx,ny,nz = size(image)

  if isempty(correctionMap)
    disturbanceTerm = zeros(Complex{T},nx,ny,nz)
  else
    if size(correctionMap) != size(image)
      error("correctionMap and image should have the same size!")
    end
    disturbanceTerm = Complex{T}.(correctionMap)
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
  # kdata = zeros(Complex{T}, size(nodes,2),nc)
  kdata = [zeros(Complex{T},size(nodes,2),nc) for echo=1:1, slice=1:1, rep=1:1]
  if verbose==true
    p = Progress(nc, 1, "Simulating data...")
  end

  E = fourierEncodingOp((nx,ny,nz),tr,opName;correctionMap=disturbanceTerm,echoImage=false,symmetrize=false)
  for c = 1:nc
    # kdata[:,c] = E*vec(sensFac[:,:,:,c].*image)
    kdata[1,1,1][:,c] .= E*vec(sensFac[:,:,:,c].*image)
    verbose && next!(p)
  end

  return AcquisitionData(tr,kdata,encodingSize=[nx,ny,nz])
end

"""
    simulation(seq::AbstractSequence, tr::Vector{Trajectory}, image::Array{T,3}
                    ; opName="fast", r1map=[], r2map=[], fmap=[], senseMaps=[]
                    , verbose=true, kargs...) where T<:Complex{<:AbstractFloat}

Simulate k-space data for all echoes of a pulse sequence.
The echo intensities are simulated using the EPG formalism
The Fourier integrals can be evaluated exactly or using NFFT

...
# Arguments
* `seq::AbstractSequence`      - pulse sequence
* `tr::Vector{Trajectory}`     - trajectories for all contrasts
* `image::Array{T,3}` - image to be transformed
* (`r1map=[]`)                 - R1 map of the object (real)
* (`r2map=[]`)                 - R2 map of the object (real) / (R2* for GRE sequences)
* (`fmap=[]`)                  - fieldmap (real)
* (`opName="fast")             - name of operator to use ("explicit" or "fast")
* (`sensmaps=[]`)              - array of coil sensitivities
* (`verbose=true`)             - prints the progress if true
* `kargs...`                   - addional keyword arguments
...

"""
function simulation(seq::AbstractSequence, tr::Vector{Trajectory{T}}
                    , image::Array{Complex{T},3}
                    ; opName="fast"
                    , r1map=[]
                    , r2map=[]
                    , fmap=[]
                    , senseMaps=[]
                    , verbose=true
                    , kargs...) where T<:AbstractFloat
  nx,ny,nz = size(image)
  ne = numContrasts(seq)

  correctionMap = zeros(Complex{T},nx,ny,nz)
  if isempty(r1map)
    r1map = zeros(nx,ny,nz)
  end
  if isempty(r2map)
    r2map = zeros(nx,ny,nz)
  else
    correctionMap = zeros(Complex{T}, nx,ny,nz) # T.(r2map)
  end
  if !isempty(fmap)
    correctionMap = correctionMap .+ 1im*fmap
  end

  # compute echo amplitudes
  ampl = zeros(Complex{T}, nx,ny,nz,ne )
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
  isempty(senseMaps) ? nc=1 : nc=size(senseMaps,4)
  out = Vector{Matrix{Complex{T}}}(undef,0)

  if verbose
    p = Progress(numContrasts(seq), 1, "Simulating data...")
  end
  # tr = [trajectory(seq,i) for i=1:numContrasts(seq)]
  for i = 1:ne
    te = echoTime(tr[i])
    nodes = kspaceNodes(tr[i])
    t = readoutTimes(tr[i])

    # include correction map and compensate for relaxation before TE,
    # which is taken into account by the EPG-Simulation
    img_scaled = ampl[:,:,:,i].*image.*exp.(r2map*te)
    if length(findall(x->isfinite(x),img_scaled)) != length(img_scaled)
      error("scaled image contains NaN")
    end

    acq = simulation(tr[i], ampl[:,:,:,i].*image, correctionMap; senseMaps=senseMaps, verbose=verbose, kargs...)
    
    append!(out, vec(acq.kdata))

    # if length(findall(x->isfinite(x),out[:,i,:])) != length(out[:,i,:])
    #   error("out contains NaN")
    # end
  end

  ndims(tr[1])==2 ? numSl=nz : numSl=1
  
  return AcquisitionData(tr, reshape(out,ne,numSl,1), encodingSize=[nx,ny,nz])
end

"""
    simulation(tr::Trajectory, image::Array{T}, correctionMap = []; opName="fast"
              , senseMaps=[], verbose=true, kargs...) where T<:Complex{<:AbstractFloat}

Transforms a given image to k-space Domain.
Dispatches whether the trajectory is 2d or 3d
The Fourier integrals can be evaluated exactly or using NFFT
Returns the demodulated signal.

...
# Arguments
* `tr::Trajectory`             - three-dimensional Trajectory
* `image::Array{T,3}` - image to be transformed
* (`correctionMap=[]`)         - sum of the field offresonance (imaginary) map and relaxation map (real)
* (`opName="fast")             - name of operator to use ("explicit" or "fast")
* (`sensmaps=[]`)              - array of coil sensitivities
* (`verbose=true`)             - prints the progress if true
* `kargs...`                   - addional keyword arguments
...

"""
function simulation(tr::Trajectory{T}
                    , image::Array{Complex{T}}
                    , correctionMap::Array{Complex{T}} = Complex{T}[]
                    ; opName="fast"
                    , senseMaps=[]
                    , kargs...) where T<:AbstractFloat

  ndims(image) > 2 ? numSlices=size(image,3) : numSlices=1
  image = reshape(image,size(image)[1],size(image)[2],numSlices)
  if !isempty(correctionMap)
    correctionMap = reshape(correctionMap,size(image))
  end

  if ndims(tr)==2
    acqData = simulation2d(tr, image, correctionMap;opName="fast", senseMaps=senseMaps,kargs...)
  else
    acqData = simulation3d(tr, image, correctionMap;opName="fast", senseMaps=senseMaps,kargs...)
  end

  return acqData
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
  signal = Array(Float64,numContrasts(seq),Nr1*Nr2)
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
* `tr::Abstract2DTrajectory` - Two-dimensional Trajectory
* `image::Matrix`            - Image to be transformed
* `correctionMap=[]`         - Correctionmap is summed up of the field offresonance and relaxation
* `alpha::Float64 = 1.75`    - Oversampling factor
* `m::Float64=4.0`           - Kenrel size
* `K::Int64=28`              - Rank how much coefficients will be used to approx. correctionterm
* `method="nfft"`            - Which Method used to get the approx. coefficients of correctionterm
* (`sensmaps=[]`)            - array of coil sensitivities
* `kargs...`                 - addional keyword arguments
...

"""
function simulation_fast(tr::Trajectory{T}
                        , image::Matrix
                        , correctionMap::Array{Complex{T}} = Complex{T}[]
                        ; alpha::Float64 = 1.75
                        , m::Float64=4.0
                        , K::Int64=16
                        , method="nfft"
                        , senseMaps = []
                        , kargs...) where T
  image = reshape(image,size(image)[1],size(image)[2],1)
  if !isempty(correctionMap)
    correctionMap = reshape(correctionMap,size(image))
  end

  return simulation(tr, Complex{T}.(image), correctionMap;opName="fast", alpha=alpha,m=m,K=K,method=method,senseMaps=senseMaps,kargs...)
end

"""
  Getting (simulated) k-data using an explicit evaluation of the encoding operator

...
# Arguments
* `tr::Abstract2DTrajectory` - Two-dimensional Trajectory
* `image::Matrix`            - Image to be transformed
* `correctionMap=[]`         - Correctionmap is summed up of the field offresonance and relaxation
* `alpha::Float64 = 1.75`    - Oversampling factor
* `m::Float64=4.0`           - Kenrel size
* `K::Int64=28`              - Rank how much coefficients will be used to approx. correctionterm
* `method="nfft"`            - Which Method used to get the approx. coefficients of correctionterm
* (`sensmaps=[]`)            - array of coil sensitivities
* `kargs...`                 - addional keyword arguments
...

"""
function simulation_explicit(tr::Trajectory{T}
                        , image::Matrix
                        , correctionMap::Array{Complex{T}} = Complex{T}[]
                        ; alpha::Float64 = 1.75
                        , m::Float64=4.0
                        , K::Int64=16
                        , method="nfft"
                        , senseMaps = []
                        , kargs...) where T

  image = reshape(image,size(image)[1],size(image)[2],1)
  if !isempty(correctionMap)
    correctionMap = reshape(correctionMap,size(image))
  end

  return simulation(tr, Complex{T}.(image), correctionMap;opName="explicit", alpha=alpha,m=m,K=K,method=method,senseMaps=senseMaps,kargs...)
end

"""
  Adds average white gaussian noise to the signal x

# Arguments
* `x::Vector`     - signal vector
* 'snr::Float64'  - target SNR

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

"""
  return AcquisitionData with white gaussian noise

# Arguments
* `acqData::AcquisitionData`  - AcquisitionData
* 'snr::Float64'              - target SNR

"""
function addNoise(acqData::AcquisitionData, snr::Float64)
  acqData2 = deepcopy(acqData)
  addNoise!(acqData2,snr)
  return acqData2
end

"""
  add white gaussian noise to AcquisitionData with (in-place)

# Arguments
* `acqData::AcquisitionData`  - AcquisitionData
* 'snr::Float64'              - target SNR

"""
function addNoise!(acqData::AcquisitionData, snr::Float64)
  for i=1:length(acqData.kdata)
    acqData.kdata[i][:] .= addNoise(vec(acqData.kdata[i]),snr,true)[:]
  end
end
