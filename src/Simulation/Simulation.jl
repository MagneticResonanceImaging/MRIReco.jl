export simulation,createSimulationProblem,simulation_fast,simulation_explicit,simulation_dynamic, simulateTempSubspace, simulation_dynamic_rigid, simulation_dynamic_rigid_parallel, AcquisitionData

include("Fieldmap.jl")
include("RelaxationMap.jl")
include("CoilSensitivity.jl")
include("ExpApproximation.jl")

#####################################################################
# simulation interface
#####################################################################
function simulation(image::Array{Float64}, simParams::Dict)
  haskey(simParams, :correctionMap) ? cmap = simParams[:correctionMap] : cmap = ComplexF64[]

  if simParams[:simulation] == "explicit"
    tr = trajectory(simParams[:trajName], simParams[:numProfiles], simParams[:numSamplingPerProfile]; simParams...)
    return simulation_explicit(tr, image, cmap; simParams...)
  elseif simParams[:simulation] == "fast"
    tr = trajectory(simParams[:trajName], simParams[:numProfiles], simParams[:numSamplingPerProfile]; simParams...)
    return simulation_fast(tr, image, cmap; simParams...)
  elseif simParams[:simulation] == "epgExplicit"
    seq = sequence(simParams[:seqName], simParams[:trajName], simParams[:numProfiles], simParams[:numSamplingPerProfile]; simParams...)
    return   simulation_explicit(seq, image; r2map=real(cmap), fmap = imag(cmap), simParams...)
  elseif simParams[:simulation] == "epgNFFT"
    seq = sequence(simParams[:seqName], simParams[:trajName], simParams[:numProfiles], simParams[:numSamplingPerProfile]; simParams...)
    return   simulation_fast(seq, image; r2map=real(cmap), fmap = imag(cmap), simParams...)
  else
    error("simulation $(simParams[:simulation]) is not known...")
  end
end

# This version stores the simulation data into a file
function simulation(image::Array{Float64}, simParams::Dict, filename::String;
                        force=false)
  if !force && isfile(filename)
    return loadIBIFile(filename)
  else
    acq = simulation(image, simParams)
    saveasIBIFile(filename, acq)
    return acq
  end
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
function simulation_explicit(tr::Abstract2DTrajectory, image::Matrix, correctionMap=[]; verbose=true, kargs...)
  nodes = kspaceNodes(tr)
  kdata = zeros(ComplexF64, size(nodes,2))

  if verbose==true
    p = Progress(size(nodes,2), 1, "Simulating data...")
  end

  if isempty(correctionMap)
    disturbanceTerm = zeros(size(image)...)
  else
    disturbanceTerm = correctionMap
  end

  t = readoutTimes(tr)

  for k=1:size(nodes,2)
    for ny=1:size(image,2)
      for nx=1:size(image,1)
        phi = (nodes[1,k]*(nx-size(image,1)/2-1)+
               nodes[2,k]*(ny-size(image,2)/2-1))
        e = exp(-2*1im*pi*phi - disturbanceTerm[nx,ny]*t[k])
        kdata[k] += image[nx,ny] * e
      end
    end
    if verbose==true
      next!(p)
    end
  end


  return AcquisitionData(tr, kdata)

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
function simulation_explicit(tr::Abstract3DTrajectory, image::Array{Float64,3}, correctionMap=[];verbose=true, kargs...)
  nodes = kspaceNodes(tr)
  out = zeros(ComplexF64, size(nodes,2))

  if verbose==true
    p = Progress(size(nodes,2), 1, "Simulating data...")
  end

  # Readout times, i.e. effects of Relaxation and Offresonce not considered in 3D
  # (for now...)

  if isempty(correctionMap)
    disturbanceTerm = zeros(size(image)...)
  else
    disturbanceTerm = correctionMap
  end

  t = readoutTimes(tr)

  for k=1:size(nodes,2)
    for nz=1:size(image,3)
      for ny=1:size(image,2)
        for nx=1:size(image,1)
          phi = (nodes[1,k]*(nx-size(image,1)/2-1)+
                 nodes[2,k]*(ny-size(image,2)/2-1)+
                 nodes[3,k]*(nz-size(image,3)/2-1))
          e = exp(-2*1im*pi*phi - disturbanceTerm[nx,ny,nz]*t[k])
          out[k] += image[nx,ny,nz] * e
        end
      end
    end
    if verbose==true
      next!(p)
    end
  end

  return AcquisitionData(tr, out)
end

function simulation_explicit(seq::AbstractSequence
                    , image::Array{Float64,3}
                    ; r1map=[]
                    , r2map=[]
                    , fmap=[]
                    , verbose=true
                    , kargs...)
  if isempty(r1map)
    r1map = zeros(size(image)...)
  end
  if isempty(r2map)
    r2map = zeros(size(image)...)
  end
  if isempty(fmap)
    fmap = zeros(size(image)...)
  end

  size_kNodes = size(kspaceNodes( trajectory(seq) ), 2)
  out = zeros( ComplexF64, size_kNodes, numEchoes(seq) )

  if verbose
    p = Progress(size(image,3)*numEchoes(seq), 1, "Simulating data...")
  end

  # compute echo amplitudes
  ampl = zeros(ComplexF64, size(image,1), size(image,2), size(image,3), numEchoes(seq) )
  for nz = 1:size(image,3)
    for ny=1:size(image,2)
      for nx=1:size(image,1)
        ampl[nx,ny,nz,:]= echoAmplitudes(seq, r1map[nx,ny,nz], r2map[nx,ny,nz] )
      end
    end
  end

  for i = 1:numEchoes(seq)
    tr = trajectory(seq,i)
    nodes = kspaceNodes(tr)
    t = readoutTimes(tr)

    for nz=1:size(image,3)
      for ny=1:size(image,2)
        for nx=1:size(image,1)
          for k=1:size_kNodes
            phi = (nodes[1,k]*(nx-size(image,1)/2-1)+
                   nodes[2,k]*(ny-size(image,2)/2-1)+
                   nodes[3,k]*(nz-size(image,3)/2-1))

            # include correction map and compensate for relaxation before TE,
            # which is taken into account by the EPG-Simulation
            e = exp( -1im*( 2*pi*phi - fmap[nx,ny,nz]*(t[k]-tr.TE)) - rmap[nx,ny,nz]*(t[k]-tr.TE) )

            out[k,i] += ampl[nx,ny,nz,i] * image[nx,ny,nz] * e
          end
        end
      end
      if verbose
        next!(p)
      end
    end
  end

  return AcquisitionData(seq, vec(out), numEchoes=numEchoes)
end

function simulation_explicit(seq::AbstractSequence
                      , image::Array{Float64,2}
                      ; r1map=[]
                      , r2map=[]
                      , fmap=[]
                      , verbose=true
                      , kargs...)
  if isempty(r1map)
    r1map = zeros(size(image)...)
  end
  if isempty(r2map)
    r2map = zeros(size(image)...)
  end
  if isempty(fmap)
    fmap = zeros(size(image)...)
  end

  tr = trajectory(seq)
  size_kNodes = size( kspaceNodes(tr), 2 )
  out = zeros( ComplexF64, size_kNodes, numEchoes(seq) )

  if verbose
    p = Progress(size(image,2)*numEchoes(seq), 1, "Simulating data...")
  end

  # compute echo amplitudes
  ampl = zeros(ComplexF64, size(image,1), size(image,2), numEchoes(seq) )
  for ny=1:size(image,2)
    for nx=1:size(image,1)
      ampl[nx,ny,:]= abs( echoAmplitudes(seq, r1map[nx,ny], r2map[nx,ny] ) )
    end
  end

  for i = 1:numEchoes(seq)
    tr = trajectory(seq,i)
    nodes = kspaceNodes(tr)
    t = readoutTimes(tr)

    for ny=1:size(image,2)
      for nx=1:size(image,1)
        for k=1:size_kNodes
          phi = (nodes[1,k]*(nx-size(image,1)/2-1)+nodes[2,k]*(ny-size(image,2)/2-1))

          # include correction map and compensate for relaxation before TE,
          # which is taken into account by the EPG-Simulation
          e = exp( -1im*( 2*pi*phi - fmap[nx,ny]*(t[k]-tr.TE)) - rmap[nx,ny]*(t[k]-tr.TE) )

          out[k,i] += ampl[nx,ny,i] * image[nx,ny] * e
        end
      end
      if verbose
        next!(p)
      end
    end
  end

  return AcquisitionData(seq, vec(out), numEchoes=numEchoes)
end

function simulation_fast(seq::AbstractSequence
                      , image::Array{Float64,2}
                      ; r1map=[]
                      , r2map=[]
                      , fmap=[]
                      , senseMaps = []
                      , verbose=true
                      , kargs...)
  if isempty(r1map)
    r1map = zeros(size(image)...)
  end
  if isempty(r2map)
    r2map = zeros(size(image)...)
  end
  if isempty(fmap)
    fmap = zeros(size(image)...)
  end

  tr = trajectory(seq)
  size_kNodes = size( kspaceNodes(tr), 2 )
  senseMaps!=[] ? numCoils= div(length(senseMaps),length(image)) : numCoils=1
  out = zeros( ComplexF64, size_kNodes, numEchoes(seq), numCoils )

  if verbose
    p = Progress(size(image,2)*numEchoes(seq), 1, "Simulating data...")
  end

  # compute echo amplitudes
  ampl = zeros(ComplexF64, size(image,1), size(image,2), numEchoes(seq) )
  for ny=1:size(image,2)
    for nx=1:size(image,1)
      ampl[nx,ny,:]= abs( echoAmplitudes(seq, r1map[nx,ny], r2map[nx,ny] ) )
    end
  end

  # simulate kspace
  for i = 1:numEchoes(seq)
    println("simulate echo $i")
    tr = trajectory(seq,i)
    # include correction map and compensate for relaxation before TE,
    # which is taken into account by the EPG-Simulation
    out[:,i,:] = simulation_fast(tr, ampl[:,:,i].*exp(r2map*tr.TE).*image, r2map+1im*fmap, senseMaps=senseMaps).kdata
  end

  println("numCoils = $numCoils")

  return AcquisitionData(seq, vec(out), numEchoes=numEchoes(seq), numCoils=numCoils)
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
  # Get kspace nodes and readouttimes
  nodes = kspaceNodes(tr)
  numOfNodes = size(nodes,2)

  times = readoutTimes(tr)
  imageShape = size(image)
  numOfPixel = prod(imageShape)

  kdata = zeros(ComplexF64,numOfNodes)
  if isempty(correctionMap)
    nfftOp = NFFTOp(imageShape, tr;symmetrize=false)
  else
    # Get coefficients for exp. of Correctionterm
    nfftOp = FieldmapNFFTOp(imageShape,tr,correctionMap; method=method, echoImage=false, symmetrize=false)
  end

  if isempty(senseMaps)
    kdata = nfftOp * vec(complex(image))
    return AcquisitionData(tr, kdata)
  else
    senseMaps = reshape(senseMaps,length(image),:)
    kdata = zeros(ComplexF64, numOfNodes, size(senseMaps,2))
    for l=1:size(senseMaps,2)
      kdata[:,l] = nfftOp * (vec(complex(image)) .* vec(senseMaps[:,l]))
    end
    return AcquisitionData(tr, vec(kdata), numCoils=size(senseMaps,2))
  end

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

function addNoise(aqData::AcquisitionData, snr::Float64)
  noisyData = addNoise(aqData.kdata, snr, true)

  return AcquisitionData(aqData.seq, noisyData, aqData.numEchoes, aqData.numCoils, aqData.numSlices, aqData.samplePointer, aqData.idx)
end

"""
  Returns EchoImages. Each echoimage is influenced by the correctionterm and the echotime
  (time offset where recording of the signal had been started)
  z   : correctionterm
  tau : echoterm

  Basically the formula of echoimage is : I_echo = I_orig .* exp(-z * tau)

"""
function simulateEchoImages(image::Matrix, echoTimes::Vector, cmap=[])
  echoImg = zeros(ComplexF64, size(image,1),size(image,2),length(echoTimes))
  if isempty(cmap)
    for i=1:length(echoTimes)
      echoImg[:,:,i] = image
    end
  else
    for i=1:length(echoTimes)
      echoImg[:,:,i] = image .* exp(-cmap .* echoTimes[i])
    end
  end
  return echoImg
end
