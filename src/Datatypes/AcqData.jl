export AcquisitionData, kData, kdataSingleSlice, convertUndersampledData,
       changeEncodingSize2D, convert3dTo2d, samplingDensity,
       numContrasts, numChannels, numSlices, numRepititions, apodization!

"""
struct describing MRI acquisition data.

# Fields
* `sequenceInfo::Dict{Symbol,Any}`          - additional information on the pulse sequence
* `traj::Vector{Trajectory}`                - trajectories for each echo/contrast
* `kdata::Array{Matrix{ComplexF64},3}`      - each matrix contains data for one trajectory
                                              (1. dim k-space nodes, 2. dim coils)
                                              the outer dims describe:
                                              1. dim echoes, 2. dim slices, 3. dim repetitions
* `subsampleIndices::Vector{Array{Int64}}`  - indices sampled for each echo/contrast
* `encodingSize::Vector{Int64}`             - size of the underlying image matrix
* `fov::Vector{Float64}`                    - field of view in m
"""
mutable struct AcquisitionData
  sequenceInfo::Dict{Symbol,Any}
  traj::Vector{Trajectory}
  kdata::Array{Matrix{ComplexF64},3}
  subsampleIndices::Vector{Vector{Int64}}
  encodingSize::Vector{Int64}
  fov::Vector{Float64}
end

"""
    numChannels(acqData::AcquisitionData)

returns the number of channels/coils in acqData
"""
numChannels(acqData::AcquisitionData) = size(acqData.kdata[1],2)

"""
    numContrasts(acqData::AcquisitionData)

returns the number of contrasts/echoes in acqData
"""
numContrasts(acqData::AcquisitionData) = size(acqData.kdata,1)

"""
    numSlices(acqData::AcquisitionData)

returns the number of slices in acqData
"""
numSlices(acqData::AcquisitionData) = size(acqData.kdata,2)

"""
    numRepititions(acqData::AcquisitionData)

returns the number of repetitions in acqData
"""
numRepititions(acqData::AcquisitionData) = size(acqData.kdata,3)

"""
    AcquisitionData(tr::T,kdata::Array{Matrix{ComplexF64},3}; seqInfo=Dict{Symbol,Any}()
                        , idx=nothing, encodingSize=Int64[0,0,0], fov=Float64[0,0,0]
                        , kargs...) where T <: Union{Trajectory,Vector{Trajectory}}

constructor for `AcquisitionData`

# Arguments
* `tr <: Union{Trajectory,Vector{Trajectory}}` - trajectories
* `kdata::Array{Matrix{ComplexF64},3}`         - k-space data

the other fields of `AcquisitionData` can be passed as keyword arguments.
"""
function AcquisitionData(tr::T,kdata::Array{Matrix{ComplexF64},3}
                        ; seqInfo=Dict{Symbol,Any}()
                        , idx=nothing
                        , encodingSize=Int64[0,0,0]
                        , fov=Float64[0,0,0]
                        , kargs...) where T <: Union{Trajectory,Vector{Trajectory}}
  tr_vec = vec(tr)
  if idx != nothing
    subsampleIndices = idx
  else
    numContr = size(kdata,1)
    if length(tr_vec) == numContr
      subsampleIndices = [collect(1:size(kspaceNodes(tr_vec[echo]),2)) for echo=1:numContr]
    else
      numSamp = size(kspaceNodes(tr_vec[1]),2)
      subsampleIndices = [collect(1:numSamp) for echo=1:numContr]
    end
  end

  return AcquisitionData(seqInfo,tr_vec,kdata,subsampleIndices,encodingSize,fov)
end

function Images.pixelspacing(acqData::AcquisitionData)
  return [1.0,1.0,1.0]*Unitful.mm
  #return fov./encodingSize*Unitful.mm  #TODO: all needs to be properly initialized
end

"""
    trajectory(acqData::AcquisitionData,i::Int64=1)

returns the `i`-th trajectory contained in `acqData`.
"""
trajectory(acqData::AcquisitionData,i::Int64=1) = acqData.traj[i]

######################
# getting k-space data
######################
"""
    kData(acqData::AcquisitionData, echo::Int64=1, coil::Int64=1, slice::Int64=1;rep::Int64=1)

returns the k-space contained in `acqData` for given `echo`, `coil`, `slice` and `rep`(etition).
"""
function kData(acqData::AcquisitionData, echo::Int64=1, coil::Int64=1, slice::Int64=1;rep::Int64=1)
  return acqData.kdata[echo,slice,rep][:,coil]
end

"""
    multiEchoData(acqData::AcquisitionData, coil::Int64, slice::Int64;rep::Int64=1)

returns the k-space contained in `acqData` for all echoes and given `coil`, `slice` and `rep`(etition).
"""
function multiEchoData(acqData::AcquisitionData, coil::Int64, slice::Int64;rep::Int64=1)
  kdata = ComplexF64[]
  for echo=1:numContrasts(acqData)
    append!(kdata,acqData.kdata[echo,slice,rep][:,coil])
  end
  return kdata
end

"""
    multiCoilData(acqData::AcquisitionData, echo::Int64, slice::Int64;rep::Int64=1)

returns the k-space contained in `acqData` for all coils and given `echo`, `slice` and `rep`(etition).
"""
function multiCoilData(acqData::AcquisitionData, echo::Int64, slice::Int64;rep::Int64=1)
  return vec(acqData.kdata[echo,slice,rep])
end

"""
    multiCoilMultiEchoData(acqData::AcquisitionData, echo::Int64, slice::Int64;rep::Int64=1)

returns the k-space contained in `acqData` for all coils, echoes and given `slice` and `rep`(etition).
"""
function multiCoilMultiEchoData(acqData::AcquisitionData,slice::Int64;rep=1)
  kdata = ComplexF64[]
  for coil=1:numChannels(acqData)
    for echo=1:numContrasts(acqData)
      append!(kdata, acqData.kdata[echo,slice,rep][:,coil])
    end
  end
  return kdata
end

"""
    profileData(acqData::AcquisitionData, echo::Int64, slice::Int64, rep::Int, prof_tr::Int)

returns the profile-data `prof_tr` contained in `acqData` for given `echo`, `coil`, `slice` and `rep`(etition).
"""
function profileData(acqData::AcquisitionData, echo::Int64, slice::Int64, rep::Int, prof_tr::Int)
  tr = trajectory(acqData,echo)
  numSamp, numSl = numSamplingPerProfile(tr), numSlices(tr)
  numChan = numChannels(acqData)
  numProf = div(length(acqData.subsampleIndices[echo]),numSamp) #numProfiles(tr)
  if dims(tr)==2 || numSl==1
    kdata = reshape(multiCoilData(acqData,echo,slice;rep=rep),numSamp,numProf,numChan)
    prof_data = kdata[:,prof_tr,:]
  else
    kdata = reshape(multiCoilData(acqData,echo,1,rep=rep),numSamp,numProf,numSl,numChan)
    prof_data = kdata[:,prof_tr,slice,:]
  end
  return prof_data
end

######################################
# utilities to convert and edit acqData
######################################
"""
    convertUndersampledData(acqData::AcquisitionData)

converts undersampled AcquisitionData, where only profiles contained in
acqData.subsampleIndices are sampled,
into a format where trajectories only contain the sampled profiles.
"""
function convertUndersampledData(acqData::AcquisitionData)

  acqDataSub = deepcopy(acqData)
  numContr = numContrasts(acqData)

  # get number of nodes and reset idx
  numNodes = size(acqData.subsampleIndices,1)
  for echo=1:numContr
    acqDataSub.subsampleIndices[echo] = collect(1:length(acqData.subsampleIndices[echo]))
  end

  # replace trajectories by Undersampled Trajectories
  for i = 1:numContr
    tr = trajectory(acqDataSub,i)
    # assume that coils and slices experience the same trajectory
    tr.nodes = tr.nodes[:,acqData.subsampleIndices[i]]
    tr.cartesian = false
  end

  return acqDataSub
end

##################
# sampling weights
##################
"""
    samplingDensity(acqData::AcquisitionData,shape::Tuple)

returns the sampling density for all trajectories contained in `acqData`.
"""
function samplingDensity(acqData::AcquisitionData,shape::Tuple)
  numContr = numContrasts(acqData)
  weights = Array{Vector{ComplexF64}}(undef,numContr)
  for echo=1:numContr
    tr = trajectory(acqData,echo)
    if isCartesian(tr)
      nodes = kspaceNodes(tr)[:,acqData.subsampleIndices[echo]]
    else
      nodes = kspaceNodes(tr)
    end
    plan = NFFTPlan(nodes, shape,3, 1.25)
    weights[echo] = sqrt.(sdc(plan, iters=3))
  end
  return weights
end

#########################################################################
# convert acqData for a reconstruction with a encodingSize (resolution)
#########################################################################
"""
    changeEncodingSize2D(acqData::AcquisitionData,newEncodingSize::Vector{Int64})

changes the encoding size of 2d encoded `acqData` to `newEncodingSize`.
Returns a new `AcquisitionData` object.
"""
function changeEncodingSize2D(acqData::AcquisitionData,newEncodingSize::Vector{Int64})
  dest = deepcopy(acqData)
  changeEncodingSize2D!(dest,newEncodingSize)
end

"""
    changeEncodingSize2D!(acqData::AcquisitionData,newEncodingSize::Vector{Int64})

does the same thing as `changeEncodingSize2D` but acts in-place on `acqData`.
"""
function changeEncodingSize2D!(acqData::AcquisitionData,newEncodingSize::Vector{Int64})
  fac = acqData.encodingSize[1:2] ./ newEncodingSize[1:2]
  numContr = numContrasts(acqData)
  numSl = numSlices(acqData)
  numReps = numRepititions(acqData)
  idx = Vector{Vector{Int64}}(undef,numContr)

  for i=1:numContr
    tr = trajectory(acqData,i)
    nodes = fac .* kspaceNodes(tr)

    # filter out nodes with magnitude > 0.5
    idxX = findall(x->(x>=-0.5)&&(x<0.5), nodes[1,:])
    idxY = findall(x->(x>=-0.5)&&(x<0.5), nodes[2,:])
    idx[i] = intersect(idxX,idxY, acqData.subsampleIndices[i])

    tr.nodes = nodes[:,idx[i]]
    times = readoutTimes(tr)
    tr.times = times[idx[i]]
    acqData.subsampleIndices[i] = acqData.subsampleIndices[i][idx[i]]
  end

  # find relevant kspace data
  kdata2 = Array{Matrix{ComplexF64}}(undef,numContr,numSl,numReps)
  for rep=1:numReps
    for slice=1:numSl
      for echo=1:numContr
        kdata2[echo,slice,rep] = 1.0/prod(fac) * acqData.kdata[echo,slice,rep][idx[echo],:]
      end
    end
  end
  acqData.kdata = kdata2

  # adjust subsample Indices
  for echo=1:numContr
    acqData.subsampleIndices[echo] = collect(1:length(idx[echo]))
  end

  return acqData
end

"""
    convert3dTo2d(acqData::AcquisitionData)

convert the 3d encoded AcquisitionData `acqData` to the equivalent 2d AcquisitionData.
"""
function convert3dTo2d(acqData::AcquisitionData)
  numContr = numContrasts(acqData)
  numChan = numChannels(acqData)
  numSl = numSlices(trajectory(acqData,1))
  numReps = numRepititions(acqData)
  # check if all trajectories are cartesian
  for i=1:numContr
    if !isCartesian(trajectory(acqData,i))
      @error "conversion to 2d is not supported for non-cartesian data"
    end
  end

  # create 2d trajectories along phase encoding directions
  tr2d = Vector{Trajectory}(undef,numContr)
  for i=1:numContr
    tr3d = trajectory(acqData,i)
    # 1. arg (numProfiles=>y), 2. arg (numSamp=>x)
    tr2d[i] = CartesianTrajectory(numSlices(tr3d),numProfiles(tr3d),TE=echoTime(tr3d),AQ=acqTimePerProfile(tr3d))
  end

  # convert k-space data and place it in the appropriate array structure
  numSamp = numSamplingPerProfile(trajectory(acqData,1)) # assume the same number of samples for all contrasts
  kdata2d = Array{Matrix{ComplexF64}}(undef,numContr,numSamp, numReps)
  for i=1:numContr
    tr = trajectory(acqData,i)
    numProf = div( size(acqData.kdata[i,1,1],1), numSamp ) #numProfiles(tr)
    kdata_i = zeros(ComplexF64, numSamp, numProf, numChan, numReps)
    #convert
    F = 1/sqrt(numSamp)*FFTOp(ComplexF64, (numSamp,))
    for r=1:numReps
      for p=1:numProf # including slices
        for c=1:numChan
          # p_idx = (s-1)*numProf+p
          kdata_i[:,p,c,r] .= adjoint(F) * acqData.kdata[i,1,r][(p-1)*numSamp+1:p*numSamp,c]
        end
      end
    end
    # place kdata in transformed Array structure
    for r=1:numReps
      for j=1:numSamp
        kdata2d[i,j,r] = kdata_i[j,:,:,r] # numSl/numSamp*kdata_i[j,:,:,r]
      end
    end
  end

  # adapt subsampleIndices
  subsampleIndices2d = Vector{Vector{Int64}}(undef, numContr)
  for i=1:numContr
    idx = div.( acqData.subsampleIndices[i] .- 1, numSamp) .+ 1
    subsampleIndices2d[i] = sort(unique(idx))
  end

  return AcquisitionData(acqData.sequenceInfo, tr2d, kdata2d, subsampleIndices2d, acqData.encodingSize, acqData.fov)
end

hann(x) = 0.5*(1-cos(2*pi*(x-0.5)))

function apodization!(acqData::AcquisitionData)
    numContr = numContrasts(acqData)
    numSl = numSlices(acqData)
    numReps = numRepititions(acqData)

    for rep=1:numReps
      for slice=1:numSl
        for echo=1:numContr
          tr = trajectory(acqData,echo)
          nodes = kspaceNodes(tr)
          for k=1:size(acqData.kdata[echo,slice,rep],1)
            weight = hann(nodes[1,k])*hann(nodes[2,k])+0.5
            acqData.kdata[echo,slice,rep][k,:] .*= weight
          end
        end
      end
    end

    return acqData
end
