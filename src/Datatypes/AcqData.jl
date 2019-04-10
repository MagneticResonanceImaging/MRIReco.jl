export AcquisitionData, kData, kdataSingleSlice, convertUndersampledData,
       weightData!, weightedData, unweightData!, unweightDataSquared!, changeEncodingSize2D

"""
 contains all relevant info for the aquired data.
 kdata is stored as a vector.
 samplePointer contains the start of each datasample (for one echo, one slice and one coil).
 datasamples are ordered like a 4d object with:
 1. dim : kspace nodes
 2. dim : echo times
 3. dim : coils
 4. slices / repetitions
"""
mutable struct AcquisitionData
  sequenceInfo::Dict{Symbol,Any}
  traj::Vector{Trajectory}
  kdata::Vector{ComplexF64}
  numEchoes::Int64
  numCoils::Int64
  numSlices::Int64
  samplePointer::Vector{Int64}
  subsampleIndices::Array{Int64}
  encodingSize::Vector{Int64}
  fov::Vector{Float64}
end

##############
# constructors
##############
function AcquisitionData(tr::Trajectory, kdata::Vector{ComplexF64}
                        ; seqInfo=Dict{Symbol,Any}(),kargs...)
  return AcquisitionData([tr], kdata; seqInfo=seqInfo, kargs...)
end

function AcquisitionData(tr::Vector{Trajectory},kdata::Vector{ComplexF64}
                        ; seqInfo=Dict{Symbol,Any}()
                        , numCoils=1
                        , numEchoes=1
                        , numSlices=1
                        , idx=nothing
                        , encodingSize=Int64[]
                        , fov = Float64[]
                        , kargs...)
  encodingDims = dims(tr[1])
  if encodingDims==2
    return AcquisitionData2d(tr,kdata,seqInfo=seqInfo,numCoils=numCoils,numEchoes=numEchoes,
              numSlices=numSlices,idx=idx,encodingSize=encodingSize,fov=fov)
  elseif encodingDims==3
    return AcquisitionData3d(tr,kdata,seqInfo=seqInfo,numCoils=numCoils,numEchoes=numEchoes,
                numSlices=numSlices,idx=idx,encodingSize=encodingSize,fov=fov)
  else
    @error "encoding $(encodingDims) not yet supported"
  end
  return AcquisitionData2d(seqInfo,tr,kdata,numCoils=numCoils,numEchoes=numEchoes,numSlices=numSlices,idx=idx,encodingSize=encodingSize,fov=fov)
end

function AcquisitionData2d(tr::Vector{Trajectory}, kdata::Vector{ComplexF64};
                        seqInfo=Dict{Symbol,Any}(), numCoils=1, numEchoes=1, numSlices=1, idx=nothing,
                        encodingSize=Int64[], fov=Float64[])
  idx==nothing && ( idx = collect(1:length(kdata)) )
  numSamplesPerShot = div(length(kdata),(numEchoes*numCoils*numSlices))
  samplePointer = collect(1:numSamplesPerShot:length(kdata)-numSamplesPerShot+1)
  return AcquisitionData(seqInfo, tr, kdata, numEchoes, numCoils, numSlices,
                         samplePointer, idx, encodingSize, fov)
end

#
# aquisition data for single sequences
#
function AcquisitionData3d(tr::Vector{Trajectory}, kdata::Vector{ComplexF64};
                        seqInfo=Dict{Symbol,Any}(), numCoils=1, numEchoes=1, numSlices=1, idx=nothing,
                        encodingSize=Int64[], fov=Float64[])
  idx==nothing && ( idx = collect(1:length(kdata)) )
  numSamplesPerShot = div(length(kdata),(numEchoes*numCoils))
  samplePointer = collect(1:numSamplesPerShot:length(kdata)-numSamplesPerShot+1)
  return AcquisitionData(seqInfo, tr, kdata, numEchoes, numCoils, numSlices,
                         samplePointer, idx, encodingSize, fov)
end

function Images.pixelspacing(acqData::AcquisitionData)
  return [1.0,1.0,1.0]*Unitful.mm
  #return fov./encodingSize*Unitful.mm  #TODO: all needs to be properly initialized
end

############
# trajectory
############
trajectory(acqData::AcquisitionData,i::Int64=1) = acqData.traj[i]

######################
# getting k-space data
######################

#
# return kdata for a given coil, echo and image slice
#
function kData(acqData::AcquisitionData, echo::Int64=1, coil::Int64=1, slice::Int64=1)
  setNumber = ((slice-1)*acqData.numCoils+coil-1)*acqData.numEchoes+echo
  numSets = acqData.numEchoes*acqData.numCoils*acqData.numSlices
  if length(acqData.samplePointer) >1
    if setNumber < numSets
      return acqData.kdata[acqData.samplePointer[setNumber] : acqData.samplePointer[setNumber+1]-1]
    else
      return acqData.kdata[acqData.samplePointer[setNumber] : end]
    end
  else
    return acqData.kdata
  end
end

#
# return multi-echo-data for a given coil and image slice
#
function multiEchoData(acqData::AcquisitionData, coil::Int64, slice::Int64)
  numNodesPerSet = div(length(acqData.kdata),acqData.numSlices*acqData.numCoils)
  kdata = reshape(acqData.kdata, numNodesPerSet, acqData.numCoils, acqData.numSlices)
  return kdata[:,coil,slice]
end

#
# return multi-coil-data for a given echo and image slice
#
function multiCoilData(acqData::AcquisitionData, echo::Int64, slice::Int64)
  numNodesPerSet = div(length(acqData.kdata),acqData.numCoils*acqData.numEchoes*acqData.numSlices)
  kdata = reshape(acqData.kdata, numNodesPerSet, acqData.numEchoes, acqData.numCoils, acqData.numSlices)
  return vec(kdata[:,echo,:,slice])

end

#
# return multi-coil-multi-echo-data for a given slice (for all coils and echoes)
#
function multiCoilMultiEchoData(acqData::AcquisitionData,slice::Int64)
  numNodesPerSlice = div(length(acqData.kdata),acqData.numSlices)
  return acqData.kdata[(slice-1)*numNodesPerSlice+1:slice*numNodesPerSlice]
end


######################################
# utilities to convert and edit acqData
######################################

"""
dims(tr::Trajectory) =
 convert undersampled AcquisitionData, where only profiles contained in
 acqData.subsampleIndices are sampled,
 into a format where trajectories only contain the sampled profiles
"""
function convertUndersampledData(acqData::AcquisitionData)

  acqDataSub = deepcopy(acqData)

  # get number of nodes and reset idx
  numNodes = size(acqData.subsampleIndices,1)
  acqDataSub.subsampleIndices = collect(1:length(acqData.kdata))

  # replace trajectories by Undersampled Trajectories
  for i = 1:acqData.numEchoes
    tr = trajectory(acqDataSub,i)
    # assume that coils and slices experience the same trajectory
    tr.nodes = tr.nodes[:,acqData.subsampleIndices[:,i,1,1]]
  end

  return acqDataSub
end

#
# weight kdata with the sampling density of the trajectories used
#
function weightedData(acqData::AcquisitionData, shape::Tuple)
  res = deepcopy(acqData)
  weightData!(res, shape)
  return res
end

function weightData!(acqData::AcquisitionData, shape::NTuple{2,Int64})
  ft = [ NFFTOp(shape,trajectory(acqData,i)) for i=1:acqData.numEchoes]
  for i = 1:acqData.numSlices
    for j = 1:acqData.numCoils
      for k = 1:acqData.numEchoes
        idx = ((i-1)*acqData.numCoils+j-1)*acqData.numEchoes+k

        if k!=acqData.numEchoes || j!=acqData.numCoils || i!=acqData.numSlices
          acqData.kdata[acqData.samplePointer[idx] : acqData.samplePointer[idx+1]-1] .*= sqrt.(ft[k].density)
        else
          acqData.kdata[acqData.samplePointer[idx] : end] .*= sqrt.(ft[k].density)
        end
      end
    end
  end
  return acqData
end

function weightData!(acqData::AcquisitionData, shape::NTuple{3,Int64})
  ft = [ NFFTOp(shape,trajectory(acqData,i)) for i=1:acqData.numEchoes]
  for j = 1:acqData.numCoils
    for k = 1:acqData.numEchoes
      idx = (j-1)*acqData.numEchoes+k
      if k!=acqData.numEchoes || j!=acqData.numCoils
        acqData.kdata[acqData.samplePointer[idx] : acqData.samplePointer[idx+1]-1] .*= sqrt.(ft[k].density)
      else
        acqData.kdata[acqData.samplePointer[idx] : end] .*= sqrt.(ft[k].density)
      end
    end
  end
  return acqData
end

function unweightData!(acqData::AcquisitionData, shape::NTuple{2,Int64})
  ft = [ NFFTOp(shape,trajectory(acqData,i)) for i=1:acqData.numEchoes]
  for i = 1:acqData.numSlices
    for j = 1:acqData.numCoils
      for k = 1:acqData.numEchoes
        idx = ((i-1)*acqData.numCoils+j-1)*acqData.numEchoes+k
        if k!=acqData.numEchoes || j!=acqData.numCoils || i!=acqData.numSlices
          acqData.kdata[acqData.samplePointer[idx] : acqData.samplePointer[idx+1]-1] ./= sqrt.(ft[k].density)
        else
          acqData.kdata[acqData.samplePointer[idx] : end] ./= sqrt.(ft[k].density)
        end
      end
    end
  end
  return acqData
end

function unweightData!(acqData::AcquisitionData, shape::NTuple{3,Int64})
  ft = [ NFFTOp(shape,trajectory(acqData,i)) for i=1:acqData.numEchoes]
  for j = 1:acqData.numCoils
    for k = 1:acqData.numEchoes
      idx = (j-1)*acqData.numEchoes+k
      if k!=acqData.numEchoes || j!=acqData.numCoils
        acqData.kdata[acqData.samplePointer[idx] : acqData.samplePointer[idx+1]-1] ./= sqrt.(ft[k].density)
      else
        acqData.kdata[acqData.samplePointer[idx] : end] ./= sqrt.(ft[k].density)
      end
    end
  end
  return acqData
end

function unweightDataSquared!(acqData::AcquisitionData, shape::Tuple)
  unweightData!(acqData, shape)
  unweightData!(acqData, shape)
  return acqData
end

#########################################################################
# convert acqData for a reconstruction with a encodingSize (resolution)
#########################################################################
function changeEncodingSize2D(acqData::AcquisitionData,newEncodingSize::Vector{Int64})
  dest = deepcopy(acqData)
  changeEncodingSize2D!(dest,newEncodingSize)
end

function changeEncodingSize2D!(acqData::AcquisitionData,newEncodingSize::Vector{Int64})
  fac = acqData.encodingSize ./ newEncodingSize
  idx = Vector{Vector{Int64}}(undef,acqData.numEchoes)

  for i=1:acqData.numEchoes
    tr = trajectory(acqData,i)
    nodes = fac .* kspaceNodes(tr)

    # filter out nodes with magnitude > 0.5
    idxX = findall(x->(x>=-0.5)&&(x<0.5), nodes[1,:])
    idxY = findall(x->(x>=-0.5)&&(x<0.5), nodes[2,:])
    idx[i] = intersect(idxX,idxY)

    tr.nodes = nodes[:,idx[i]]
    times = readoutTimes(tr)
    tr.times = times[idx[i]]
  end

  # find relevant kspace data
  kdata = reshape(acqData.kdata,:,acqData.numEchoes,acqData.numCoils, acqData.numSlices)
  kdata2 = zeros(ComplexF64,length(idx[1]),acqData.numEchoes,acqData.numCoils, acqData.numSlices)
  for i=1:acqData.numEchoes
    kdata2[:,i,:,:] = 1.0/prod(fac) * kdata[idx[i],i,:,:]
  end
  acqData.kdata = vec(kdata2)

  # update sample pointer
  numSamplesPerShot = div(length(kdata2),(acqData.numEchoes*acqData.numCoils*acqData.numSlices))
  acqData.samplePointer = collect(1:numSamplesPerShot:length(kdata2)-numSamplesPerShot+1)

  return acqData
end
