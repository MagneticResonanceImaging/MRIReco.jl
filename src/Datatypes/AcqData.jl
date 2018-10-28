export AcquisitionData, kData, kdataSingleSlice, convertUndersampledData,weightData!, weightedData

#=
 contains all relevant info for the aquired data.
 kdata is stored as a vector.
 samplePointer contains the start of each datasample (for one echo, one slice and one coil).
 datasamples are ordered like a 4d object with:
 1. dim : kspace nodes
 2. dim : echo times
 3. dim : coils
 4. slices
=#
mutable struct AcquisitionData{S<:AbstractSequence}
  seq::S
  kdata::Vector{ComplexF64}
  numEchoes::Int64
  numCoils::Int64
  numSlices::Int64
  samplePointer::Vector{Int64}
  idx::Array{Int64}
end

##############
# constructors
##############
#
# constructor function for single trajectories
#
function AcquisitionData(tr::AbstractTrajectory, kdata::Vector{ComplexF64};
                        numCoils=1, numEchoes=1, numSlices=1, idx=nothing)
  seq = SESequence(tr)
  return AcquisitionData(seq, kdata, numCoils=numCoils, numEchoes=numEchoes, numSlices=numSlices)
end

#
# aquisition data for single sequences
#
function AcquisitionData(seq::AbstractSequence, kdata::Vector{ComplexF64};
                        numCoils=1, numEchoes=1, numSlices=1, idx=nothing)
  idx==nothing && ( idx = collect(1:length(kdata)) )
  numSamplesPerShot = div(length(kdata),(numEchoes*numCoils*numSlices))
  samplePointer = collect(1:numSamplesPerShot:length(kdata)-numSamplesPerShot+1)
  return AcquisitionData(seq, kdata, numEchoes, numCoils, numSlices, samplePointer, idx)
end

function Images.pixelspacing(acqData::AcquisitionData)
  return [1.0,1.0,1.0]*Unitful.mm
end

######################
# getting k-space data
######################

#
# return kdata for a given coil, echo and image slice
#
function kData(aqData::AcquisitionData, echo::Int64=1, coil::Int64=1, slice::Int64=1)
  setNumber = ((slice-1)*aqData.numCoils+coil-1)*aqData.numEchoes+echo
  numSets = aqData.numEchoes*aqData.numCoils*aqData.numSlices
  if length(aqData.samplePointer) >1
    if setNumber < numSets
      return aqData.kdata[aqData.samplePointer[setNumber] : aqData.samplePointer[setNumber+1]-1]
    else
      return aqData.kdata[aqData.samplePointer[setNumber] : end]
    end
  else
    return aqData.kdata
  end
end

#
# return multi-echo-data for a given coil and image slice
#
function multiEchoData(aqData::AcquisitionData, coil::Int64, slice::Int64)
  numNodesPerSet = div(length(aqData.kdata),aqData.numSlices*aqData.numCoils)
  kdata = reshape(aqData.kdata, numNodesPerSet, aqData.numCoils, aqData.numSlices)
  return kdata[:,coil,slice]
end

#
# return multi-coil-data for a given echo and image slice
#
function multiCoilData(aqData::AcquisitionData, echo::Int64, slice::Int64)
  numNodesPerSet = div(length(aqData.kdata),aqData.numCoils*aqData.numEchoes*aqData.numSlices)
  kdata = reshape(aqData.kdata, numNodesPerSet, aqData.numEchoes, aqData.numCoils, aqData.numSlices)
  return vec(kdata[:,echo,:,slice])

end

#
# return multi-coil-multi-echo-data for a given slice (for all coils and echoes)
#
function multiCoilMultiEchoData(aqData::AcquisitionData,slice::Int64)
  numNodesPerSlice = div(length(aqData.kdata),aqData.numSlices)
  return aqData.kdata[(slice-1)*numNodesPerSlice+1:slice*numNodesPerSlice]
end


######################################
# utilities to convert and edit aqData
######################################

#
# convert undersampled AcquisitionData, where only profiles contained in aqData.idx are sampled,
# into a format where trajectories only contain the sampled profiles
#
function convertUndersampledData(aqData::AcquisitionData)

  aqDataCopy = deepcopy(aqData)

  # get number of nodes and reset idx
  numNodes = size(aqData.idx,1)
  aqDataCopy.idx = collect(1:length(aqData.kdata))

  # replace trajectories by Undersampled Trajectories
  for i = 1:aqData.numEchoes
    tr = trajectory(aqData.seq,i)
    # assume that coils and slices experience the same trajectory
    setTrajectory!(aqDataCopy.seq, UndersampledTrajectory(tr, aqData.idx[:,i,1,1]), i)
  end

  return aqDataCopy
end

#
# weight kdata with the sampling density of the trajectories used
#
function weightedData(aqData::AcquisitionData, shape::Tuple)
  res = deepcopy(aqData)
  weightData!(res, shape)
  return res
end

function weightData!(aqData::AcquisitionData, shape::Tuple)
  ft = [ NFFTOp(shape,trajectory(aqData.seq,i)) for i=1:aqData.numEchoes]
  for i = 1:aqData.numSlices
    for j = 1:aqData.numCoils
      for k = 1:aqData.numEchoes
        idx = ((i-1)*aqData.numCoils+j-1)*aqData.numEchoes+k
        if k!=aqData.numEchoes || j!=aqData.numCoils || i!=aqData.numSlices
          aqData.kdata[aqData.samplePointer[idx] : aqData.samplePointer[idx+1]-1] .*= sqrt.(ft[k].density)
        else
          aqData.kdata[aqData.samplePointer[idx] : end] .*= sqrt.(ft[k].density)
        end
      end
    end
  end
end
