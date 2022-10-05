export AcquisitionData, kData, kDataCart, kdataSingleSlice, convertUndersampledData,
       changeEncodingSize2D, convert3dTo2d, samplingDensity,
       numContrasts, numChannels, numSlices, numRepetitions,
       encodingSize, fieldOfView, multiCoilData

"""
struct describing MRI acquisition data.

# Fields
* `sequenceInfo::Dict{Symbol,Any}`          - additional information on the pulse sequence
* `traj::Vector{Trajectory}`                - trajectories for each echo/contrast
* `kdata::Array{Matrix{Complex{<:AbstractFloat}},3}`      - each matrix contains data for one trajectory
                                              (1. dim k-space nodes, 2. dim coils)
                                              the outer dims describe:
                                              1. dim echoes, 2. dim slices, 3. dim repetitions
* `subsampleIndices::Vector{Array{Int64}}`  - indices sampled for each echo/contrast
* `encodingSize::Vector{Int64}`             - size of the underlying image matrix
* `fov::Vector{Float64}`                    - field of view in m
"""
mutable struct AcquisitionData{T <: AbstractFloat}
  sequenceInfo::Dict{Symbol,Any}
  traj::Vector{Trajectory{T}}
  kdata::Array{Matrix{Complex{T}},3}
  subsampleIndices::Vector{Vector{Int64}}
  encodingSize::Vector{Int64}
  fov::Vector{Float64}
end

fieldOfView(acq::AcquisitionData) = acq.fov
encodingSize(acq::AcquisitionData) = acq.encodingSize

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
    numRepetitions(acqData::AcquisitionData)

returns the number of repetitions in acqData
"""
numRepetitions(acqData::AcquisitionData) = size(acqData.kdata,3)

"""
    AcquisitionData(tr::K,kdata::Array{Matrix{Complex{T}},3}; seqInfo=Dict{Symbol,Any}()
                        , idx=nothing, encodingSize=Int64[0,0,0], fov=Float64[0,0,0]
                        , kargs...) where {K <: Union{Trajectory,Vector{Trajectory}},T <: AbstractFloat}

constructor for `AcquisitionData`

# Arguments
* `tr <: Union{Trajectory,Vector{Trajectory}}` - trajectories
* `kdata::Array{Matrix{Complex{<:AbstractFloat}},3}` - k-space data

the other fields of `AcquisitionData` can be passed as keyword arguments.
"""
function AcquisitionData(tr::K, kdata::Array{Matrix{Complex{T}},3}
                        ; seqInfo=Dict{Symbol,Any}()
                        , idx=nothing
                        , encodingSize=Int64[0,0,0]
                        , fov=Float64[0,0,0]
                        , kargs...) where {T <: AbstractFloat, K <: Union{Trajectory,Vector{Trajectory{T}}}}
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

"""
    AcquistionData(kdata::Array{T,6})

Returns an AcquisitionData struct created from a zero-filled kspace with 6 dimensions :
    [x,y,z,channels,echoes,repetitions]

Different undersampling patterns can be handled **only** along the echo dimension.

This methods assumes the data to correspond to a single volume
encoded with a 3d Cartesian sequence.
"""
function AcquisitionData(kspace::Array{Complex{T},6}) where T
    sx,sy,sz,nCh,nEchos,nReps = size(kspace)

    if sz == 1
        tr = MRIBase.CartesianTrajectory(T,sx,sy,TE=T(0),AQ=T(0))
    else
        tr = MRIBase.CartesianTrajectory3D(T,sx,sy,numSlices=sz,TE=T(0),AQ=T(0))
    end

    kdata = [reshape(kspace,:,nCh) for i=1:nEchos, j=1:1, k=1:nReps]
    acq = AcquisitionData(tr,kdata)

    acq.encodingSize = [sx,sy,sz]


    for echo in 1:nEchos
        I = findall(x->x!=0,abs.(kspace[:,:,:,:,echo,1]))
        subsampleInd = LinearIndices((sx,sy,sz))[I]

        acq.subsampleIndices[echo]=subsampleInd
        for rep in 1:nReps
            acq.kdata[echo,1,rep] = acq.kdata[echo,1,rep][subsampleInd,:]
        end
    end
    return acq
end

#function Images.pixelspacing(acqData::AcquisitionData)
#  return [1.0,1.0,1.0]*Unitful.mm
#  #return fov./encodingSize*Unitful.mm  #TODO: all needs to be properly initialized
#end

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
    kDataCart(acqData::AcquisitionData)

Returns the cartesian k-space contained in `acqData` for all `echo`, `coil`, `slice` and `rep`(etition)
with dimension :
[x,y,z*slices,channels,echoes,repetitions]
"""
function kDataCart(acqData::AcquisitionData)
nx, ny, nz = acqData.encodingSize[1:3]
numChan, numSl = numChannels(acqData), numSlices(acqData)

if nz > 1 && numSl > 1
    @warn "Multi slab 3D acquisitions are concatenate along the 3rd dimension"
end

numEcho = length(acqData.traj)
numRep = numRepetitions(acqData)
kdata = zeros(ComplexF64, nx * ny * nz,numSl,numChan, numEcho,numRep)
for rep = 1:numRep
    for sl =1:numSl
        for echo = 1:numEcho
            for coil = 1:numChan
                kdata[acqData.subsampleIndices[echo],sl,coil,echo,rep] .= kData(acqData, echo, coil, sl,rep=rep)
            end
        end
    end
end
kdata = reshape(kdata, nx, ny, nz*numSl, numChan,numEcho,numRep)
return kdata
end

"""
    multiEchoData(acqData::AcquisitionData, coil::Int64, slice::Int64;rep::Int64=1)

returns the k-space contained in `acqData` for all echoes and given `coil`, `slice` and `rep`(etition).
"""
function multiEchoData(acqData::AcquisitionData{T}, coil::Int64, slice::Int64;rep::Int64=1) where T
  kdata = T[]
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
function multiCoilMultiEchoData(acqData::AcquisitionData{T},slice::Int64;rep=1) where T
  kdata = Complex{T}[]
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
  if ndims(tr)==2 || numSl==1
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
function samplingDensity(acqData::AcquisitionData{T},shape::Tuple) where T
  numContr = numContrasts(acqData)
  weights = Array{Vector{Complex{T}}}(undef,numContr)
  for echo=1:numContr
    tr = trajectory(acqData,echo)
    if isCartesian(tr)
      nodes = kspaceNodes(tr)[:,acqData.subsampleIndices[echo]]
      weights[echo] = [1.0/sqrt(prod(shape)) for node=1:size(nodes,2)]
    else
      nodes = kspaceNodes(tr)
      plan = plan_nfft(nodes, shape, m=2, σ=2)
      weights[echo] = sqrt.(sdc(plan, iters=10))
    end

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

  if !isCartesian(trajectory(acqData, 1))
    @error "changing the encoding size does not yet support non-cartesian sampling. Argument unchanged."
    return acqData
  end

  dest = deepcopy(acqData)
  changeEncodingSize2D!(dest,newEncodingSize)
end

"""
    changeEncodingSize2D!(acqData::AcquisitionData,newEncodingSize::Vector{Int64})

does the same thing as `changeEncodingSize2D` but acts in-place on `acqData`.
"""
function changeEncodingSize2D!(acqData::AcquisitionData{T},newEncodingSize::Vector{Int64}) where T
  fac = acqData.encodingSize[1:2] ./ newEncodingSize[1:2]
  numContr = numContrasts(acqData)
  numSl = numSlices(acqData)
  numReps = numRepetitions(acqData)
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

    # adjust subsampleIndices if the trajectory is cartesian (this matters for espirit)
    if tr.cartesian
      acqData.subsampleIndices[i] = findIndices(newEncodingSize,acqData.encodingSize)[idx[i]]
    else
      acqData.subsampleIndices[i] = acqData.subsampleIndices[i][idx[i]]
    end
  end

  # find relevant kspace data
  kdata2 = Array{Matrix{Complex{T}}}(undef,numContr,numSl,numReps)
  for rep=1:numReps
    for slice=1:numSl
      for echo=1:numContr
        kdata2[echo,slice,rep] = 1.0/prod(fac) * acqData.kdata[echo,slice,rep][idx[echo],:]
      end
    end
  end
  acqData.kdata = kdata2

  # # adjust subsample Indices
  # for echo=1:numContr
  #   acqData.subsampleIndices[echo] = collect(1:length(idx[echo]))
  # end

  # formally change the encoding size in the acqData 
  acqData.encodingSize = newEncodingSize

  return acqData
end

"""
    convert3dTo2d(acqData::AcquisitionData)

convert the 3d encoded AcquisitionData `acqData` to the equivalent 2d AcquisitionData.
"""
function convert3dTo2d(acqData::AcquisitionData{T}) where T
  numContr = numContrasts(acqData)
  numChan = numChannels(acqData)
  numSl = numSlices(trajectory(acqData,1))
  numReps = numRepetitions(acqData)
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
  kdata2d = Array{Matrix{Complex{T}}}(undef,numContr,numSamp, numReps)
  for i=1:numContr
    tr = trajectory(acqData,i)
    numProf = div( size(acqData.kdata[i,1,1],1), numSamp ) #numProfiles(tr)
    kdata_i = zeros(T, numSamp, numProf, numChan, numReps)
    #convert
    F = 1/sqrt(numSamp)*FFTOp(T, (numSamp,))
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

#= Is this still in use? And why pre-weight the k-space data?

function NFFT.apodization!(acqData::AcquisitionData)
    numContr = numContrasts(acqData)
    numSl = numSlices(acqData)
    numReps = numRepetitions(acqData)

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
end=#

"""
    findIndices()

puts the subsamplingIndices in the right range after converting encoding size for Cartesian trajectories if the encoding size increases. Leaves non-cartesian alone...
"""
function findIndices(newEnc::Vector{Int64},oldEnc::Vector{Int64})

    shiftX = floor((newEnc[1]-oldEnc[1])÷2)
    shiftY = floor((newEnc[2]-oldEnc[2])÷2)

    oldCart = CartesianIndices(((shiftX+1):(shiftX+oldEnc[1]),(shiftY+1):(shiftY+oldEnc[2])))
    newCart = CartesianIndices((1:newEnc[1],1:newEnc[2]))
  
  if newEnc[1] > oldEnc[1] || newEnc[2] > oldEnc[2]
    
    return LinearIndices(newCart)[oldCart][:]

  else

    return LinearIndices(oldCart)[:]
  end

end