export BrukerFile, recoData

include("Jcampdx.jl")


function latin1toutf8(str::AbstractString)
  buff = Char[]
  for c in Vector{UInt8}(str)
    push!(buff,c)
  end
  string(buff...)
end

function latin1toutf8(str::Nothing)
  error("Can not convert noting to UTF8!")
end

mutable struct BrukerFile <: MRIFile
  path::String
  params::JcampdxFile
  paramsProc::JcampdxFile
  methodRead
  acqpRead
  visupars_globalRead
  recoRead
  methrecoRead
  visuparsRead
  mpiParRead
  maxEntriesAcqp
end

function BrukerFile(path::String; maxEntriesAcqp=2000)
  params = JcampdxFile()
  paramsProc = JcampdxFile()

  return BrukerFile(path, params, paramsProc, false, false, false,
               false, false, false, false, maxEntriesAcqp)
end

function getindex(b::BrukerFile, parameter)#::String
  if !b.acqpRead && ( parameter=="NA" || parameter=="NR" || parameter=="NI" || parameter=="NSLICES" ||
                      parameter=="SW" || parameter[1:2] == "GO" || parameter[1:3] == "ACQ" )
    acqppath = joinpath(b.path, "acqp")
    read(b.params, acqppath, maxEntries=b.maxEntriesAcqp)
    b.acqpRead = true
  elseif !b.visupars_globalRead && length(parameter) >= 4 &&
         parameter[1:4] == "Visu"
    visupath = joinpath(b.path, "visu_pars")
    if isfile(visupath)
      read(b.params, visupath)
      b.visupars_globalRead = true
    end
  elseif !b.mpiParRead && length(parameter) >= 6 &&
         parameter[1:6] == "CONFIG"
    mpiParPath = joinpath(b.path, "mpi.par")
    if isfile(mpiParPath)
      read(b.params, mpiParPath)
      b.mpiParRead = true
    end
  else !b.methodRead 
      methodpath = joinpath(b.path, "method")
      read(b.params, methodpath)
      b.methodRead = true
  end

  if haskey(b.params, parameter)
    return b.params[parameter]
  else
    return ""
  end
end

# for pdata

function getindex(b::BrukerFile, parameter, procno::Int64)#::String
  if !b.recoRead && lowercase( parameter[1:4] ) == "reco"
    recopath = joinpath(b.path, "pdata", string(procno), "reco")
    read(b.paramsProc, recopath, maxEntries=13)
    b.recoRead = true
  elseif !b.methrecoRead && parameter[1:3] == "PVM"
    methrecopath = joinpath(b.path, "pdata", string(procno), "methreco")
    read(b.paramsProc, methrecopath)
    b.methrecoRead = true
  elseif !b.visuparsRead && parameter[1:4] == "Visu"
    visuparspath = joinpath(b.path, "pdata", string(procno), "visu_pars")
    if isfile(visuparspath)
      read(b.paramsProc, visuparspath)
      b.visuparsRead = true
    end
  end

  return b.paramsProc[parameter]
end

function Base.show(io::IO, b::BrukerFile)
  print(io, "BrukerFile: ", b.path)
end

# acquisition parameters
function acqStartTime(b::BrukerFile)
  m = match(r"<(.+)\+",b["ACQ_time"])
  timeString = replace(m.captures[1],"," => ".")
  return DateTime( timeString )
end
function acqNumFrames(b::BrukerFile)
  M = Int64(b["ACQ_jobs"][1][8])
  return div(M,acqNumPeriodsPerFrame(b))
end

acqNumAverages(b::BrukerFile) = parse(Int,b["NA"])
acqNumSlices(b::BrukerFile) = parse(Int,b["NSLICES"])
acqNumInterleaves(b::BrukerFile) = parse(Int,b["NI"])

function acqSize(b::BrukerFile)
  N = parse.(Int,b["ACQ_size"])
  N[1] = div(N[1],2)
  return N
end
function acqFov(b::BrukerFile)
  N = parse.(Float64,b["ACQ_fov"])./100
  if length(N) == 3
    return N
  else
    return push!(N, parse(Float64,b["ACQ_slice_sepn"][1])./100)
  end
end
##$PVM_EncAvailReceivers=1

acqNumCoils(b::BrukerFile) = parse(Int,b["PVM_EncNReceivers"])
acqNumEchos(b::BrukerFile) = parse(Int,b["ACQ_n_echo_images"])
acqPhaseFactor(b::BrukerFile) = parse(Int,b["ACQ_phase_factor"])
acqRareFactor(b::BrukerFile) = parse(Int,b["ACQ_phase_factor"])
acqSpatialSize1(b::BrukerFile) = parse(Int,b["ACQ_spatial_size_1"])
acqNumRepetitions(b::BrukerFile) = parse(Int,b["NR"])
acqObjOrder(b::BrukerFile) = parse.(Int,b["ACQ_obj_order"])
acqReadOffset(b::BrukerFile) = parse.(Float64,b["ACQ_read_offset"])
acqPhase1Offset(b::BrukerFile) = parse.(Float64,b["ACQ_phase1_offset"])
acqPhase2Offset(b::BrukerFile) = parse.(Float64,b["ACQ_phase2_offset"])
acqGradMatrix(b::BrukerFile) = parse.(Float64,b["ACQ_grad_matrix"])
acqSliceOffset(b::BrukerFile) = parse.(Float64,b["ACQ_slice_offset"])
acqFlipAngle(b::BrukerFile) = parse(Float64,b["ACQ_flip_angle"])
acqProtocolName(b::BrukerFile) = b["ACQ_protocol_name"]
acqInterEchoTime(b::BrukerFile) = parse(Float64,b["ACQ_inter_echo_time"][1])
acqEchoTime(b::BrukerFile) = parse(Float64,b["ACQ_echo_time"][1])
acqRepetitionTime(b::BrukerFile) = parse(Float64,b["ACQ_repetition_time"][1])

Base.ndims(b::BrukerFile) = parse(Int, b["ACQ_dim"])

pvmEncNReceivers(b::BrukerFile) =parse.(Int,b["PVM_EncNReceivers"])
pvmEncAvailReceivers(b::BrukerFile) =parse.(Int,b["PVM_EncAvailReceivers"])
pvmEncSteps1(b::BrukerFile) = parse.(Int,b["PVM_EncSteps1"])
pvmEncSteps2(b::BrukerFile) = parse.(Int,b["PVM_EncSteps2"])
pvmEncValues1(b::BrukerFile) = parse.(Float32,b["PVM_EncSteps1"])
pvmEncValues2(b::BrukerFile) = parse.(Float32,b["PVM_EncSteps2"])
pvmMatrix(b::BrukerFile) = parse.(Int,b["PVM_Matrix"])
pvmSpiralMode(b::BrukerFile) = b["PVM_SpiralMode"]
pvmSpiralNbOfInterleaves(b::BrukerFile) = parse.(Int,b["PVM_SpiralNbOfInterleaves"])
pvmSpiralNbOfGradientPoints(b::BrukerFile) = parse.(Int,b["PVM_SpiralNbOfGradientPoints"])
pvmSpiralNbOfAcqPoints(b::BrukerFile) = parse.(Int,b["PVM_SpiralNbOfAcqPoints"])
pvmSpiralEchoTime(b::BrukerFile) = parse.(Int,b["PVM_SpiralEchoTime"])
pvmSpiralAcqDwellTime(b::BrukerFile) = parse(Float32, b["PVM_SpiralAcqDwellTime"])
pvmSpiralAcquisitionTime(b::BrukerFile) = parse(Float32, b["PVM_SpiralAcquisitionTime"])
pvmSpiralNavSize(b::BrukerFile) = parse.(Int,b["PVM_SpiralNavSize"])
pvmSpiralPreSize(b::BrukerFile) = parse.(Int,b["PVM_SpiralPreSize"])
pvmSpiralSize(b::BrukerFile) = parse.(Int,b["PVM_SpiralSize"])
pvmSpiralPostSize(b::BrukerFile) = parse.(Int,b["PVM_SpiralPostSize"])
pvmSpiralPointsPerRotation(b::BrukerFile) = parse.(Int,b["PVM_SpiralPointsPerRotation"])
pvmSpiralShape1(b::BrukerFile) = parse.(Float32,b["PVM_SpiralShape1"])
pvmSpiralShape2(b::BrukerFile) = parse.(Float32,b["PVM_SpiralShape2"])

pvmTrajInterleaves(b::BrukerFile) = parse(Int,b["PVM_TrajInterleaves"])
pvmTrajSamples(b::BrukerFile) = parse(Int,b["PVM_TrajSamples"])
pvmTrajResultSize(b::BrukerFile) = parse(Int,b["PVM_TrajResultSize"])

pvmTrajKx(b::BrukerFile) = parse.(Float32,b["PVM_TrajKx"])
pvmTrajKy(b::BrukerFile) = parse.(Float32,b["PVM_TrajKy"])


"""
brukerParams(b::BrukerFile)

Generate a standard parameter dictionary from Bruker acquistion 
Some fields needs to be filled outside lik

"""
function brukerParams(b::BrukerFile)
  params = Dict{String,Any}()
  params["trajectory"] = "cartesian"

  params["encodedSize"] = pvmMatrix(b)
  F = acqFov(b)
  params["encodedFOV"] = F
  params["receiverChannels"] = pvmEncNReceivers(b)
  params["H1resonanceFrequency_Hz"] = round(Int, parse(Float64,b["SW"])*1000000)
  params["studyID"] = b["VisuStudyId"]
  #params["studyDescription"] = b["ACQ_scan_name"]
  #params["studyInstanceUID"] =
  params["referringPhysicianName"] = latin1toutf8(b["ACQ_operator"])

  params["patientName"] = b["VisuSubjectName"]

  params["measurementID"] = parse(Int64,b["VisuExperimentNumber"])
  params["seriesDescription"] = b["ACQ_scan_name"]

  params["institutionName"] = latin1toutf8(b["ACQ_institution"])
  params["stationName"] = b["ACQ_station"]
  params["systemVendor"] = "Bruker"

  params["TR"] = acqRepetitionTime(b)
  params["TE"] = acqEchoTime(b)
  #params["TI"] = ???
  params["flipAngle_deg"] = acqFlipAngle(b)
  params["sequence_type"] = acqProtocolName(b)
  params["echo_spacing"] = acqInterEchoTime(b)

  return params
end

function acqDataType(b::BrukerFile)
  format = b["GO_raw_data_format"]
  if format == "GO_32BIT_SGN_INT"
    return Int32
  elseif format == "GO_16BIT_SGN_INT"
    return Int16
  elseif format == "GO_32BIT_FLOAT"
  else
    @error "Data type unknown: $(format)"
  end
  return Int8
end

function MRIBase.RawAcquisitionData(b::BrukerFile)
  if isfile(joinpath(b.path, "fid"))
    protocolName = acqProtocolName(b)
    @info "Bruker protocol name : $protocolName"
    if acqProtocolName(b) == "UTE3D"
      return RawAcquisitionData_3DUTE(b)
    elseif acqProtocolName(b) == "SPIRAL"
      return RawAcquisitionDataRawDataSpiral(b)
    else
      return RawAcquisitionDataFid(b) #standard cartesian cases
    end
  end
end

include("BrukerSequence.jl")

function RawAcquisitionDataRawDataSpiral(b::BrukerFile)
  # Have not a way to read out this. Not sure if its true
  T = Complex{Int32}
  
  filename = joinpath(b.path, "rawdata.job0")

  profileLength = pvmTrajSamples(b)
  phaseFactor = acqPhaseFactor(b)
  numSlices = acqNumSlices(b)
  numEchos = acqNumEchos(b)
  numRep = acqNumRepetitions(b)
  numProfiles = pvmTrajInterleaves(b)
  numChannels = acqNumCoils(b)

  I = open(filename,"r") do fd
    read!(fd,Array{T,7}(undef, profileLength,
                                   numChannels,
                                   numEchos,
                                   phaseFactor,
                                   numSlices,
                                   numProfiles,
                                   numRep))
  end

  # TODO the following of course is just for debugging
  return I
end

function RawAcquisitionDataFid(b::BrukerFile)
    T = Complex{acqDataType(b)}

    filename = joinpath(b.path, "fid")

    N = acqSize(b)
    # The data is padded in case it is not a multiple of 1024 
    # For multi-channel acquisition data at concatenate then padded to a multiple of 1024 bytes
    numChannel = pvmEncNReceivers(b)
    profileLength = Int((ceil(N[1]*numChannel*sizeof(T)/1024))*1024/sizeof(T))
    numAvailableChannel = pvmEncAvailReceivers(b)
    phaseFactor = acqPhaseFactor(b)
    numSlices = acqNumSlices(b)
    numEchos = acqNumEchos(b)
    numEncSteps2 = length(N) == 3 ? N[3] : 1
    numRep = acqNumRepetitions(b)

    I = open(filename,"r") do fd
      read!(fd,Array{T,7}(undef, profileLength,
                                     numEchos,
                                     phaseFactor,
                                     numSlices,
                                     div(N[2], phaseFactor),
                                     numEncSteps2,
                                     numRep))[1:N[1]*numChannel,:,:,:,:,:,:]
    end

    encSteps1 = pvmEncSteps1(b)
    encSteps1 = encSteps1.-minimum(encSteps1)

    encSteps2 = pvmEncSteps2(b)
    encSteps2 = encSteps2.-minimum(encSteps2)

    objOrd = acqObjOrder(b)
    objOrd = objOrd.-minimum(objOrd)

    gradMatrix = acqGradMatrix(b)

    offset1 = acqReadOffset(b)
    offset2 = acqPhase1Offset(b)
    offset3 = ndims(b) == 2 ? acqSliceOffset(b) : acqPhase2Offset(b)

    profiles = Profile[]
    for nR = 1:numRep
      for nEnc2 = 1:numEncSteps2
        for nPhase2 = 1:div(N[2], phaseFactor)
          for nSl = 1:numSlices
            for nPhase1 = 1:phaseFactor
              for nEcho=1:numEchos
                  counter = EncodingCounters(kspace_encode_step_1=encSteps1[nPhase1+phaseFactor*(nPhase2-1)],
                                             kspace_encode_step_2=encSteps2[nEnc2],
                                             average=0,
                                             slice=objOrd[nSl],
                                             contrast=nEcho-1,
                                             phase=0,
                                             repetition=nR-1,
                                             set=0,
                                             segment=0 )

                  G = gradMatrix[:,:,nSl]
                  read_dir = (G[1,1],G[2,1],G[3,1])
                  phase_dir = (G[1,2],G[2,2],G[3,2])
                  slice_dir = (G[1,3],G[2,3],G[3,3])

                  # Not sure if the following is correct...
                  pos = offset1[nSl]*G[:,1] +
                        offset2[nSl]*G[:,2] +
                        offset3[nSl]*G[:,3]

                  position = (pos[1], pos[2], pos[3])

                  head = AcquisitionHeader(number_of_samples=N[1], idx=counter,
                                           read_dir=read_dir, phase_dir=phase_dir,
                                           slice_dir=slice_dir, position=position,
                                           center_sample=div(N[1],2),
                                           available_channels = numAvailableChannel, #TODO
                                           active_channels = numChannel)
                  traj = Matrix{Float32}(undef,0,0)
                  dat = map(T, reshape(I[:,nEcho,nPhase1,nSl,nPhase2,nEnc2,nR],:,numChannel))
                  push!(profiles, Profile(head,traj,dat) )
              end
            end
          end
        end
      end
    end

    params = brukerParams(b)
    params["trajectory"] = "cartesian"
    N = acqSize(b)
    if length(N) < 3
      N_ = ones(Int,3)
      N_[1:length(N)] .= N
      N = N_
    end
    params["encodedSize"] = N

    return RawAcquisitionData(params, profiles)
end


##### Reco
function recoData(f::BrukerFile)
  recoFilename = joinpath(f.path,"pdata", "1", "2dseq")
  N = recoSize(f)

  #if f["RECO_wordtype",1] != "_16BIT_SGN_INT"
  #  @error "Not yet implemented!"
  #end

  I = open(recoFilename,"r") do fd
    read!(fd,Array{Int16,length(N)}(undef,N...))
  end
  return map(Float32,I)
end

recoFov(f::BrukerFile) = push!(parse.(Float64,f["RECO_fov",1])./100,
                                parse(Float64,f["ACQ_slice_sepn"][1])./100)
recoFovCenter(f::BrukerFile) = zeros(3)
recoSize(f::BrukerFile) = push!(parse.(Int,f["RECO_size",1]),
                                parse(Int,f["RecoObjectsPerRepetition",1]))
