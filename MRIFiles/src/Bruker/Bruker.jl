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

include("BrukerSequencePV6.jl")
include("BrukerSequence360.jl")

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
    visupaths = joinpath.(readdir(joinpath(b.path,"pdata"),join=true),"visu_pars")
    visupath = visupaths[isfile.(visupaths)]
    if !isempty(visupath)
      read(b.params, visupath[1])
      b.visupars_globalRead = true
    end
  elseif !b.mpiParRead && length(parameter) >= 6 &&
         parameter[1:6] == "CONFIG"
    mpiParPath = joinpath(b.path, "mpi.par")
    if isfile(mpiParPath)
      read(b.params, mpiParPath)
      b.mpiParRead = true
    end
  elseif !b.methodRead
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
  elseif b.visuparsRead != procno && parameter[1:4] == "Visu"
    visuparspath = joinpath(b.path, "pdata", string(procno), "visu_pars")
    if isfile(visuparspath)
      read(b.paramsProc, visuparspath)
      b.visuparsRead = procno
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
  N = parse.(Float64,b["ACQ_fov"]) .* 10 #cm to mm
  if length(N) == 3
    return N
  else
    return push!(N, parse(Float64,b["PVM_SPackArrSliceDistance"][1]) .* 10) #cm to mm
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
acqReadOffset(b::BrukerFile) = parse.(Float64,b["ACQ_read_offset"]) #mm
acqPhase1Offset(b::BrukerFile) = parse.(Float64,b["ACQ_phase1_offset"]) #mm
acqPhase2Offset(b::BrukerFile) = parse.(Float64,b["ACQ_phase2_offset"]) #mm
acqGradMatrix(b::BrukerFile) = parse.(Float64,b["ACQ_grad_matrix"])
acqSliceOffset(b::BrukerFile) = parse.(Float64,b["ACQ_slice_offset"]) #mm
acqFlipAngle(b::BrukerFile) = parse(Float64,b["ACQ_flip_angle"])
acqProtocolName(b::BrukerFile) = b["ACQ_protocol_name"]
acqInterEchoTime(b::BrukerFile) = parse(Float64,b["ACQ_inter_echo_time"][1])
acqEchoTime(b::BrukerFile) = parse(Float64,b["ACQ_echo_time"][1])
acqRepetitionTime(b::BrukerFile) = parse(Float64,b["ACQ_repetition_time"][1])

Base.ndims(b::BrukerFile) = parse(Int, b["ACQ_dim"])

pvmEffPhase2Offset(b::BrukerFile) =parse.(Float32,b["PVM_EffPhase2Offset"]) #mm
pvmEncNReceivers(b::BrukerFile) =parse.(Int,b["PVM_EncNReceivers"])
pvmEncAvailReceivers(b::BrukerFile) =parse.(Int,b["PVM_EncAvailReceivers"])
pvmEncSteps1(b::BrukerFile) = parse.(Int,b["PVM_EncSteps1"])
pvmEncSteps2(b::BrukerFile) = parse.(Int,b["PVM_EncSteps2"])
pvmEncValues1(b::BrukerFile) = parse.(Float32,b["PVM_EncSteps1"])
pvmEncValues2(b::BrukerFile) = parse.(Float32,b["PVM_EncSteps2"])
pvmMatrix(b::BrukerFile) = parse.(Int,b["PVM_Matrix"])
pvmEncMatrix(b::BrukerFile) = parse.(Int,b["PVM_EncMatrix"])
pvmEncPpi(b::BrukerFile) = parse.(Int,b["PVM_EncPpi"]) 
pvmEncPft(b::BrukerFile) = parse.(Float32,b["PVM_EncPft"])
pvmAntiAlias(b::BrukerFile) = parse.(Float32,b["PVM_AntiAlias"])
pvmEncZf(b::BrukerFile) = parse.(Float32,b["PVM_EncZf"])
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

visuExperimentNumber(b::BrukerFile) = isempty(b["VisuExperimentNumber"]) ? 0 : parse(Int64,b["VisuExperimentNumber"])
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
  params["H1resonanceFrequency_Hz"] = round(Int, parse(Float64,b["PVM_FrqWork"][1])*1000000)
  params["studyID"] = b["VisuStudyId"]
  #params["studyDescription"] = b["ACQ_scan_name"]
  #params["studyInstanceUID"] =
  params["referringPhysicianName"] = latin1toutf8(b["ACQ_operator"])

  params["patientName"] = b["VisuSubjectName"]

  params["measurementID"] = visuExperimentNumber(b)
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

function acqWordSize(b::BrukerFile)
  format = b["ACQ_word_size"]
  if format == "_32_BIT"
    return Int32
  elseif format == "_16_BIT"
    return Int16
  else
    @error "Data type unknown: $(format)"
  end
  return Int8
end

function MRIBase.RawAcquisitionData(b::BrukerFile)
  protocolName = acqProtocolName(b)
  @info "Bruker protocol name : $protocolName"

  if occursin("PV-360",b["ACQ_sw_version"]) ## PV360
    return RawAcquisitionDataFid_360(b) #standard cartesian cases
  else ## PV6
    if protocolName == "UTE3D"
      return RawAcquisitionData_3DUTE(b)
    elseif protocolName == "SPIRAL"
      return RawAcquisitionDataRawDataSpiral(b)
    else
      return RawAcquisitionDataFid(b) #standard cartesian cases
    end
  end
end


##### Reco
function recoData(f::BrukerFile,procno::Int=1)
  recoFilename = joinpath(f.path,"pdata", string(procno), "2dseq")
  if !isfile(recoFilename) @error "file : $recoFilename does not exist"; return end

  nFrame = parse.(Int64,f["VisuCoreFrameCount",procno])
  N = parse.(Int64,f["VisuCoreSize",procno])

  if(f["VisuCoreWordType",procno] == "_16BIT_SGN_INT")
    T = Int16
  elseif (f["VisuCoreWordType",procno] == "_32BIT_SGN_INT")
    T = Int32
  elseif (f["VisuCoreWordType",procno] == "_32BIT_FLOAT")
    T = Float32
  elseif (f["VisuCoreWordType",procno] == "_8BIT_UNSGN_INT")
    T = UInt8
  end

  I = open(recoFilename,"r") do fd
    read!(fd,Array{T,length(N)+1}(undef,N...,nFrame))
  end

  if(f["VisuCoreFrameType",procno] == ["REAL_IMAGE", "IMAGINARY_IMAGE"])
    @warn "1st half of image are the REAL part and 2nd half is the imaginary"
  end

  return map(Float32,I)
end

recoFov(f::BrukerFile) = push!(parse.(Float64,f["RECO_fov",1])./100,
                                parse(Float64,f["ACQ_slice_sepn"][1])./100)
recoFovCenter(f::BrukerFile) = zeros(3)
recoSize(f::BrukerFile) = push!(parse.(Int,f["RECO_size",1]),
                                parse(Int,f["RecoObjectsPerRepetition",1]))
