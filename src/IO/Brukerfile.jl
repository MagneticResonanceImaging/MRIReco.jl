include("Jcampdx.jl")

export BrukerFile, recoData


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
  if !b.acqpRead && ( parameter=="NA" || parameter=="NR" || parameter=="NI" ||
                      parameter=="SW" || parameter[1:2] == "GO" || parameter[1:3] == "ACQ" )
    acqppath = joinpath(b.path, "acqp")
    read(b.params, acqppath, maxEntries=b.maxEntriesAcqp)
    b.acqpRead = true
  elseif !b.methodRead && length(parameter) >= 3 &&
         (parameter[1:3] == "PVM" || parameter[1:3] == "MPI")
    methodpath = joinpath(b.path, "method")
    read(b.params, methodpath)
    b.methodRead = true
  elseif !b.visupars_globalRead && length(parameter) >= 4 &&
         parameter[1:4] == "Visu"
    visupath = joinpath(b.path, "visu_pars")
    if isfile(visupath)
      keylist = ["VisuStudyId","VisuStudyNumber","VisuExperimentNumber","VisuSubjectName"]
      read(b.params, visupath,keylist)
      b.visupars_globalRead = true
    end
  elseif !b.mpiParRead && length(parameter) >= 6 &&
         parameter[1:6] == "CONFIG"
    mpiParPath = joinpath(b.path, "mpi.par")
    if isfile(mpiParPath)
      read(b.params, mpiParPath)
      b.mpiParRead = true
    end
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

pvmEncSteps1(b::BrukerFile) = parse.(Int,b["PVM_EncSteps1"])
pvmEncSteps2(b::BrukerFile) = parse.(Int,b["PVM_EncSteps2"])
pvmEncValues1(b::BrukerFile) = parse.(Float32,b["PVM_EncSteps1"])
pvmEncValues2(b::BrukerFile) = parse.(Float32,b["PVM_EncSteps2"])
pvmMatrix(b::BrukerFile) = parse.(Int,b["PVM_Matrix"])


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

function RawAcquisitionData(b::BrukerFile)
    dtype = Complex{acqDataType(b)}

    filename = joinpath(b.path, "fid")

    N = acqSize(b)
    # The data is padded in case it is not a multiple of 1024
    profileLength = div((div(N[1]*sizeof(dtype),1024)+1)*1024,sizeof(dtype))
    phaseFactor = acqPhaseFactor(b)
    numSlices = acqNumSlices(b)
    numEchos = acqNumEchos(b)
    numEncSteps2 = length(N) == 3 ? N[3] : 1
    numRep = acqNumRepetitions(b)

    I = open(filename,"r") do fd
      read!(fd,Array{dtype,7}(undef, profileLength,
                                     numEchos,
                                     phaseFactor,
                                     numSlices,
                                     div(N[2], phaseFactor),
                                     numEncSteps2,
                                     numRep))[1:N[1],:,:,:,:,:,:]
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
                                           available_channels = 1, #TODO
                                           active_channels = 1)
                  traj = Matrix{Float32}(undef,0,0)
                  dat = map(ComplexF64, reshape(I[:,nEcho,nPhase1,nSl,nPhase2,nEnc2,nR],:,1))
                  push!(profiles, Profile(head,traj,dat) )
              end
            end
          end
        end
      end
    end

    params = Dict{String,Any}()
    params["trajectory"] = "cartesian"
    N = acqSize(b)
    if length(N) < 3
      N_ = ones(Int,3)
      N_[1:length(N)] .= N
      N = N_
    end
    params["encodedSize"] = N
    F = acqFov(b)
    params["encodedFOV"] = F
    params["receiverChannels"] = 1
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
    read!(fd,Array{Int16,3}(undef,N...))
  end
  return map(Float32,I)
end

recoFov(f::BrukerFile) = push!(parse.(Float64,f["RECO_fov",1])./100,
                                parse(Float64,f["ACQ_slice_sepn"][1])./100)
recoFovCenter(f::BrukerFile) = zeros(3)
recoSize(f::BrukerFile) = push!(parse.(Int,f["RECO_size",1]),
                                parse(Int,f["RecoObjectsPerRepetition",1]))
