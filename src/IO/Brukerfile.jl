include("Jcampdx.jl")

export BrukerFile, studyName, studyNumber, experimentName, experimentNumber,
       scannerFacility, scannerOperator, scannerName, acqStartTime, acqNumFrames, acqNumAverages


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
  if !b.acqpRead && ( parameter=="NA" || parameter[1:3] == "ACQ" )
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

# study parameters
studyName(b::BrukerFile) = string(experimentSubject(b),"_",
                                  latin1toutf8(b["VisuStudyId"]),"_",
                                  b["VisuStudyNumber"])
studyNumber(b::BrukerFile) = parse(Int64,b["VisuStudyNumber"])
experimentName(b::BrukerFile) = latin1toutf8(b["ACQ_scan_name"])
experimentNumber(b::BrukerFile) = parse(Int64,b["VisuExperimentNumber"])

# scanner parameters
scannerFacility(b::BrukerFile) = latin1toutf8(b["ACQ_institution"])
scannerOperator(b::BrukerFile) = latin1toutf8(b["ACQ_operator"])
scannerName(b::BrukerFile) = b["ACQ_station"]

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
##$PVM_EncAvailReceivers=1

acqNumCoils(b::BrukerFile) = parse(Int,b["PVM_EncNReceivers"])
acqNumEchos(b::BrukerFile) = parse(Int,b["ACQ_n_echo_images"])
acqPhaseFactor(b::BrukerFile) = parse(Int,b["ACQ_phase_factor"])
acqRareFactor(b::BrukerFile) = parse(Int,b["ACQ_phase_factor"])
acqSpatialSize1(b::BrukerFile) = parse(Int,b["ACQ_spatial_size_1"])
acqNumRepetitions(b::BrukerFile) = parse(Int,b["NR"])
acqObjOrder(b::BrukerFile) = parse.(Int,b["ACQ_obj_order"])

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
    @error "Data type curr"
  end
  return Int8
end

# interface of MRIFile
function rawdata(b::BrukerFile)
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

  I_ = reshape(permutedims(I, (1,3,5,6,2,4,7) ), N[1], :,
                              numEncSteps2, numEchos, numSlices, numRep)
  encSteps1 = pvmEncSteps1(b)
  idxE1 = collect(1:length(encSteps1))
  permE1 = invpermute!(idxE1, encSteps1.-minimum(encSteps1).+1)

  encSteps2 = pvmEncSteps2(b)
  idxE2 = collect(1:length(encSteps2))
  permE2 = invpermute!(idxE2, encSteps2.-minimum(encSteps2).+1)

  objOrd = acqObjOrder(b)
  idxO = collect(1:length(objOrd))
  permO = invpermute!(idxO, objOrd.-minimum(objOrd).+1)

  I_ = I_[:,permE1,permE2,:,permO,:]

  return vec( map(ComplexF64, I_) )
end


function trajectory(b::BrukerFile)
  N = acqSize(b)
  if length(N) == 2
    return trajectory("Cartesian", N[1], N[2])
  elseif length(N) == 3
    return trajectory("Cartesian3D", N[1], N[2]; numSlices=N[3])
  end
  #return nothing
end


function sequence(f::BrukerFile)
  # TODO
end

function acquisitionData(b::BrukerFile)
  return AcquisitionData(trajectory(b), rawdata(b),
                          numCoils=acqNumCoils(b),
                          numEchoes=acqNumEchos(b),
                          numSlices=acqNumSlices(b),
                          encodingSize=acqSize(b))
end

##### Reco
function recoData(f::BrukerFile)
  recoFilename = joinpath(f.path,"pdata", "1", "2dseq")
  N = recoSize(f)

  #if f["RECO_wordtype",1] != "_16BIT_SGN_INT"
  #  @error "Not yet implemented!"
  #end

  I = open(recoFilename,"r") do fd
    read!(fd,Array{Int16,3}(undef,1,prod(N),1))
  end
  return map(Float32,I)
end

recoFov(f::BrukerFile) = push!(parse.(Float64,f["RECO_fov",1])./1000,
                                parse(Float64,f["ACQ_slice_sepn"][1])./1000)
recoFovCenter(f::BrukerFile) = zeros(3)
recoSize(f::BrukerFile) = push!(parse.(Int,f["RECO_size",1]),
                                parse(Int,f["RecoObjectsPerRepetition",1]))
#recoOrder(f::BrukerFile) = f["/reconstruction/order"]
#recoPositions(f::BrukerFile) = f["/reconstruction/positions"]
