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


# interface of MRIFile

function rawdata(f::BrukerFile, repetition=1)
  return map(ComplexF64,vec(permutedims(f.data[:,:,findIndices(f,repetition)],(1,3,2))))
end

function acquisitionData(f::BrukerFile)
  return AcquisitionData(trajectory(f), rawdata(f),
                          numCoils=numChannels(f), numEchoes=1, numSlices=1)
end

#=
function measData(b::BrukerFile, frames=1:acqNumFrames(b), periods=1:acqNumPeriodsPerFrame(b),
                  receivers=1:rxNumChannels(b))

  dataFilename = joinpath(b.path,"rawdata.job0")
  dType = acqNumAverages(b) == 1 ? Int16 : Int32

  s = open(dataFilename)

  if numSubPeriods(b) == 1
    raw = Mmap.mmap(s, Array{dType,4},
             (rxNumSamplingPoints(b),rxNumChannels(b),acqNumPeriodsPerFrame(b),acqNumFrames(b)))
  else
    raw = Mmap.mmap(s, Array{dType,5},
             (rxNumSamplingPoints(b),numSubPeriods(b),rxNumChannels(b),acqNumPeriodsPerFrame(b),acqNumFrames(b)))
    raw = dropdims(sum(raw,dims=2),dims=2)
  end
  data = raw[:,receivers,periods,frames]
  close(s)

  return reshape(data, rxNumSamplingPoints(b), length(receivers),length(periods),length(frames))
end
=#

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
