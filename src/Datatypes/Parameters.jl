export Limit, MeasurementDependency, CoilDescription, minimalHeader

function GeneralParameters(str::String)
  return GeneralParameters( parse_string(str) )
end

function addToDict!(params::Dict, e::XMLElement, paramName::String, ::Type{T},
                    dictParamName=paramName) where T
  u = get_elements_by_tagname(e, paramName)
  if !isempty(u)
    if length(u) == 1
      params[dictParamName] = customParse(T,u[1])
    else
      tmp = T[]
      for l=1:length(u)
        push!(tmp, customParse(T,u[l]))
      end
      params[dictParamName] = tmp
    end
  end
end

struct Limit
  minimum::Int
  maximum::Int
  center::Int
end

struct MeasurementDependency
  dependencyType::String
  measurementID::String
end

struct CoilDescription
  number::Int
  name::String
end

customParse(::Type{T}, node) where T = parse(T,content(node))
customParse(::Type{String}, node) where T = content(node)
function customParse(::Type{Vector{T}}, node) where T
  x = content(get_elements_by_tagname(node, "x")[1])
  y = content(get_elements_by_tagname(node, "y")[1])
  z = content(get_elements_by_tagname(node, "z")[1])
  return [parse(T,x),parse(T,y),parse(T,z)]
end
function customParse(::Type{MeasurementDependency}, node)
  a = content(get_elements_by_tagname(node, "dependencyType")[1])
  b = content(get_elements_by_tagname(node, "measurementID")[1])
  return MeasurementDependency(a, b)
end
function customParse(::Type{CoilDescription}, node)
  a = content(get_elements_by_tagname(node, "coilNumber")[1])
  b = content(get_elements_by_tagname(node, "coilName")[1])
  return CoilDescription(parse(Int,a), b)
end
function customParse(::Type{Limit}, node)
  a = content(get_elements_by_tagname(node, "minimum")[1])
  b = content(get_elements_by_tagname(node, "maximum")[1])
  c = content(get_elements_by_tagname(node, "center")[1])
  return Limit(parse(Int,a), parse(Int,b), parse(Int,c))
end

function GeneralParameters(xdoc::XMLDocument)
    params = Dict{String,Any}()

    # SubjectInformation
    e = get_elements_by_tagname(LightXML.root(xdoc),"subjectInformation")
    if !isempty(e)
      addToDict!(params, e[1], "patientName", String)
      addToDict!(params, e[1], "patientWeight_kg", Float64)
      addToDict!(params, e[1], "patientID", String)
      addToDict!(params, e[1], "patientBirthdate", String)
      addToDict!(params, e[1], "patientGender", String)
    end

    # StudyInformation
    e = get_elements_by_tagname(LightXML.root(xdoc),"studyInformation")
    if !isempty(e)
      addToDict!(params, e[1], "studyDate", String)
      addToDict!(params, e[1], "studyDate", String)
      addToDict!(params, e[1], "studyID", String)
      addToDict!(params, e[1], "accessionNumber", Int)
      addToDict!(params, e[1], "referringPhysicianName", String)
      addToDict!(params, e[1], "studyDescription", String)
      addToDict!(params, e[1], "studyInstanceUID", String)
    end

    # MeasurementInformation
    e = get_elements_by_tagname(LightXML.root(xdoc),"measurementInformation")
    if !isempty(e)
      addToDict!(params, e[1], "measurementID", String)
      addToDict!(params, e[1], "seriesDate", String)
      addToDict!(params, e[1], "seriesTime", String)
      addToDict!(params, e[1], "patientPosition", String) # non-optional?
      addToDict!(params, e[1], "initialSeriesNumber", Int)
      addToDict!(params, e[1], "protocolName", String)
      addToDict!(params, e[1], "seriesDescription", String)
      addToDict!(params, e[1], "seriesInstanceUIDRoot", String)
      addToDict!(params, e[1], "frameOfReferenceUID", String)

      addToDict!(params, e[1], "measurementDependency", MeasurementDependency)
      addToDict!(params, e[1], "referencedImageSequence", String)
    end

    # AcquisitionSystemInformation
    e = get_elements_by_tagname(LightXML.root(xdoc),"acquisitionSystemInformation")
    if !isempty(e)
      addToDict!(params, e[1], "systemVendor", String)
      addToDict!(params, e[1], "systemModel", String)
      addToDict!(params, e[1], "systemFieldStrength_T", Float64)
      addToDict!(params, e[1], "relativeReceiverNoiseBandwidth", Float64)
      addToDict!(params, e[1], "receiverChannels", Int)
      addToDict!(params, e[1], "institutionName", String)
      addToDict!(params, e[1], "stationName", String)
      addToDict!(params, e[1], "coilLabel", CoilDescription)
    end

    # ExperimentalConditions
    e = get_elements_by_tagname(LightXML.root(xdoc),"experimentalConditions")
    if !isempty(e)
      addToDict!(params, e[1], "H1resonanceFrequency_Hz", Int)
    end

    # Encoding
    e = get_elements_by_tagname(LightXML.root(xdoc),"encoding")
    if !isempty(e)
      if !isempty(e[1]["encodedSpace"])
        params["encodedSize"] = customParse(Vector{Int}, e[1]["encodedSpace"][1]["matrixSize"][1])
        params["encodedFOV"] = customParse(Vector{Float64}, e[1]["encodedSpace"][1]["fieldOfView_mm"][1])
      end
      if !isempty(e[1]["reconSpace"])
        params["reconSize"] = customParse(Vector{Int}, e[1]["reconSpace"][1]["matrixSize"][1])
        params["reconFOV"] = customParse(Vector{Float64}, e[1]["reconSpace"][1]["fieldOfView_mm"][1])
      end

      d = e[1]["encodingLimits"]
      if !isempty(d)
        addToDict!(params, d[1], "kspace_encoding_step_0", Limit, "enc_lim_kspace_encoding_step_0")
        addToDict!(params, d[1], "kspace_encoding_step_1", Limit, "enc_lim_kspace_encoding_step_1")
        addToDict!(params, d[1], "kspace_encoding_step_2", Limit, "enc_lim_kspace_encoding_step_2")
        addToDict!(params, d[1], "average", Limit, "enc_lim_average")
        addToDict!(params, d[1], "slice", Limit, "enc_lim_slice")
        addToDict!(params, d[1], "contrast", Limit, "enc_lim_contrast")
        addToDict!(params, d[1], "phase", Limit, "enc_lim_phase")
        addToDict!(params, d[1], "repetition", Limit, "enc_lim_repetition")
        addToDict!(params, d[1], "set", Limit, "enc_lim_set")
        addToDict!(params, d[1], "segment", Limit, "enc_lim_segment")
      end

      addToDict!(params, e[1], "trajectory", String)

      d = e[1]["trajectoryDescription"]
      if !isempty(d)
        y = Dict{String,Any}()
        y["identifier"] = content(d[1]["identifier"][1])
        for q in d[1]["userParameterLong"]
          y[content(q["name"][1])] = parse(Int,content(q["value"][1]))
        end
        for q in d[1]["userParameterDouble"]
          y[content(q["name"][1])] = parse(Float64,content(q["value"][1]))
        end
        for q in d[1]["userParameterString"]
          y[content(q["name"][1])] = content(q["value"][1])
        end
        if !isempty(d[1]["comment"])
          y["comment"] = content(d[1]["comment"][1])
        end
        params["trajectoryDescription"] = y
      end

      d = e[1]["parallelImaging"]
      if !isempty(d)
        a = parse(Int,content(d[1]["accelerationFactor"][1]["kspace_encoding_step_1"][1]))
        b = parse(Int,content(d[1]["accelerationFactor"][1]["kspace_encoding_step_2"][1]))
        params["accelerationFactor"] = (a,b)

        addToDict!(params, d[1], "calibrationMode", String)
        addToDict!(params, d[1], "interleavingDimension", String)
      end
      addToDict!(params, e[1], "echoTrainLength", Int)
    end

    # Support wrongly created files (which are downloaded in the tests):
    e = get_elements_by_tagname(LightXML.root(xdoc),"parallelImaging")
    if !isempty(e)
      a = parse(Int,content(e[1]["accelerationFactor"][1]["kspace_encoding_step_1"][1]))
      b = parse(Int,content(e[1]["accelerationFactor"][1]["kspace_encoding_step_2"][1]))
      params["accelerationFactor"] = (a,b)

      addToDict!(params, e[1], "calibrationMode", String)
      addToDict!(params, e[1], "interleavingDimension", String)
    end

    # SequenceParameters
    e = get_elements_by_tagname(LightXML.root(xdoc),"sequenceParameters")
    if !isempty(e)
      addToDict!(params, e[1], "TR", Float64)
      addToDict!(params, e[1], "TE", Float64)
      addToDict!(params, e[1], "TI", Float64)
      addToDict!(params, e[1], "flipAngle_deg", Float64)
      addToDict!(params, e[1], "sequence_type", String)
      addToDict!(params, e[1], "echo_spacing", Float64)
    end

    # waveformInformation
    e = get_elements_by_tagname(LightXML.root(xdoc),"waveformInformation")
    if !isempty(e)
      addToDict!(params, e[1], "waveformName", Float64)
      addToDict!(params, e[1], "waveformType", Float64)

      if !isempty(e[1]["userParameters"])
          d = e[1]["userParameters"]
          y = Dict{String,Any}()
          for q in d["userParameterLong"]
            y[content(q["name"][1])] = parse(Int,content(q["value"][1]))
          end
          for q in d["userParameterDouble"]
            y[content(q["name"][1])] = parse(Float64,content(q["value"][1]))
          end
          for q in d["userParameterString"]
            y[content(q["name"][1])] = content(q["value"][1])
          end
          params["waveformUserParameters"] = y
      end
    end

    # UserParameters
    e = get_elements_by_tagname(LightXML.root(xdoc),"userParameters")
    if !isempty(e)
        d = e[1]
        y = Dict{String,Any}()
        for q in d["userParameterLong"]
          y[content(q["name"][1])] = parse(Int,content(q["value"][1]))
        end
        for q in d["userParameterDouble"]
          y[content(q["name"][1])] = parse(Float64,content(q["value"][1]))
        end
        for q in d["userParameterString"]
          y[content(q["name"][1])] = content(q["value"][1])
        end
        params["userParameters"] = y
    end

    return params
end


function generateGroup(params, paramVec, node, groupName, prefix="")
  if any( [haskey(params,prefix*n) for n in paramVec] )
    xs = new_child(node, groupName)
    for p in paramVec
      p_ = prefix*p
      if haskey(params,p_)
        insertNode(xs, p, params[p_])
      end
    end
  end
end


function insertNode(xs, paramName::String, param::T) where T
  xsp = new_child(xs, paramName)
  add_text(xsp, string(param))
end

function insertNode(xs, paramName::String, param::Vector{T}) where T <: Number
  xsp = new_child(xs, paramName)
  insertNode(xsp, "x", param[1])
  insertNode(xsp, "y", param[2])
  insertNode(xsp, "z", param[3])
end

function insertNode(xs, paramName::String, param::Vector{T}) where T
  for p in param
    insertNode(xs, paramName, p)
  end
end

function insertNode(xs, paramName::String, param::Tuple{T,T}) where T <: Integer
  xsp = new_child(xs, "accelerationFactor")
  insertNode(xsp, "kspace_encoding_step_1", param[1])
  insertNode(xsp, "kspace_encoding_step_2", param[2])
end

function insertNode(xs, paramName::String, param::MeasurementDependency)
  xsp = new_child(xs, paramName)
  insertNode(xsp, "dependencyType", param.dependencyType)
  insertNode(xsp, "measurementID", param.measurementID)
end

function insertNode(xs, paramName::String, param::CoilDescription)
  xsp = new_child(xs, paramName)
  insertNode(xsp, "coilNumber", param.number)
  insertNode(xsp, "coilName", param.name)
end

function insertNode(xs, paramName::String, param::Limit)
  xsp = new_child(xs, paramName)
  insertNode(xsp, "minimum", param.minimum)
  insertNode(xsp, "maximum", param.maximum)
  insertNode(xsp, "center", param.center)
end

function insertNode(xs, paramName::String, param_::Dict{String,Any})
  xsp = new_child(xs, paramName)
  param = deepcopy(param_)
  if haskey(param, "comment")
    insertNode(xsp, "comment", param["comment"])
    delete!(param,"comment")
  end
  if haskey(param, "identifier")
    insertNode(xsp, "identifier", param["identifier"])
    delete!(param,"identifier")
  end
  for (key,val) in param
    if typeof(val) <: Integer
      xspp = new_child(xsp, "userParameterLong")
      insertNode(xspp, "name", key)
      insertNode(xspp, "value", val)
    elseif typeof(val) <: AbstractFloat
      xspp = new_child(xsp, "userParameterDouble")
      insertNode(xspp, "name", key)
      insertNode(xspp, "value", val)
    elseif typeof(val) <: AbstractString
      xspp = new_child(xsp, "userParameterString")
      insertNode(xspp, "name", key)
      insertNode(xspp, "value", val)
    end
  end
end


function GeneralParametersToXML(params::Dict{String,Any})

  xdoc = XMLDocument()
  xroot = create_root(xdoc, "ismrmrdHeader")
  set_attributes(xroot, Dict{String,String}("xmlns"=>"http://www.ismrm.org/ISMRMRD",
                         "xmlns:xsi"=>"http://www.w3.org/2001/XMLSchema-instance",
                         "xmlns:xs"=>"http://www.w3.org/2001/XMLSchema",
                         "xsi:schemaLocation"=>"http://www.ismrm.org/ISMRMRD ismrmrd.xsd"))

  p = ["patientName","patientWeight_kg", "patientID", "patientBirthdate", "patientGender"]
  generateGroup(params, p, xroot, "subjectInformation")

  p = ["studyDate", "studyDate", "studyID", "accessionNumber","referringPhysicianName",
       "studyDescription", "studyInstanceUID"]
  generateGroup(params, p, xroot, "studyInformation")

  p = ["measurementID", "seriesDate", "seriesTime", "patientPosition", "initialSeriesNumber",
       "protocolName", "seriesDescription", "seriesInstanceUIDRoot", "frameOfReferenceUID",
       "measurementDependency", "referencedImageSequence"]
  generateGroup(params, p, xroot, "measurementInformation")

  p = ["systemVendor", "systemModel", "systemFieldStrength_T",
       "relativeReceiverNoiseBandwidth", "receiverChannels",
       "institutionName", "stationName", "coilLabel"]
  generateGroup(params, p, xroot, "acquisitionSystemInformation")


  p = ["H1resonanceFrequency_Hz"]
  generateGroup(params, p, xroot, "experimentalConditions")

  pEnc = ["encodedSize", "encodedFOV"]
  pRecon = ["reconSize", "reconFOV"]
  pEncLim = ["kspace_encoding_step_0", "kspace_encoding_step_1",
             "kspace_encoding_step_2", "average", "slice",
             "contrast", "phase", "repetition", "set",
             "segment"]

  pElse = ["trajectory", "trajectoryDescription", "echoTrainLength"]

  pPI = ["accelerationFactor", "calibrationMode", "interleavingDimension"]

  p = union(pEnc, pRecon, pEncLim, pElse, pPI)

  if any( [haskey(params,n) for n in p] )
    xs = new_child(xroot, "encoding")

    if haskey(params, "encodedSize") || haskey(params, "encodedFOV")
      xsp = new_child(xs, "encodedSpace")
      insertNode(xsp, "matrixSize", params["encodedSize"])
      insertNode(xsp, "fieldOfView_mm", params["encodedFOV"])
    end

    if haskey(params, "reconSize") || haskey(params, "reconFOV")
      xsp = new_child(xs, "reconSpace")
      insertNode(xsp, "matrixSize", params["reconSize"])
      insertNode(xsp, "fieldOfView_mm", params["reconFOV"])
    end

    generateGroup(params, pEncLim, xs, "encodingLimits", "enc_lim_")

    for n in pElse
      if haskey(params, n)
        insertNode(xs, n, params[n])
      end
    end

    generateGroup(params, pPI, xs, "parallelImaging")
  end

  p = ["TR", "TE", "TI", "flipAngle_deg", "sequence_type", "echo_spacing"]
  generateGroup(params, p, xroot, "sequenceParameters")

  p = ["waveformName", "waveformType", "waveformUserParameters"]
  generateGroup(params, p, xroot, "waveformInformation")

  p = ["userParameters"]
  generateGroup(params, p, xroot, "userParameters")

  return string(xdoc)
end

function minimalHeader(encodingSize::NTuple{3,Int},fov::NTuple{3,AbstractFloat};f_res::Integer=1,tr_name::AbstractString="cartesian",numChannels::Int=1)
  params = Dict{String,Any}()
  params["H1resonanceFrequency_Hz"] = f_res
  params["encodedSize"] = collect(encodingSize)
  params["encodedFOV"] = collect(fov)
  params["trajectory"] = tr_name
  params["receiverChannels"] = numChannels

  return params
end
