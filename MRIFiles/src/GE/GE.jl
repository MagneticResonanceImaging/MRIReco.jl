# TODO

struct GEFile <: MRIFile
  filename::String
end

function MRIBase.RawAcquisitionData(f::GEFile)
  error("Not yet implemented!")
end