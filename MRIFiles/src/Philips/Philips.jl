# TODO

struct PhilipsFile <: MRIFile
  filename::String
end

function MRIBase.RawAcquisitionData(f::PhilipsFile)
  error("Not yet implemented!")
end