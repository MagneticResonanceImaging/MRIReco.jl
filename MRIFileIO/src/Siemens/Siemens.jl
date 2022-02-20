# TODO


struct SiemensFile <: MRIFile
  filename::String
end

function MRIBase.RawAcquisitionData(f::SiemensFile)
  error("Not yet implemented!")
end