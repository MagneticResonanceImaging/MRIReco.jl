# TODO

include("SiemensRaw.jl")

struct SiemensFile <: MRIFile
  filename::String
end

function MRIBase.RawAcquisitionData(f::SiemensFile)
  error("Not yet implemented!")
end

export readRawSiemens
function readRawSiemens(filename)

  # See https://github.com/ismrmrd/siemens_to_ismrmrd/blob/master/main.cpp

  open(filename) do fd
    buffer = Vector{UInt8}(undef, sizeof(MrParcRaidFileHeader))
    readbytes!(fd, buffer)
    ParcRaidHead =  reinterpret(MrParcRaidFileHeader, buffer)[1]

    @info ParcRaidHead.hdSize_
    @info ParcRaidHead.count_


    VBFILE = false

    if ParcRaidHead.hdSize_ > 32
        VBFILE = true

        # Rewind, we have no raid file header.
        # siemens_dat.seekg(0, std::ios::beg);
        seekstart(fd)

        ParcRaidHead.hdSize_ = ParcRaidHead.count_
        ParcRaidHead.count_ = 1
    elseif ParcRaidHead.hdSize_ != 0
        # This is a VB line data file
        error("Only VD line files with MrParcRaidFileHeader.hdSize_ == 0 (MR_PARC_RAID_ALLDATA) supported.")
    end

  end

end