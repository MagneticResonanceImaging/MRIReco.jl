
const MDH_NUMBEROFEVALINFOMASK = 2
const MDH_NUMBEROFICEPROGRAMPARA_VB = 4
const MDH_NUMBEROFICEPROGRAMPARA_VD = 24
const MDH_FREEHDRPARA_VB = 4
const MDH_DMA_LENGTH_MASK = 0x01FFFFFF#L
const MDH_PACK_BIT_MASK = 0x02000000#L
const MDH_ENABLE_FLAGS_MASK = 0xFC000000#L
const MDH_SYNCDATA = 0x00000020#L

@enum PMU_Type begin
	END =  0x01FF0000
	ECG1 = 0x01010000
	ECG2 = 0x01020000
	ECG3 = 0x01030000
	ECG4 = 0x01040000
	PULS = 0x01050000
	RESP = 0x01060000
	EXT1 = 0x01070000
	EXT2 = 0x01080000
end

#=
enum class Trajectory {
  TRAJECTORY_CARTESIAN = 0x01,
  TRAJECTORY_RADIAL    = 0x02,
  TRAJECTORY_SPIRAL    = 0x04,
  TRAJECTORY_BLADE     = 0x08
}; =#

struct mdhLC 
  ushLine::UInt16
  ushAcquisition::UInt16
  ushSlice::UInt16
  ushPartition::UInt16
  ushEcho::UInt16
  ushPhase::UInt16
  ushRepetition::UInt16
  ushSet::UInt16
  ushSeg::UInt16
  ushIda::UInt16
  ushIdb::UInt16
  ushIdc::UInt16
  ushIdd::UInt16
  ushIde::UInt16
end

struct mdhCutOff
  ushPre::UInt16
  ushPost::UInt16
end

struct mdhSlicePosVec
  flSag::Float32
  flCor::Float32
  flTra::Float32
end


struct PMUdata 
	data::UInt16
	trigger::UInt16
end

struct mdhSliceData
  sSlicePosVec::mdhSlicePosVec
  aflQuaternion::NTuple{4,Float32}
end

# This is the VB line header 
struct sMDH
  ulFlagsAndDMALength::UInt32
  lMeasUID::UInt32
  ulScanCounter::UInt32
  ulTimeStamp::UInt32
  ulPMUTimeStamp::UInt32
  aulEvalInfoMask::NTuple{2,UInt32}
  ushSamplesInScan::UInt16
  ushUsedChannels::UInt16
  sLC::mdhLC
  sCutOff::mdhCutOff

  ushKSpaceCentreColumn::UInt16
  ushCoilSelect::UInt16
  fReadOutOffcentre::Float32
  ulTimeSinceLastRF::UInt32
  ushKSpaceCentreLineNo::UInt16
  ushKSpaceCentrePartitionNo::UInt16
  aushIceProgramPara::NTuple{4,UInt16}
  aushFreePara::NTuple{4,UInt16}

  sSliceData::mdhSliceData

  ushChannelId::UInt16
  ushPTABPosNeg::UInt16
end

# This is the VD line header 
struct sScanHeader
	  ulFlagsAndDMALength::UInt32
	  lMeasUID::Int32
	  ulScanCounter::UInt32
	  ulTimeStamp::UInt32
	  ulPMUTimeStamp::UInt32
	  ushSystemType::UInt16
	  ulPTABPosDelay::UInt16
	  lPTABPosX::Int32
	  lPTABPosY::Int32
	  lPTABPosZ::Int32
	  ulReserved1::UInt32
	  aulEvalInfoMask::NTuple{2,UInt32}
	  ushSamplesInScan::UInt16
	  ushUsedChannels::UInt16
	  sLC::mdhLC
	  sCutOff::mdhCutOff
	  ushKSpaceCentreColumn::UInt16
	  ushCoilSelect::UInt16
	  fReadOutOffcentre::Float32
	  ulTimeSinceLastRF::UInt32
	  ushKSpaceCentreLineNo::UInt16
	  ushKSpaceCentrePartitionNo::UInt16
	  sSliceData::mdhSliceData
	  aushIceProgramPara::NTuple{24,UInt16}
	  aushReservedPara::NTuple{4,UInt16}
	  ushApplicationCounter::UInt16
	  ushApplicationMask::UInt16
	  ulCRC::UInt32
end

struct sChannelHeader
  ulTypeAndChannelLength::UInt32
  lMeasUID::Int32
  ulScanCounter::UInt32
  ulReserved1::UInt32
  ulSequenceTime::UInt32
  ulUnused2::UInt32
  ulChannelId::UInt16
  ulUnused3::UInt16
  ulCRC::UInt32
end

struct MrParcRaidFileEntry
  measId_::UInt32
  fileId_::UInt32
  off_::UInt64
  len_::UInt64
  patName_::NTuple{64,UInt8}
  protName_::NTuple{64,UInt8}
end

struct MrParcRaidFileHeader
  hdSize_::UInt32
  count_::UInt32
end

#=
typedef struct
{
  sMDH mdh;
  void* previous;
  void* next;
  float* data;
} SiemensMdhNode;=#

struct SiemensBaseParameters
  matrix_size::NTuple{3,UInt32}
  pat_matrix_size::NTuple{3,UInt32}
  base_resolution::UInt32
  phase_encoding_lines::UInt32
  partitions::UInt32
  dimensions::UInt32
  phase_resolution::Float32
  slice_resolution::Float32
  dwell_time_us::Float32
  acceleration_factor_pe::UInt32
  acceleration_factor_3d::UInt32
end