# # [Load and Save Raw Data](@id load-save)
#
# MRIFiles reads raw data from [ISMRMRD](https://ismrmrd.readthedocs.io/) files and from Bruker datasets, and it
# writes ISMRMRD files. This guide shows how to load and save ISMRMRD files, how to save acquisition data, and how to
# load a Bruker dataset and convert it into an ISMRMRD file. [Data Model](@ref data-model) explains the two data types
# that appear here, `RawAcquisitionData` and `AcquisitionData`.
#
# ## Loading an ISMRMRD File
#
# An `ISMRMRDFile` refers to a file by its name, and `AcquisitionData` reads it into the form that the reconstruction
# needs. The example loads the knee data from [Direct Reconstruction](@ref direct-reconstruction):

using MRIReco, MRIFiles

using Pkg.Artifacts
artifacts = joinpath(pkgdir(MRIReco), "docs", "Artifacts.toml")
datadir = ensure_artifact_installed("MRIRecoDocsData", artifacts; io = devnull)
file = ISMRMRDFile(joinpath(datadir, "knee_3dFSE_slice170.h5"))
acqData = AcquisitionData(file)
nothing #hide

# `AcquisitionData(file)` reads the file into a `RawAcquisitionData` object, which holds the header and the profiles
# as they are stored in the file, and then converts this object. To look at the raw data or to change it before the
# conversion, do the two steps separately:

rawData = RawAcquisitionData(file)
acqData = AcquisitionData(rawData)
length(rawData.profiles), size(acqData.kdata[1])

# The file contains 320 profiles, and the converted data has 102400 samples from 8 coils.
# [Bin Profiles into Contrasts](@ref binning) is an example of changing raw data before the conversion.
#
# Files from a scanner can also contain readouts that are not part of the image, such as noise measurements. The
# conversion does not skip them, so remove them first with [`remove_raw_by_flags`](@ref), for example
# `rawData = remove_raw_by_flags(rawData)`. [Data Model](@ref data-model) describes what the conversion does and what
# it expects.
#
# An ISMRMRD file stores the data of an acquisition in an HDF5 group, and MRIFiles reads the group `dataset` by
# default. If the group of a file has another name, loading fails with a `KeyError`, and the name has to be passed as
# the second argument, for example `AcquisitionData(file, "scan1")`.
#
# Raw data from Siemens, GE and Philips scanners can be converted into ISMRMRD files with the converters of the
# ISMRMRD project: [`siemens_to_ismrmrd`](https://github.com/ismrmrd/siemens_to_ismrmrd),
# [`ge_to_ismrmrd`](https://github.com/ismrmrd/ge_to_ismrmrd) and
# [`philips_to_ismrmrd`](https://github.com/ismrmrd/philips_to_ismrmrd).
#
# ## Saving Raw Data
#
# `save` writes a `RawAcquisitionData` object to an ISMRMRD file:

outfile = joinpath(mktempdir(), "knee.h5")
MRIFiles.save(ISMRMRDFile(outfile), rawData)

# The file then contains the header parameters and the profiles of `rawData`. An existing file is overwritten, and
# the name of the group can be passed as a third argument. MRIFiles does not export `save`
# ([#315](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/315)), so it is called with the module name.
#
# ## Saving Acquisition Data
#
# `save` accepts only raw data. Convert `AcquisitionData`, for example from a simulation, with `RawAcquisitionData`
# first. For non-Cartesian data that was loaded from a file, this conversion can fail
# ([#346](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/346)). Here, the knee data is undersampled
# with MRISampling and saved:

using MRISampling

acqDataSub = sample_kspace(acqData, 4.0, "poisson", calsize = 30, seed = 1234)
subfile = joinpath(mktempdir(), "knee_undersampled.h5")
MRIFiles.save(ISMRMRDFile(subfile), RawAcquisitionData(acqDataSub))

# For this Cartesian data, loading the file returns the same k-space data and sampling pattern:

acqDataLoaded = AcquisitionData(ISMRMRDFile(subfile))
acqDataLoaded.kdata == acqDataSub.kdata,
    acqDataLoaded.subsampleIndices == acqDataSub.subsampleIndices

# A file written from acquisition data does not keep everything:
#
# - The samples and the k-space positions are stored in single precision, so double-precision data, for example
#   from [`simulation`](@ref), is loaded in single precision. The field of view is stored with the precision of the
#   data.
# - Only `"Cartesian"` and `"Cartesian3D"` trajectories are loaded as Cartesian trajectories again. Any other
#   trajectory is loaded as a non-Cartesian trajectory with the name `"Custom"` and the same nodes. This includes
#   `"EPI"` trajectories, which MRIReco otherwise treats as Cartesian, so their reconstruction changes. Data whose
#   trajectory already had the name `"Custom"` is loaded with additional profiles of zeros if it has several slices
#   or repetitions ([#331](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/331)).
# - The loaded `"Custom"` trajectory is always marked as circular. If the original trajectory was not,
#   reconstructions such as the direct one now set the corners of the image to zero. The `subsampleIndices` of the
#   loaded data can differ as well, which is one reason why converting such data into raw data again can fail.
# - After loading, the echo time is zero, `sequenceInfo` is empty and the readout duration returned by
#   `acqTimePerProfile` is the default of 1 ms. The sampling times of a non-Cartesian trajectory are rebuilt from the
#   original readout duration, those of a Cartesian trajectory from the default. This matters if the reconstruction
#   uses these times, for example for an off-resonance correction.
# - The header has only five entries, among them the matrix size and the field of view. Elements that the ISMRMRD
#   schema requires, such as `reconSpace` and `encodingLimits`, are missing, so other ISMRMRD software may not accept
#   the file ([#33](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/33)).
#
# ## Loading a Bruker Dataset
#
# A Bruker dataset is the folder of one scan, which ParaVision names by its number, inside the folder of the study.
# MRIFiles reads several files of this folder, among them `acqp`, `method`, the raw data (`fid` for ParaVision 6,
# `rawdata.job0` for ParaVision 360) and the file `visu_pars` of a reconstruction in the subfolder `pdata`, so copy
# the folder as a whole. A `BrukerFile` is created from the path of this folder. The data is loaded in two steps,
# because `AcquisitionData` does not accept a `BrukerFile`. The example loads a Cartesian dataset that is part of the
# documentation data:

file = BrukerFile(joinpath(datadir, "brukerfileCart"))
rawData = RawAcquisitionData(file)
acqData = AcquisitionData(rawData)
encodingSize(acqData), numSlices(acqData)

# The data can then be reconstructed as in [Direct Reconstruction](@ref direct-reconstruction). Leave out
# `params[:reconSize]`, which the tutorial sets to the size of the knee data. The image then has the matrix size of
# the data, here 96 × 96 for each of the 15 slices.
#
# MRIFiles reads Cartesian acquisitions from ParaVision 6 and ParaVision 360. The reader is tested with FLASH and
# RARE data from ParaVision 6 and with FLASH data from ParaVision 360. MRIFiles also reads 3D UTE acquisitions from
# ParaVision 6. It recognizes them by the protocol name `UTE3D`. For a UTE scan with another protocol name, call
# `RawAcquisitionData_3DUTE(file)` instead of `RawAcquisitionData(file)`. For a UTE scan with several repetitions,
# the conversion into `AcquisitionData` currently adds profiles of zeros to each repetition and the same number of
# profiles at the center of k-space to the trajectory
# ([#331](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/331)).
#
# If the field of view was shifted along a phase encoding direction when the scan was planned, the image
# reconstructed with MRIReco is shifted relative to the image from the scanner. For Cartesian data,
# `AcquisitionData(rawData, OffsetBruker = true)` corrects this shift.
#
# `recoData(file)` returns the images that the scanner reconstructed (`pdata/1/2dseq`), and a second argument selects
# another reconstruction by its number in `pdata`. The matrix size and orientation of these images can differ from
# those of an image reconstructed with MRIReco. Parameters of the `method` file can be read by indexing the
# `BrukerFile` with their name. `file["PVM_Fov"]`, for example, returns the field of view as a vector of strings, and
# indexing with a name that does not exist returns an empty string.
#
# ## Converting a Bruker Dataset into an ISMRMRD File
#
# Saving the raw data of a Bruker dataset with `save` converts the dataset into an ISMRMRD file, which can be read
# like any other:

brukerOutfile = joinpath(mktempdir(), "brukerfileCart.h5")
MRIFiles.save(ISMRMRDFile(brukerOutfile), rawData)
length(RawAcquisitionData(ISMRMRDFile(brukerOutfile)).profiles)


# As with a file written from acquisition data, the header of this file lacks elements that the ISMRMRD schema
# requires ([#33](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/33)). MRIFiles does not write Bruker
# datasets.
