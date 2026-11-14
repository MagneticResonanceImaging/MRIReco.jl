export softwareCoilCompression, geometricCoilCompression, applyCoilCompression, applyCoilCompressionSensitivityMaps

#
# Coil compression replaces the channels by a smaller number of virtual coils, which are linear
# combinations of the channels. The combination is given by a compression matrix `ccMat` of size
# channels × virtual coils, in one of three forms:
#
# - a `Matrix`: one matrix for all the data (software coil compression of an array),
# - a `Vector` of matrices: one matrix per slice (software coil compression of an `AcquisitionData`),
# - a 3D `Array` (positions × channels × virtual coils): one matrix per position along the fully
#   sampled dimension `dim` (geometric coil compression).
#
# Sensitivity maps (x, y, z, channels) and the k-space returned by `kDataCart`
# (x, y, z, channels, contrasts, repetitions) both store the channels along dimension 4,
# so the same functions apply the compression to both.
#

const PerSliceCCMat = AbstractVector{<:AbstractMatrix}
const GeometricCCMat = AbstractArray{<:Number, 3}

#
# Software coil compression (SCC)
#

"""
    softwareCoilCompression(kdata::Matrix{T}, numVC::Int = size(kdata, 2)) where T <: Complex

Compress the k-space matrix `kdata` (samples × channels) to `numVC` virtual coils.

The compression matrix `ccMat` (channels × `numVC`) consists of the right singular vectors of `kdata`
that belong to the `numVC` largest singular values. Return `(kdataCC, ccMat)` with
`kdataCC = kdata * ccMat` (samples × `numVC`).

# Reference
- Huang F, Vijayakumar S, Li Y, Hertel S, Duensing GR. A software channel compression technique
  for faster reconstruction with many channels. Magnetic Resonance Imaging 2008;26:133–141.
  doi: 10.1016/j.mri.2007.04.010
"""
function softwareCoilCompression(kdata::Matrix{T}, numVC::Int = size(kdata, 2)) where T <: Complex
    checkNumVC(numVC, size(kdata, 2))
    ccMat = svd(kdata).V[:, 1:numVC]
    return kdata * ccMat, ccMat
end

"""
    softwareCoilCompression(kspace::Array{T,6}, numVC::Int = size(kspace, 4); sContr = 1, sRep = 1) where T <: Complex

Compress the Cartesian `kspace` (x, y, z, channels, contrasts, repetitions), as returned by
[`kDataCart`](@ref), to `numVC` virtual coils.

A single compression matrix is computed from contrast `sContr` and repetition `sRep` and applied to
all contrasts and repetitions with [`applyCoilCompression`](@ref).
Return `(kspaceCC, ccMat)`, with `ccMat` of size channels × `numVC`.

# Reference
- Huang F, Vijayakumar S, Li Y, Hertel S, Duensing GR. A software channel compression technique
  for faster reconstruction with many channels. Magnetic Resonance Imaging 2008;26:133–141.
  doi: 10.1016/j.mri.2007.04.010
"""
function softwareCoilCompression(kspace::Array{T,6}, numVC::Int = size(kspace, 4); sContr::Int = 1, sRep::Int = 1) where T <: Complex
    nChan = size(kspace, 4)
    _, ccMat = softwareCoilCompression(reshape(kspace[:, :, :, :, sContr, sRep], :, nChan), numVC)
    return applyCoilCompression(kspace, ccMat), ccMat
end

function softwareCoilCompression(kspace::Array{T,6}, numVC::Int, sContr::Int, sRep::Int = 1) where T <: Complex
    Base.depwarn("`softwareCoilCompression(kspace, numVC, sContr, sRep)` is deprecated, use " *
                 "`softwareCoilCompression(kspace, numVC; sContr = sContr, sRep = sRep)` instead.",
                 :softwareCoilCompression)
    return softwareCoilCompression(kspace, numVC; sContr, sRep)
end

"""
    softwareCoilCompression(acqData::AcquisitionData, numVC::Int = numChannels(acqData); sContr = 1, sRep = 1)

Compress `acqData` to `numVC` virtual coils. Works for any trajectory and for undersampled data.

One compression matrix is computed per slice, from contrast `sContr` and repetition `sRep`, and
applied to all contrasts and repetitions of that slice with [`applyCoilCompression`](@ref).
Return `(acqDataCC, ccMat)`, where `ccMat` is a vector with one matrix (channels × `numVC`) per slice.
It can be passed to [`applyCoilCompressionSensitivityMaps`](@ref) to compress sensitivity maps.

# Reference
- Huang F, Vijayakumar S, Li Y, Hertel S, Duensing GR. A software channel compression technique
  for faster reconstruction with many channels. Magnetic Resonance Imaging 2008;26:133–141.
  doi: 10.1016/j.mri.2007.04.010
"""
function softwareCoilCompression(acqData::AcquisitionData{T,D}, numVC::Int = numChannels(acqData); sContr::Int = 1, sRep::Int = 1) where {T,D}
    ccMat = Matrix{Complex{T}}[softwareCoilCompression(acqData.kdata[sContr, sl, sRep], numVC)[2]
                               for sl in 1:numSlices(acqData)]
    return applyCoilCompression(acqData, ccMat), ccMat
end

#
# Geometric coil compression (GCC)
#

"""
    geometricCoilCompression(kspace::Array{T,6}, numVC::Int = size(kspace, 4); dim = 1, sContr = 1, sRep = 1) where T <: Complex

Compress the Cartesian `kspace` (x, y, z, channels, contrasts, repetitions), as returned by
[`kDataCart`](@ref), to `numVC` virtual coils with geometric coil compression.

Dimension `dim` (1, 2 or 3) must be fully sampled, typically the readout direction. The k-space is
inverse Fourier transformed along `dim`, and one compression matrix is computed for each position
along `dim` from contrast `sContr` and repetition `sRep`. The matrices are aligned so that the
virtual coils vary smoothly between positions, and applied to all contrasts and repetitions with
[`applyCoilCompression`](@ref).

Return `(kspaceCC, ccMat)`, where `ccMat` has size positions × channels × `numVC`.

# Reference
- Zhang T, Pauly JM, Vasanawala SS, Lustig M. Coil compression for accelerated imaging with Cartesian
  sampling. Magnetic Resonance in Medicine 2013;69:571–582. doi: 10.1002/mrm.24267
"""
function geometricCoilCompression(kspace::Array{T,6}, numVC::Int = size(kspace, 4); dim::Int = 1, sContr::Int = 1, sRep::Int = 1) where T <: Complex
    ccMat = gccMatrices(kspace, numVC, dim, sContr, sRep)
    return applyCoilCompression(kspace, ccMat; dim), ccMat
end

"""
    geometricCoilCompression(acqData::AcquisitionData, numVC::Int = numChannels(acqData); dim = 1, sContr = 1, sRep = 1)

Compress the Cartesian `acqData` to `numVC` virtual coils with geometric coil compression.
See the method for arrays for the meaning of `dim`, `sContr` and `sRep`.

The compression matrices are computed from the k-space returned by [`kDataCart`](@ref) and shared by
all slices. For 2D data, `dim` must be 1 or 2. The result keeps the sampling pattern, the
trajectories and the other fields of `acqData`, so undersampled data stay undersampled.
Return `(acqDataCC, ccMat)`, where `ccMat` has size positions × channels × `numVC`.
"""
function geometricCoilCompression(acqData::AcquisitionData{T,D}, numVC::Int = numChannels(acqData); dim::Int = 1, sContr::Int = 1, sRep::Int = 1) where {T,D}
    checkGCCDim(acqData, dim)
    ccMat = gccMatrices(kDataCart(acqData), numVC, dim, sContr, sRep)
    return applyCoilCompression(acqData, ccMat; dim), ccMat
end

# Compute the aligned GCC matrices (positions along `dim` × channels × numVC)
function gccMatrices(kspace::Array{T,6}, numVC::Int, dim::Int, sContr::Int, sRep::Int) where T <: Complex
    dim in 1:3 || throw(ArgumentError("the compression dimension `dim` must be 1, 2 or 3, got $dim"))
    nChan = size(kspace, 4)
    checkNumVC(numVC, nChan)

    calib = hybridSpace(kspace[:, :, :, :, sContr, sRep], dim)
    ccMat = zeros(T, size(kspace, dim), nChan, numVC)
    for p in axes(ccMat, 1)
        A = reshape(selectdim(calib, dim, p), :, nChan) # (other positions) × channels
        ccMat[p, :, :] = svd(A).V[:, 1:numVC]
    end
    return alignCCMtx(ccMat, numVC)
end

"""
    alignCCMtx(ccMat::AbstractArray{T,3}, numVC::Int = size(ccMat, 3))

Align the geometric coil compression matrices `ccMat` (positions × channels × virtual coils).

The SVD determines the virtual coils of each position only up to a unitary transformation.
Starting from the central position, the matrix of each position is rotated to be as close as
possible to the matrix of its neighbor (orthogonal Procrustes problem), so that the virtual
coils vary smoothly. `ccMat` is modified in place and returned.

# Reference
- Zhang T, Pauly JM, Vasanawala SS, Lustig M. Coil compression for accelerated imaging with Cartesian
  sampling. Magnetic Resonance in Medicine 2013;69:571–582. doi: 10.1002/mrm.24267
"""
function alignCCMtx(ccMat::AbstractArray{T,3}, numVC::Int = size(ccMat, 3)) where T
    nPos = size(ccMat, 1)
    n0 = max(1, nPos ÷ 2) # align everything based on the central position

    for range in (n0-1:-1:1, n0+1:nPos) # backwards to the first, then forward to the last position
        A0 = ccMat[n0, :, 1:numVC]
        for n in range
            A1 = ccMat[n, :, 1:numVC]
            usv = svd(A1' * A0)
            ccMat[n, :, 1:numVC] = A1 * usv.U * usv.V'
            A0 = ccMat[n, :, 1:numVC]
        end
    end
    return ccMat
end

#
# Applying compression matrices
#

"""
    applyCoilCompression(kspace::Array{T,6}, ccMat; dim = 1) where T <: Complex

Apply the compression matrix `ccMat` to the Cartesian `kspace` (x, y, z, channels, contrasts,
repetitions), as returned by [`kDataCart`](@ref), for all contrasts and repetitions.

`ccMat` can be
- a matrix (channels × virtual coils), as returned by [`softwareCoilCompression`](@ref) for an array,
- a vector with one matrix per slice (dimension 3), as returned by
  [`softwareCoilCompression`](@ref) for an `AcquisitionData`,
- a 3D array (positions × channels × virtual coils), as returned by
  [`geometricCoilCompression`](@ref). The k-space is then inverse Fourier transformed along `dim`,
  compressed position by position and transformed back.

Return the compressed k-space.
"""
applyCoilCompression(kspace::Array{T,6}, ccMat::Union{AbstractMatrix, PerSliceCCMat}) where T <: Complex =
    compressChannels(kspace, ccMat)

function applyCoilCompression(kspace::Array{T,6}, ccMat::GeometricCCMat; dim::Int = 1) where T <: Complex
    hybrid = compressChannels(hybridSpace(kspace, dim), ccMat, dim)
    return fftshift(fft(ifftshift(hybrid, dim), dim), dim)
end

"""
    applyCoilCompression(acqData::AcquisitionData, ccMat; dim = 1)

Apply the compression matrix `ccMat` to all contrasts, slices and repetitions of `acqData`.

`ccMat` can be
- a matrix (channels × virtual coils), used for every slice,
- a vector with one matrix per slice, as returned by [`softwareCoilCompression`](@ref),
- a 3D array (positions × channels × virtual coils), as returned by
  [`geometricCoilCompression`](@ref). This requires Cartesian data, with `dim` fully sampled.

Return a new `AcquisitionData`. The sampling pattern, the trajectories and the other fields are
copied from `acqData`.
"""
function applyCoilCompression(acqData::AcquisitionData{T,D}, ccMat::PerSliceCCMat) where {T,D}
    nContr, nSl, nRep = size(acqData.kdata)
    length(ccMat) == nSl || throw(DimensionMismatch("acqData has $nSl slices but $(length(ccMat)) compression matrices were given"))

    kdataCC = [Matrix{Complex{T}}(acqData.kdata[contr, sl, rep] * ccMat[sl])
               for contr in 1:nContr, sl in 1:nSl, rep in 1:nRep]
    return withKData(acqData, kdataCC)
end

applyCoilCompression(acqData::AcquisitionData, ccMat::AbstractMatrix) =
    applyCoilCompression(acqData, fill(ccMat, numSlices(acqData)))

function applyCoilCompression(acqData::AcquisitionData{T,D}, ccMat::GeometricCCMat; dim::Int = 1) where {T,D}
    checkGCCDim(acqData, dim)
    kspaceCC = applyCoilCompression(kDataCart(acqData), ccMat; dim)

    # Write the compressed k-space back at the sampled positions (inverse of kDataCart).
    # Missing profiles are zero for all positions along `dim` and stay zero after the compression.
    nContr, nSl, nRep = size(acqData.kdata)
    numVC = size(kspaceCC, 4)
    kspaceCC = reshape(kspaceCC, :, nSl, numVC, nContr, nRep)
    kdataCC = [Matrix{Complex{T}}(kspaceCC[acqData.subsampleIndices[contr], sl, :, contr, rep])
               for contr in 1:nContr, sl in 1:nSl, rep in 1:nRep]
    return withKData(acqData, kdataCC)
end

"""
    applyCoilCompressionSensitivityMaps(smaps::AbstractArray{T,4}, ccMat; dim = 1) where T <: Complex

Compress the sensitivity maps `smaps` (x, y, z, channels) with the compression matrix `ccMat` used
for the data, so that they are the sensitivities of the virtual coils:
``S_v(\\mathbf{r}) = \\sum_c S_c(\\mathbf{r})\\,A_{cv}``.

`ccMat` can be a matrix, a vector with one matrix per slice (dimension 3 of `smaps`) or a 3D array
from [`geometricCoilCompression`](@ref), see [`applyCoilCompression`](@ref). For geometric coil
compression, the matrices are applied position by position along `dim`; no Fourier transform is
needed because the maps are already in image space.
"""
applyCoilCompressionSensitivityMaps(smaps::AbstractArray{T,4}, ccMat::Union{AbstractMatrix, PerSliceCCMat}) where T <: Complex =
    compressChannels(smaps, ccMat)

applyCoilCompressionSensitivityMaps(smaps::AbstractArray{T,4}, ccMat::GeometricCCMat; dim::Int = 1) where T <: Complex =
    compressChannels(smaps, ccMat, dim)

#
# Internal helpers
#

# Multiply `data` along its 4th (channel) dimension by `ccMat` (channels × virtual coils).
function compressChannels(data::AbstractArray{<:Number}, ccMat::AbstractMatrix)
    ndims(data) >= 4 || throw(ArgumentError("the channels must be stored along dimension 4, got an array with $(ndims(data)) dimensions"))
    sz = size(data)
    nChan, numVC = size(ccMat)
    sz[4] == nChan || throw(DimensionMismatch("data has $(sz[4]) channels but ccMat has $nChan rows"))

    # (dims 1:3) × channels × (dims 5:end)
    n1, n2 = prod(sz[1:3]), prod(sz[5:end])
    A = reshape(data, n1, nChan, n2)
    out = Array{promote_type(eltype(data), eltype(ccMat))}(undef, n1, numVC, n2)
    for k in 1:n2
        mul!(view(out, :, :, k), view(A, :, :, k), ccMat)
    end
    return reshape(out, sz[1:3]..., numVC, sz[5:end]...)
end

# One matrix per block along dimension `dim`. Dimension `dim` is split into `length(ccMats)` blocks
# of equal size: one block per slice for 2D multi-slice data, a single block for a 3D volume,
# or one block per slab for 3D multi-slab data.
function compressChannels(data::AbstractArray{<:Number}, ccMats, dim::Int)
    nBlocks = size(ccMats, 1)
    blockMatrix(b) = ccMats isa AbstractVector ? ccMats[b] : ccMats[b, :, :]
    nBlocks == 1 && return compressChannels(data, blockMatrix(1))

    n = size(data, dim)
    n % nBlocks == 0 || throw(DimensionMismatch("size of dimension $dim ($n) is not a multiple of the number of compression matrices ($nBlocks)"))
    blk = n ÷ nBlocks
    sz = size(data)
    numVC = size(blockMatrix(1), 2)
    out = Array{promote_type(eltype(data), eltype(blockMatrix(1)))}(undef, sz[1:3]..., numVC, sz[5:end]...)
    for b in 1:nBlocks
        idx = (b - 1) * blk + 1 : b * blk
        # copy the block so that compressChannels works on a contiguous array (BLAS)
        selectdim(out, dim, idx) .= compressChannels(collect(selectdim(data, dim, idx)), blockMatrix(b))
    end
    return out
end

compressChannels(data::AbstractArray{<:Number}, ccMat::PerSliceCCMat) = compressChannels(data, ccMat, 3)

# Inverse Fourier transform along `dim` only (centered)
hybridSpace(kspace::AbstractArray, dim::Int) = fftshift(ifft(ifftshift(kspace, dim), dim), dim)

# New AcquisitionData with the k-space data `kdata` and the other fields copied from `acqData`
withKData(acqData::AcquisitionData{T,D}, kdata) where {T,D} =
    AcquisitionData{T,D}(deepcopy(acqData.sequenceInfo), deepcopy(acqData.traj), kdata,
                         deepcopy(acqData.subsampleIndices), acqData.encodingSize, acqData.fov)

checkNumVC(numVC, nChan) = 1 <= numVC <= nChan ||
    throw(ArgumentError("the number of virtual coils must be between 1 and the number of channels ($nChan), got $numVC"))

function checkGCCDim(acqData::AcquisitionData{T,D}, dim::Int) where {T,D}
    dim in 1:3 || throw(ArgumentError("the compression dimension `dim` must be 1, 2 or 3, got $dim"))
    D == 2 && dim == 3 && throw(ArgumentError("for 2D data, `dim` must be 1 or 2 (dimension 3 contains the slices)"))
    all(isCartesian, acqData.traj) || throw(ArgumentError("geometric coil compression requires Cartesian data"))

    # every sampled line along `dim` must be complete
    N = acqData.encodingSize
    for idx in acqData.subsampleIndices
        mask = falses(prod(N))
        mask[idx] .= true
        mask = reshape(mask, N)
        all(mask .== any(mask, dims = dim)) || throw(ArgumentError(
            "dimension $dim is not fully sampled; geometric coil compression requires a fully sampled dimension, usually the readout"))
    end
    return nothing
end
