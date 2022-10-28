export CoilCompression, softwareCoilCompression, applyCoilCompressionSensitivityMaps
function CoilCompression(data::Union{AcquisitionData{T,D},Matrix{T}},numVC::Int64 = nothing; CCMode = "SCC") where {T<:Complex,D}
    if isnothing(numVC)
        if typeof(data) <: AcquisitionData
            numVC = size(data.kdata,2)
        else
            numVC = size(kdata,2)
        end
    end

    if (CCMode == "SCC")
        return softwareCoilCompression(data,numVC)
    elseif (CCMode == "GCC")
        return geometricCoilCompression(data,numVC) #perform gcc along readout
    elseif (CCMode == "ECC")
        @error "ECC currently not implemented"
    else
        @error "Supported mode are : SCC or GCC (ECC not implemented)"
    end
end

"""
    geometricCoilCompression(kdata::Matrix{T}, numVC::Int64 = size(kdata, 2)) where T <: Complex

return a SVD-based Geometric Coil Compression (GCC) matrix for `numVC` virtual coils.
Coil compression is performed for each position along the readout (kx) direction.

Input :
    -   kspace::Array{T,4} : dimension : kx,ky,kz,channel

Reference :
- Zhang T, Pauly JM, Vasanawala SS, Lustig M. Coil compression for accelerated imaging with Cartesian sampling. Magnetic Resonance in Medicine 2013;69:571–582 doi: 10.1002/mrm.24267.
"""
function geometricCoilCompression(kspace::Array{T,6}, numVC::Int64 = size(kdata, 2)) where T <: Complex
   @error "not implemented"
end

"""
    softwareCoilCompression(kdata::Matrix{T}, numVC::Int64 = size(kdata, 2);) where T <: Complex

Perform SVD-based coil compression for the given `kdata` and `smaps`

Reference :
- Huang F, Vijayakumar S, Li Y, Hertel S, Duensing GR. A software channel compression technique
for faster reconstruction with many channels. Magnetic Resonance Imaging 2008;26:133–141 doi: 10.1016/j.mri.2007.04.010.
"""
function softwareCoilCompression(kdata::Matrix{T}, numVC::Int64 = size(kdata, 2);) where T <: Complex
  usv = svd(kdata)
  ccMat = usv.Vt[:, 1:numVC]
  kdataCC = kdata * ccMat

  return kdataCC, ccMat
end

"""
softwareCoilCompression(acqData::AcquisitionData{T,D}, numVC::Int64 = size(acqData.kdata[1,1,1],2);contr_::Int = 1,rep_::Int = 1) where {T,D}

Perform SVD-based coil compression for the given `kdata` and `smaps`.
By default the first repetition/echoes is used to estimates the coil compression matrix

Reference :
- Huang F, Vijayakumar S, Li Y, Hertel S, Duensing GR. A software channel compression technique
for faster reconstruction with many channels. Magnetic Resonance Imaging 2008;26:133–141 doi: 10.1016/j.mri.2007.04.010.
"""
function softwareCoilCompression(acqData::AcquisitionData{T,D}, numVC::Int64 = size(acqData.kdata[1,1,1],2);sContr::Int = 1,sRep::Int = 1) where {T,D}
    acqDataCC = deepcopy(acqData)
    nSl = numSlices(acqData)
    ccMat_vec = Matrix{Complex{T}}[]

    for sl = 1:nSl
      acqDataCC.kdata[sContr, sl, sRep],ccMat = softwareCoilCompression(acqData.kdata[sContr, sl, sRep],numVC)
      push!(ccMat_vec,ccMat)
      # compress kdata-slice
      for rep = 1:numRepetitions(acqData), contr = 1:numContrasts(acqData)
        if(rep == sRep) && (contr == sContr) # already performed above
            continue
        end
        acqDataCC.kdata[contr, sl, rep] = acqData.kdata[contr, sl, rep] * ccMat
      end
    end

    return acqDataCC, ccMat_vec
  end

"""
    applyCoilCompressionSensitivityMaps(smaps::Array{T,4}, ccMat::Array{T}) where T <: Complex

Apply coil compression for the given sensitivty maps `smaps` using the Coil Compression Matrix ccMat
"""
function applyCoilCompressionSensitivityMaps(smaps::Array{T,4}, ccMat::Array{T}) where T <: Complex
  nx, ny, nz, nc = size(smaps)

  smapsCC = zeros(T, nx, ny, nz, numVC)
  for j = 1:nz, i = 1:ny
    smapsCC[:, i, j, :] .= smaps[:, i, j, :] * ccMat
  end

  return smapsCC
end

"""
    applyCoilCompressionSensitivityMaps(smaps::Array{T,4}, ccMat_vec::Vector{Matrix{T}}) where T <: Complex

Apply coil compression to multislice 2D sensitivity maps `smaps` using a vector of Coil Compression Matrix ccMat_vec. Each indices corrsponds to the number of slice
"""
function applyCoilCompressionSensitivityMaps(smaps::Array{T,4}, ccMat_vec::Vector{Matrix{T}}) where T <: Complex
    nx, ny, nz, nc = size(smaps)

    if nz != length(ccMat_vec)
        @error "Number of Coil compression matrix is not equal to the number of slice"
    end

    smapsCC = zeros(T, nx, ny, nz, numVC)
    for j = 1:nz, i = 1:ny
      smapsCC[:, i, j, :] .= smaps[:, i, j, :] * ccMat_vec[j]
    end

    return smapsCC
  end
