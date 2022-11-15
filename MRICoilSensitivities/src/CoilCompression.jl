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
    geometricCoilCompression(kspace::Matrix{T}, numVC::Int64 = size(kdata, 2); ws::Int = 1) where T <: Complex

return a SVD-based Geometric Coil Compression (GCC) matrix for `numVC` virtual coils.
Coil compression is performed for each position along the readout (kx) direction.

Input :
    -   kspace::Array{T,6} : dimension : kx,ky,kz,channel,contrasts,repetitions
    -   numVS::Int

Optional :
    -   ws::Int = 1 : sliding window

Reference :
- Zhang T, Pauly JM, Vasanawala SS, Lustig M. Coil compression for accelerated imaging with Cartesian sampling. Magnetic Resonance in Medicine 2013;69:571–582 doi: 10.1002/mrm.24267.
"""
function geometricCoilCompression(kspace::Array{T,6}, numVC::Int = size(kdata, 4); ws::Int = 1) where T <: Complex
    @error "not implemented yet"
    #=
    # round sliding window size to nearest odd number
    ws = floor(ws/2)*2 + 1;

    # perform ifft along 1st dimension
    im_data = ifftshift(ifft(ifftshift(kspace),1))
    for kx in axes(kspace,1)

    end

    #=
    #Calculate compression matrices for each position in the readout over a sliding window of size ws
    kspaceCC = zeros(T,Nc,min(Nc,ws*Ny*Nz),Nx);
    zpim = zpad(im,[Nx + ws-1,Ny,Nz,Nc]);
    for n = [1:Nx]
        %tmpc = squeeze(im(n,:,:));
        tmpc = reshape(zpim(n:n+ws-1,:,:,:),ws*Ny*Nz,Nc);
        [U,S,V] = svd(tmpc,'econ');
        mtx(:,:,n) = V;

        usv = svd(kdata)
        ccMat = usv.Vt[:, 1:numVC]
        kdataCC = kdata * ccMat
    end
    =#
    =#
end

"""
    softwareCoilCompression(kdata::Matrix{T}, numVC::Int64 = size(kdata, 2);) where T <: Complex

Perform SVD-based coil compression for the given `kdata` and `smaps`

Reference :
- Huang F, Vijayakumar S, Li Y, Hertel S, Duensing GR. A software channel compression technique
for faster reconstruction with many channels. Magnetic Resonance Imaging 2008;26:133–141 doi: 10.1016/j.mri.2007.04.010.
"""
function softwareCoilCompression(kdata::Matrix{T}, numVC::Int64 = size(kdata, 2)) where T <: Complex
  usv = svd(kdata)
  ccMat = usv.Vt[:, 1:numVC]
  kdataCC = kdata * ccMat

  return kdataCC, ccMat
end

"""
    softwareCoilCompression(kspace::Array{T,6}, numVC::Int64 = size(kspace, 4),sContr::Int = 1,sRep::Int = 1) where T <: Complex

Perform SVD-based coil compression for the given `kspace` [sx,sy,sz,channels, contrasts,repetitions].
Default extraction from kDataCart.

## To do :
-   implement multislice reconstruction

## Reference :
- Huang F, Vijayakumar S, Li Y, Hertel S, Duensing GR. A software channel compression technique
for faster reconstruction with many channels. Magnetic Resonance Imaging 2008;26:133–141 doi: 10.1016/j.mri.2007.04.010.
"""

function softwareCoilCompression(kspace::Array{T,6}, numVC::Int64 = size(kspace, 4),sContr::Int = 1,sRep::Int = 1) where T <: Complex
    sx,sy,sz,nChan,nContr,nRep = size(kspace)
    kspaceCC = zeros(T,sx*sy*sz,numVC,nContr,nRep)

    kspaceCC[:,:,sContr,sRep],ccMat = softwareCoilCompression(reshape(kspace[:,:,:,:,sContr,sRep],:,nChan),numVC)
    for rep in 1:nRep, contr in 1:nContr
    if(rep == sRep) && (contr == sContr) # already performed above
        continue
    end
        kspaceCC[:,:,sContr,sRep],ccMat = softwareCoilCompression(reshape(kspace[:,:,:,:,contr,rep],:,nChan),numVC) * ccMat
    end
    kspaceCC = reshape(kspaceCC,sx,sy,sz,numVC,nContr,nRep)
    return kspaceCC, ccMat
end

"""
softwareCoilCompression(acqData::AcquisitionData{T,D}, numVC::Int64 = size(acqData.kdata[1,1,1],2);contr_::Int = 1,rep_::Int = 1) where {T,D}

Perform SVD-based coil compression for the given `acqData`
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
