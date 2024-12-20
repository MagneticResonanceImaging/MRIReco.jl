export softwareCoilCompression, applyCoilCompressionSensitivityMaps, geometricCoilCompression

"""
    geometricCoilCompression(kspace::Matrix{T}, numVC::Int64 = size(kdata, 2)) where T <: Complex

return a SVD-based Geometric Coil Compression (GCC) matrix for `numVC` virtual coils.
Coil compression is performed for each position along the readout (kx) direction.

Input :
    -   kspace::Array{T,6} : dimension : kx,ky,kz,channel,contrasts,repetitions
    -   numVS::Int

Optional :
    -   dim = 1 : dimension to perform the FFT and apply the coil compression
    -   sContr = 1 : Contrast to use for calibration
    -   sRep = 1 : repetition to use for calibration

Reference :
- Zhang T, Pauly JM, Vasanawala SS, Lustig M. Coil compression for accelerated imaging with Cartesian sampling. Magnetic Resonance in Medicine 2013;69:571–582 doi: 10.1002/mrm.24267.
"""
function geometricCoilCompression(kspace::Array{T,6}, numVC::Int = size(kdata, 4); dim::Int=1,sContr::Int = 1,sRep::Int = 1) where T <: Complex
    # permute data if compression is done on 2nd or 3rd dimensions
    # done for simple reusible code.

    if dim==2
        kspace = permutedims(kspace,[2,1,3,4,5,6]);
    end

    if dim==3
        kspace = permutedims(kspace,[3,1,2,4,5,6]);
    end

    if dim>3 | dim <1
        @error "Error, compression dimension is 1 2 or 3"
    end

    sx,sy,sz,nChan,nContr,nRep = size(kspace)
    # perform ifft along 1st dimension
    calibData = ifftshift(ifft(ifftshift(kspace),1))
    calibData = reshape(calibData,sx,sy*sz,nChan,nContr,nRep)

    kspaceCC = zeros(T,sx,sy*sz,numVC,nContr,nRep)
    ccMat = zeros(T,sx,nChan,numVC)

    for kx in axes(kspace,1)
        usv = svd(calibData[kx,:,:,sContr,sRep])
        ccMat[kx,:,:]  = usv.V[:, 1:numVC]
    end

    # align ccMat
    ccMat = alignCCMtx(ccMat,numVC)
    # apply ccMat
    for kx in axes(kspace,1)
        for rep in 1:nRep, contr in 1:nContr
            kspaceCC[kx,:,:,contr,rep] = calibData[kx,:,:,contr,rep] * ccMat[kx,:,:]
        end
    end

    kspaceCC = reshape(kspaceCC,sx,sy,sz,numVC,nContr,nRep)
    kspaceCC = fftshift(fft(fftshift(kspaceCC),1))

    if dim==2
        kspaceCC = permutedims(kspaceCC,[2,1,3,4,5,6]);
    end

    if dim==3
        kspaceCC = permutedims(kspaceCC,[3,1,2,4,5,6]);
    end

    return kspaceCC, ccMat
end


function geometricCoilCompression(acq::AcquisitionData{T,D}, numVC::Int = size(acq.kdata,2); dim::Int=1,sContr::Int = 1,sRep::Int = 1) where {T,D}
    kdata = kDataCart(acq)
    D == 2 ? enc2D = true : enc2D = false

    kdataCC,ccMat2 = geometricCoilCompression(kdata, numVC, dim=dim,sContr=sContr,sRep=sRep)
    return AcquisitionData(kdataCC;enc2D),ccMat2
end

"""
     mtx = alignCCmtx(mtx, [ ncc)

Align coil compression matrices based on nearest spanning vectors in
subspaces. This is an implementation based on Zhang et. al MRM 2013;69(2):571-82.

 Inputs:

 Outputs:
           mtx - aligned compression matrices.

"""
function alignCCMtx(ccMat,numVC::Int=size(ccMat,3))
    sx,_ = size(ccMat);

    # align everything based on the middle slice.
    n0 = floor(Int32,sx/2);
    A00 = ccMat[n0,:,1:numVC];

    # Align backwards to first slice
    A0 = copy(A00)
    for n = n0-1:-1:1
        A1 = ccMat[n,:,1:numVC];
        C = A1'*A0;
        usv= svd(C);
        P = usv.V*usv.U';
        ccMat[n,:,1:numVC] = A1*P';
        A0 = ccMat[n,:,1:numVC];
    end

    # Align forward to end slice
    A0 = copy(A00);
    for n = n0+1:sx
        A1 = ccMat[n,:,1:numVC];
        C = A1'*A0;
        usv= svd(C);
        P = usv.V*usv.U';
        ccMat[n,:,1:numVC] = A1*P';
        A0 = ccMat[n,:,1:numVC];
    end
    return ccMat
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
  ccMat = usv.V[:, 1:numVC]
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
        kspaceCC[:,:,sContr,sRep],ccMat = reshape(kspace[:,:,:,:,contr,rep],:,nChan) * ccMat
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
