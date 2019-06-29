export estimateCoilSensitivities, mergeChannels, espirit

"""
    `s = estimateCoilSensitivities(I::AbstractArray{T,5})`

Estimates the coil sensitivity based on a reconstruction where the data
from each coil has been reconstructed individually.
Returns a 5D array.
"""
function estimateCoilSensitivities(I::AbstractArray{T,5}) where T
  numChan = size(I,5)

  I_sum = sqrt.( sum(abs.(I).^2, dims=5) )

  s = similar(I)
  for i=1:numChan
    s[:,:,:,:,i] = I[:,:,:,:,i] ./ I_sum
  end

  return s
end


"""
    `I4 = mergeChannels(I::AbstractArray{T,5})`

Merge the channels of a multi-coil reconstruction.
Returns a 4D array.
"""
mergeChannels(I::AbstractArray{T,5}) where T = sqrt.(sum(abs.(I).^2,dims=5))


"""
    `maps = espirit(acqData::AcquisitionData, ksize::NTuple{2,Int64}, ncalib::Int64
               ; eigThresh_1::Float64=0.02, eigThresh_2::Float64=0.95)`

obtains coil sensitivities from a calibration area using ESPIRiT
adapted from the MATLAB code by Uecker et al. for the paper'
M. Uecker, P. Lai, MJ Murphy, P. Virtue, M Elad, JM Pauly, SS Vasanawala and M Lustig, "ESPIRiT- an
eigenvalue approach to autocalibrating parallel MRI: Where SENSE meets GRAPPA", Magn Reson Med, 2013

the matlab code can be found at: [http://people.eecs.berkeley.edu/~mlustig/Software.html]

# Arguments
* `acqData::AcquisitionData`  - AcquisitionData
* `ksize::NTuple{2,Int64}`    - size of the k-space matrix
* `ncalib::Int64`             - number of calibration points in each dimension
* `eigThresh_1::Float64=0.02` - threshold for the singular values of the calibration matrix (relative to the largest value)
* `eigThresh_2::Float64=0.95` - threshold of the image space kernels (if no singular value > `eigThresh_2` exists)
                                , the corresponding pixel has a sensitivity of 0.

Returns a 4D array.
"""
function espirit(acqData::AcquisitionData, ksize::NTuple{2,Int64}, ncalib::Int64
                ; eigThresh_1::Float64=0.02, eigThresh_2::Float64=0.95)

  if !isCartesian(trajectory(acqData,1))
    @error "espirit does not yet support non-cartesian sampling"
  end

  nx,ny = acqData.encodingSize[1:2]
  numChan, numSl = numChannels(acqData), numSlices(acqData)
  maps = zeros(ComplexF64,acqData.encodingSize[1],acqData.encodingSize[2],numSl,numChan)

  for slice = 1:numSl
    # form zeropadded array with kspace data
    kdata = zeros(ComplexF64,nx*ny,numChan)
    for coil = 1:numChan
      kdata[acqData.subsampleIndices[1],coil] .= kData(acqData,1,coil,slice)
    end
    kdata = reshape(kdata,nx,ny,numChan)

    calib = crop(kdata,(ncalib,ncalib,numChan))

    maps[:,:,slice,:] .= espirit(calib,(nx,ny),ksize,eigThresh_1=eigThresh_1,eigThresh_2=eigThresh_2)
  end
  return maps
end


"""
`maps = espirit(calibData::Array{ComplexF64,3}, imsize::NTuple{2,Int64}, ksize::NTuple{2,Int64}
                ; eigThresh_1::Float64=0.02, eigThresh_2::Float64=0.95)`
Returns a 3D array.
"""
function espirit(calibData::Array{ComplexF64,3}, imsize::NTuple{2,Int64}, ksize::NTuple{2,Int64}
                ; eigThresh_1::Float64=0.02, eigThresh_2::Float64=0.95)
  sx,sy = imsize
  nc = size(calibData,3)

  # compute calibration matrix, perform SVD and convert right singular vectors
  # into k-space kernels
  k, S = dat2Kernel(calibData,ksize)
  idx = findlast(x->x>=S[1]*eigThresh_1,S) # last singular vector to keep

  # crop kernels and compute eigenvalue decomposition in images space
  # to get sensitivity maps
  M, W = kernelEig(k[:,:,:,1:idx],(sx,sy))

  # sensitivity maps correspond to eigen vectors with a singular value of 1
  msk = zeros(size(W,1),size(W,2),1)
  msk[findall(x->x>eigThresh_2, abs.(W[:,:,end]))] .= 1
  maps = M[:,:,:,end] .* repeat(msk,1,1,nc)

  return maps
end

#
# computes ESPIRiT step I:
#   form calibration matrix from data, perform SVD and convert right singular vectors
#   into k-space kernels
#
function dat2Kernel(data::Array{T,3}, ksize::NTuple{2,Int64}) where T
  sx,sy,nc = size(data)
  imSize=(sx,sy)

  tmp = im2row(data,ksize)
  tsx,tsy,tsz = size(tmp)
  A = reshape(tmp, tsx, tsy*tsz)
  usv = svd(A)
  kernel = reshape(usv.V,ksize[1],ksize[2],nc,size(usv.V,2))
  return kernel, usv.S
end

# computes ESPIRiT step II:
#   eigenvalue decomposition of a k-space kernel in image space.
#
# inputs:
#         kernel - kspace kernels (4D)
#         imSize - size of the image to compute maps
#
# outputs :
#         eigenvecs: images representing the eigenvectors (sx,sy,numChan,numChan)
#         eigenvals: images representing th eigenvalues (sx,sy,numChan)
function kernelEig(kernel::Array{T,4},imsize::NTuple{2,Int64}) where T
  nx,ny,nc,nv = size(kernel)
  ksize = (nx,ny)

  # rotate kernel to order by maximum variance
  k = permutedims(kernel,[1,2,4,3])
  k = reshape(k,nx*ny*nv,nc)

  if size(k,1) < size(k,2)
    usv = svd(k,full=true)
  else
    usv = svd(k)
  end

  k = k*usv.V
  kernel = reshape(k, nx,ny,nv,nc)
  kernel = permutedims(kernel,[1,2,4,3])
  kern2 = zeros(ComplexF64, imsize[1],imsize[2],size(kernel,3),size(kernel,4))
  for n = 1:size(kernel,4)
    kern2[:,:,:,n] .= fft2c( zeropad( conj.(kernel[end:-1:1,end:-1:1,:,n])*sqrt(imsize[1]*imsize[2]), (imsize[1],imsize[2],size(kernel,3)) ) )
  end
  kern2 = kern2/sqrt(prod(ksize))

  eigenVecs = zeros(ComplexF64,imsize[1],imsize[2],nc,min(nc,nv))
  eigenVals = zeros(ComplexF64,imsize[1],imsize[2],min(nc,nv))

  for n = 1:prod(imsize)
    x,y = Tuple( CartesianIndices((imsize[1],imsize[2]))[n] )
    mtx = kern2[x,y,:,:]

    cdv = svd(mtx)
    ph = transpose( repeat( exp.(-1im*angle.(cdv.U[1,:])), 1,size(cdv.U,1) ) )
    C = usv.V*(cdv.U .* ph)
    D = real.(cdv.S)
    eigenVals[x,y,:] .= D[end:-1:1]
    eigenVecs[x,y,:,:] .= C[:,end:-1:1]
  end

  return eigenVecs, eigenVals
end

#
# rearrange data from the image data to form kernels for dat2Kernel
#
function im2row(img::Array{T,3}, winSize::NTuple{2,Int64}) where T
  sx,sy,sz = size(img)
  res = zeros(T,(sx-winSize[1]+1)*(sy-winSize[2]+1),prod(winSize),sz)
  cnt = 0
  for y = 1:winSize[2]
    for x = 1:winSize[1]
      cnt += 1
      res[:,cnt,:] .= reshape(img[x:sx-winSize[1]+x,y:sy-winSize[2]+y,:]
                              ,(sx-winSize[1]+1)*(sy-winSize[2]+1),sz)
    end
  end
  return res
end

#
# crop the central area (of size s) of an array A
#
function crop(A::Array{T,3},s::NTuple{3,Int64}) where T
  nx,ny,nz = size(A)
  idx_x = div(nx,2)-div(s[1],2)+1:div(nx,2)-div(s[1],2)+s[1]
  idx_y = div(ny,2)-div(s[2],2)+1:div(ny,2)-div(s[2],2)+s[2]
  idx_z = div(nz,2)-div(s[3],2)+1:div(nz,2)-div(s[3],2)+s[3]
  return A[idx_x,idx_y,idx_z]
end

#
# zeropoad an array A to a target size s
#
function zeropad(A::Array{T,3},s::NTuple{3,Int64}) where T
  nx,ny,nz = size(A)
  idx_x = div(s[1],2)-div(nx,2)+1:div(s[1],2)-div(nx,2)+nx
  idx_y = div(s[2],2)-div(ny,2)+1:div(s[2],2)-div(ny,2)+ny
  idx_z = div(s[3],2)-div(nz,2)+1:div(s[3],2)-div(nz,2)+nz

  res = zeros(T,s[1],s[2],s[3])
  res[idx_x,idx_y,idx_z] .= A

  return res
end

#
#  apply 2d fft to all slices of a higher dimensional array
#
function fft2c(x::Array{T}) where T
  s = size(x)
  x = reshape(x,s[1],s[2],:)
  res = zeros(ComplexF64,size(x))
  for n=1:size(x,3)
    res[:,:,n] = 1/sqrt(s[1]*s[2])*fftshift(fft(ifftshift(x[:,:,n])))
  end
  return reshape(res,s)
end

#
#  apply 2d ifft to all slices of a higher dimensional array
#
function ifft2c(x::Array{T}) where T
  s = size(x)
  x = reshape(x,s[1],s[2],:)
  res = zeros(ComplexF64,size(x))
  for n=1:size(x,3)
    res[:,:,n] = sqrt(s[1]*s[2])*fftshift(fft(ifftshift(x[:,:,n])))
  end
  return reshape(res,s)
end
