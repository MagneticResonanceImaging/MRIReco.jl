export estimateCoilSensitivities, mergeChannels, espirit, estimateCoilSensitivitiesFixedPoint, geometricCC_2d

"""
    `s = estimateCoilSensitivities(I::AbstractArray{T,5})`

Estimates the coil sensitivity based on a reconstruction where the data
from each coil has been reconstructed individually.
Returns a 5D array.
"""
function estimateCoilSensitivities(I::AbstractArray{T,5}, thresh = 1.e-2) where {T}
  nx, ny, nz, ne, numChan = size(I)

  I_sum = sqrt.(sum(abs.(I) .^ 2, dims = 5)) .+ eps()
  I_max = maximum(abs.(I_sum))
  msk = zeros(size(I_sum))
  msk[findall(x -> x > thresh * I_max, I_sum)] .= 1

  s = zeros(eltype(I), size(I))
  for i = 1:numChan
    s[:, :, :, :, i] = msk .* I[:, :, :, :, i] ./ I_sum
  end

  return s
end


"""
    `I4 = mergeChannels(I::AbstractArray{T,5})`

Merge the channels of a multi-coil reconstruction.
Returns a 4D array.
"""
mergeChannels(I::AbstractArray{T,5}) where {T} = sqrt.(sum(abs.(I) .^ 2, dims = 5))


"""
    espirit(acqData::AcquisitionData, ksize::NTuple{2,Int} = (6,6), ncalib::Int = 24
               ; eigThresh_1::Number=0.02, eigThresh_2::Number=0.95, nmaps = 1)

    espirit(calibData::Array{T}, imsize::NTuple{N,Int}, ksize::NTuple{N,Int} = (6,6[,6])
               ; eigThresh_1::Number = 0.02, eigThresh_2::Number = 0.95, nmaps = 1)

Obtains coil sensitivities from a calibration area using ESPIRiT. The code is adapted from the MATLAB code by Uecker et al. (cf. [Uecker et al. "ESPIRiT—an eigenvalue approach to autocalibrating parallel MRI: Where SENSE meets GRAPPA"](https://doi.org/10.1002/mrm.24751)). The matlab code can be found at: [http://people.eecs.berkeley.edu/~mlustig/Software.html]

# Method 1
  The first method of this function works with 2D/multi-slice data in the MRIReco.jl data format:

    espirit(acqData::AcquisitionData, ksize::NTuple{2,Int} = (6,6), ncalib::Int = 24
              ; eigThresh_1::Number=0.02, eigThresh_2::Number=0.95)

## Required Arguments
* `acqData::AcquisitionData`  - AcquisitionData

## Optional Arguments
* `ksize::NTuple{2,Int64}`    - size of the k-space kernel; `default = (6,6)`
* `ncalib::Int64`             - number of calibration points in each dimension; `default = 30`

## Keyword Arguments
  * `eigThresh_1::Number=0.02`  - threshold for the singular values of the calibration matrix (relative to the largest value); reduce for more accuracy, increase for saving memory and computation time.
  * `eigThresh_2::Number=0.95`  - threshold to mask the final maps: for each voxel, the map will be set to 0, if, for this voxel, no singular value > `eigThresh_2` exists.
  * `nmaps = 1`                 - Number of maps that are calcualted. Set to 1 for regular SENSE; set to 2 for soft-SENSE (cf. [Uecker et al. "ESPIRiT—an eigenvalue approach to autocalibrating parallel MRI: Where SENSE meets GRAPPA"](https://doi.org/10.1002/mrm.24751)).

# Method 2
  The second method of this function works with single slice 2D or 3D data; the first argument is the calibration data, i.e. the cropped center of k-space:

    espirit(calibData::Array{T}, imsize::NTuple{N,Int}, ksize::NTuple{N,Int} = (6,6[,6])
              ; eigThresh_1::Number = 0.02, eigThresh_2::Number = 0.95, nmaps = 1)


# Required Arguments
  * `calibData::Array{T}`       - Center of k-space in the format `kx × ky [× kz] × coils`, where the kz dimension is optional. The type `T` of the input data determines the type of the calculated maps. Reasonable choises are in the range of `Nx = Ny [= Nz] = 24`.
  * `imsize::NTuple{N,Int}`     - matrix size of the final maps

# Optional Arguments
  * `ksize::NTuple{N,Int}`      - number of calibration points in each dimension; the default is `(6,6)` for 2D and `(6,6,6)` for 3D.

# Keyword Arguments
  * `eigThresh_1::Number=0.02`  - threshold for the singular values of the calibration matrix (relative to the largest value); reduce for more accuracy, increase for saving memory and computation time.
  * `eigThresh_2::Number=0.95`  - threshold to mask the final maps: for each voxel, the map will be set to 0, if, for this voxel, no singular value > `eigThresh_2` exists.
  * `nmaps = 1`                 - Number of maps that are calcualted. Set to 1 for regular SENSE; set to 2 for soft-SENSE (cf. [Uecker et al. "ESPIRiT—an eigenvalue approach to autocalibrating parallel MRI: Where SENSE meets GRAPPA"](https://doi.org/10.1002/mrm.24751)).
"""
function espirit(acqData::AcquisitionData, ksize::NTuple{2,Int} = (6,6), ncalib::Int = 24
  ; eigThresh_1::Number = 0.02, eigThresh_2::Number = 0.95, nmaps::Int = 1)

  if !isCartesian(trajectory(acqData, 1))
    @error "espirit does not yet support non-cartesian sampling"
  end

  nx, ny = acqData.encodingSize[1:2]
  numChan, numSl = numChannels(acqData), numSlices(acqData)
  maps = zeros(ComplexF64, acqData.encodingSize[1], acqData.encodingSize[2], numSl, numChan, nmaps)

  for slice = 1:numSl
    # form zeropadded array with kspace data
    kdata = zeros(ComplexF64, nx * ny, numChan)
    for coil = 1:numChan
      kdata[acqData.subsampleIndices[1], coil] .= kData(acqData, 1, coil, slice)
    end
    kdata = reshape(kdata, nx, ny, numChan)

    calib = crop(kdata, (ncalib, ncalib, numChan))

    maps[:, :, slice, :, :] .= espirit(calib, (nx, ny), ksize, eigThresh_1 = eigThresh_1, eigThresh_2 = eigThresh_2, nmaps = nmaps)
  end

  if nmaps == 1
    maps = maps[:,:,:,:,1]
  end
  return maps
end


function espirit(calibData::Array{T}, imsize::NTuple{N,Int}, ksize::NTuple{N,Int} = ntuple(_ -> 6, N)
  ; eigThresh_1::Number = 0.02, eigThresh_2::Number = 0.95, nmaps = 1) where {T,N}
  nc = size(calibData)[end]

  # compute calibration matrix, perform SVD and convert right singular vectors
  # into k-space kernels
  k, S = dat2Kernel(calibData, ksize)
  idx = findlast(x -> x >= S[1] * eigThresh_1, S) # last singular vector to keep

  # crop kernels and compute eigenvalue decomposition in images space
  # to get sensitivity maps
  maps, W = kernelEig(k[CartesianIndices((ksize..., nc)), 1:idx], imsize, nmaps)

  # sensitivity maps correspond to eigen vectors with a singular value of 1
  msk = falses(size(W))
  msk[findall(x -> x > eigThresh_2, abs.(W))] .= true
  maps .*= msk
  maps = fftshift(maps, 1:length(imsize))

  return maps
end

#
# computes ESPIRiT step I:
#   form calibration matrix from data, perform SVD and convert right singular vectors
#   into k-space kernels
#
function dat2Kernel(data::Array{T,M}, ksize::NTuple{N,Int64}) where {T,N,M}
  M != N+1 && error("data must have 1 dimension more than the length of ksize")

  nc = size(data)[end]

  tmp = im2row(data, ksize)
  tsx, tsy, tsz = size(tmp)
  A = reshape(tmp, tsx, tsy * tsz)

  nblas = BLAS.get_num_threads()
  usv = try
    BLAS.set_num_threads(Threads.nthreads())
    svd(A)
  finally
    BLAS.set_num_threads(nblas)
  end

  kernel = reshape(usv.V, ksize..., nc, size(usv.V, 2))
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
#         eigenvals: images representing the eigenvalues (sx,sy,numChan)
function kernelEig(kernel::Array{T}, imsize::Tuple, nmaps=1) where {T}

  kern_size = size(kernel)
  ksize = kern_size[1:end-2]
  nc = kern_size[end-1]
  nv = kern_size[end]
  flip_nc_nv  = [1:length(ksize); length(ksize) + 2; length(ksize) + 1]
  flip_nc_nv2 = [length(ksize) + 2; length(ksize) + 1; 1:length(ksize)]

  nmaps = min(nc, nv, nmaps)

  # rotate kernel to order by maximum variance
  kernel = permutedims(kernel, flip_nc_nv)
  kernel = reshape(kernel, :, nc)

  if size(kernel, 1) < size(kernel, 2)
    usv = svd(kernel, full = true)
  else
    usv = svd(kernel)
  end

  kernel = kernel * usv.V
  kernel = reshape(kernel, ksize..., nv, nc)
  kernel = permutedims(kernel, flip_nc_nv2)

  kern2 = zeros(T, nc, nc, imsize...)
  kernel_k = zeros(T, nc, imsize...)
  kernel_i = similar(kernel_k)
  fftplan = plan_fft(kernel_i, 2:length(ksize)+1; flags=FFTW.MEASURE)

  for iv ∈ axes(kernel, 2)
    @views kernel_k[:,CartesianIndices(ksize)] .= kernel[:,iv,CartesianIndices(ksize)]
    mul!(kernel_i, fftplan, kernel_k)
    for ix ∈ CartesianIndices(imsize)
      @views mul!(kern2[:,:,ix], kernel_i[:,ix], kernel_i[:,ix]', 1, 1)
    end
  end

  kern2 ./= sqrt(prod(ksize))
  kern2 .= conj.(kern2)

  eigenVecs = Array{T}(undef, imsize..., nc, nmaps)
  eigenVals = Array{T}(undef, imsize...,  1, nmaps)

  nblas = BLAS.get_num_threads()
  try
    BLAS.set_num_threads(1)
    Threads.@threads for n ∈ CartesianIndices(imsize)
      S, U = eigen!(@view kern2[:, :, n]; sortby = (λ) -> -abs(λ))

      @views eigenVals[n, 1, :] .= real.(S[1:nmaps])

      U = @view U[:,1:nmaps]
      U .*= transpose(exp.(-1im .* angle.(@view U[1, :])))
      @views mul!(eigenVecs[n, :, :], usv.V, U)
    end
  finally
    BLAS.set_num_threads(nblas)
  end

  return eigenVecs, eigenVals
end

#
# rearrange data from the image data to form kernels for dat2Kernel
#
function im2row(img::Array{T,M}, winSize::NTuple{N,Int}) where {N,M,T}

  simg = size(img)
  sk = simg[1:end-1]
  ncoils = simg[end]

  res = Array{T}(undef, prod(sk .- winSize .+ 1), prod(winSize), ncoils)

  m = CartesianIndices(sk .- winSize .+ 1)
  m = m .- CartesianIndex(ntuple(_ -> 1, N))
  cnt = 0
  for ms ∈ CartesianIndices(winSize)
    cnt += 1
    res[:, cnt, :] .= reshape(img[m .+ ms, :], :, ncoils)
  end
  return res
end

#
# crop the central area (of size s) of an array A
#
function crop(A::Array{T,3}, s::NTuple{3,Int64}) where {T}
  nx, ny, nz = size(A)
  idx_x = div(nx, 2)-div(s[1], 2)+1:div(nx, 2)-div(s[1], 2)+s[1]
  idx_y = div(ny, 2)-div(s[2], 2)+1:div(ny, 2)-div(s[2], 2)+s[2]
  idx_z = div(nz, 2)-div(s[3], 2)+1:div(nz, 2)-div(s[3], 2)+s[3]
  return A[idx_x, idx_y, idx_z]
end




"""
return SVD-based coil compression matrix for `numVC` virtual coils
"""
function geometricCCMat(kdata::Matrix{ComplexF64}, numVC::Int64 = size(kdata, 2))
  return svd(kdata).Vt[:, 1:numVC]
end

"""
perform SVD-based coil compression for the given `kdata` and `smaps`
"""
function geometricCC(kdata::Matrix{ComplexF64}, smaps::Array{ComplexF64,4}, numVC::Int64 = size(kdata, 2))
  nx, ny, nz, nc = size(smaps)
  usv = svd(kdata)
  kdataCC = kdata * usv.Vt[:, 1:numVC]
  smapsCC = zeros(ComplexF64, nx, ny, nz, numVC)
  for j = 1:nz, i = 1:ny
    smapsCC[:, i, j, :] .= smaps[:, i, j, :] * usv.Vt[:, 1:numVC]
  end

  return kdataCC, smapsCC
end

"""
perform SVD-based coil compression for the given 2d-encoded `acqData` and `smaps`
"""
function geometricCC_2d(acqData::AcquisitionData, smaps::Array{ComplexF64,4}, numVC::Int64 = size(smaps, 4))
  nx, ny, nz, nc = size(smaps)
  acqDataCC = deepcopy(acqData)
  smapsCC = zeros(ComplexF64, nx, ny, nz, numVC)
  for sl = 1:numSlices(acqData)
    # use first echo and first repetition to determine compression matrix
    ccMat = geometricCCMat(acqData.kdata[1, sl, 1], numVC)
    # compress slice of the smaps
    for i = 1:ny
      smapsCC[:, i, sl, :] .= smaps[:, i, sl, :] * ccMat
    end
    # compress kdata-slice
    for rep = 1:numRepetitions(acqData), contr = 1:numContrasts(acqData)
      acqDataCC.kdata[contr, sl, rep] = acqData.kdata[contr, sl, rep] * ccMat
    end
  end
  return acqDataCC, smapsCC
end

function estimateCoilSensitivitiesFixedPoint(acqData::AcquisitionData;
  iterations = 1, outerIterations = 3, cmap = nothing)

  N = acqData.encodingSize

  params = Dict{Symbol,Any}()
  params[:reco] = "standard" # "multiCoil"
  #params[:solver] = "cgnr" #solver
  #params[:iterations] = iterations
  params[:reconSize] = (N, N)
  params[:alpha] = 1.25
  #params[:iterationsInner] = 5
  params[:normalizeReg] = true


  if cmap != nothing
    params[:correctionMap] = cmap
    params[:K] = 15
    params[:m] = 3.0
  end

  Ireco = reconstruction(acqData, params)
  s = estimateCoilSensitivities(IrecoCorr).data

  #if !isempty(coilsens)
  #  params[:senseMaps] = coilsens
  #end



end
