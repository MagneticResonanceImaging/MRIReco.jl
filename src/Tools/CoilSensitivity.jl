export estimateCoilSensitivities, mergeChannels, espirit, estimateCoilSensitivitiesFixedPoint, geometricCC_2d
using Polyester
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
    `maps = espirit(acqData::AcquisitionData, ksize::NTuple{2,Int64}, ncalib::Int64
               ; eigThresh_1::Float64=0.02, eigThresh_2::Float64=0.95)`

obtains coil sensitivities from a calibration area using ESPIRiT
adapted from the MATLAB code by Uecker et al. for the paper'
M. Uecker, P. Lai, MJ Murphy, P. Virtue, M Elad, JM Pauly, SS Vasanawala and M Lustig, "ESPIRiT- an
eigenvalue approach to autocalibrating parallel MRI: Where SENSE meets GRAPPA", Magn Reson Med, 2013

the matlab code can be found at: [http://people.eecs.berkeley.edu/~mlustig/Software.html]

# Arguments
* `acqData::AcquisitionData`  - AcquisitionData
* `ksize::NTuple{2,Int64}`    - size of the k-space kernel
* `ncalib::Int64`             - number of calibration points in each dimension
* `eigThresh_1::Float64=0.02` - threshold for the singular values of the calibration matrix (relative to the largest value)
* `eigThresh_2::Float64=0.95` - threshold of the image space kernels (if no singular value > `eigThresh_2` exists)
                                , the corresponding pixel has a sensitivity of 0.

Returns a 4D array.
"""
function espirit(acqData::AcquisitionData, ksize::NTuple{2,Int64}, ncalib::Int
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


"""
`maps = espirit(calibData::Array{ComplexF64,3}, imsize::NTuple{2,Int64}, ksize::NTuple{2,Int64}
                ; eigThresh_1::Float64=0.02, eigThresh_2::Float64=0.95)`
Returns a 3D array.
"""
function espirit(calibData::Array{T}, imsize::NTuple{N,Int64}, ksize::NTuple{N,Int64}
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
  BLAS.set_num_threads(Threads.nthreads())
  usv = svd(A)
  BLAS.set_num_threads(nblas)

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
  flip_nc_nv = [1:length(ksize); length(ksize) + 2; length(ksize) + 1]

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
  kernel = permutedims(kernel, flip_nc_nv)

  kern2 = zeros(T, imsize..., nc, nv)
  kern2[CartesianIndices(kernel)] .= kernel
  @views fft!(kern2, 1:length(ksize))
  kern2 ./= sqrt(prod(ksize))
  kern2 .= conj.(kern2)

  eigenVecs = Array{T}(undef, imsize..., nc, nmaps)
  eigenVals = Array{T}(undef, imsize...,  1, nmaps)

  nblas = BLAS.get_num_threads()
  BLAS.set_num_threads(1)
  @batch for n ∈ CartesianIndices(imsize)
    mtx = @view kern2[n, :, :]

    cdv = svd(mtx)
    ph = transpose(exp.(-1im * angle.(cdv.U[1, :])))

    @views eigenVals[n, 1, :] .= real.(cdv.S[1:nmaps])
    D = usv.V * (cdv.U .* ph)
    @views eigenVecs[n, :, :] .= D[:,1:nmaps]
  end
  BLAS.set_num_threads(nblas)

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
    for rep = 1:numRepititions(acqData), contr = 1:numContrasts(acqData)
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
