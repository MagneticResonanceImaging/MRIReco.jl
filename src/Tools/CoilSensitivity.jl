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
  * `use_poweriterations = true` - flag to determine if power iterations are used; power iterations are only used if `nmaps == 1`. They provide speed benefits over the full eigen decomposition, but are an approximation.
"""
function espirit(acqData::AcquisitionData{T}, ksize::NTuple{2,Int} = (6,6), ncalib::Int = 24
  ; eigThresh_1::Number = 0.02, eigThresh_2::Number = 0.95, nmaps::Int = 1, use_poweriterations::Bool = true) where T

  if !isCartesian(trajectory(acqData, 1))
    @error "espirit does not yet support non-cartesian sampling"
  end

  nx, ny = acqData.encodingSize[1:2]
  numChan, numSl = numChannels(acqData), numSlices(acqData)
  maps = zeros(Complex{T}, acqData.encodingSize[1], acqData.encodingSize[2], numSl, numChan, nmaps)

  for slice = 1:numSl
    # form zeropadded array with kspace data
    kdata = zeros(Complex{T}, nx * ny, numChan)
    for coil = 1:numChan
      kdata[acqData.subsampleIndices[1], coil] .= kData(acqData, 1, coil, slice)
    end
    kdata = reshape(kdata, nx, ny, numChan)

    calib = crop(kdata, (ncalib, ncalib, numChan))

    maps[:, :, slice, :, :] .= espirit(calib, (nx, ny), ksize, eigThresh_1 = eigThresh_1, eigThresh_2 = eigThresh_2, nmaps = nmaps, use_poweriterations = use_poweriterations)
  end

  if nmaps == 1
    maps = maps[:,:,:,:,1]
  end
  return maps
end


function espirit(calibData::Array{T}, imsize::NTuple{N,Int}, ksize::NTuple{N,Int} = ntuple(_ -> 6, N)
  ; eigThresh_1::Number = 0.02, eigThresh_2::Number = 0.95, nmaps = 1, use_poweriterations = true) where {T,N}
  nc = size(calibData)[end]

  # compute calibration matrix, perform SVD and convert right singular vectors
  # into k-space kernels
  k, S = dat2Kernel(calibData, ksize)
  idx = findlast(x -> x >= S[1] * eigThresh_1, S) # last singular vector to keep

  # crop kernels and compute eigenvalue decomposition in images space
  # to get sensitivity maps
  maps, W = kernelEig(k[CartesianIndices((ksize..., nc)), 1:idx], imsize, nmaps; use_poweriterations = use_poweriterations)

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

  A = im2row(data, ksize)
  A = reshape(A, size(A,1), :)

  nblas = BLAS.get_num_threads()
  usv = try
    BLAS.set_num_threads(Threads.nthreads())
    svd!(A)
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
function kernelEig(kernel::Array{T}, imsize::Tuple, nmaps=1; use_poweriterations=true) where {T}

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

  usv = svd(kernel, full = (size(kernel, 1) < size(kernel, 2)))

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
    @floop for ix ∈ CartesianIndices(imsize), j ∈ axes(kernel_i,1), i ∈ axes(kernel_i,1)
      kern2[i,j,ix] += kernel_i[i,ix] * conj(kernel_i[j,ix])
    end
  end

  kern2 ./= prod(ksize)
  kern2 .= conj.(kern2)

  eigenVecs = Array{T}(undef, imsize..., nc, nmaps)
  eigenVals = Array{T}(undef, imsize...,  1, nmaps)

  nblas = BLAS.get_num_threads()
  try
    BLAS.set_num_threads(1)
    b    = [randn(T, nc)  for _=1:Threads.nthreads()]
    bᵒˡᵈ = [similar(b[1]) for _=1:Threads.nthreads()]

    @floop for n ∈ CartesianIndices(imsize)
      if use_poweriterations && nmaps==1
        S, U = power_iterations!(view(kern2,:, :, n), b=b[Threads.threadid()], bᵒˡᵈ=bᵒˡᵈ[Threads.threadid()])
        U .*= transpose(exp.(-1im .* angle.(U[1])))
        @views eigenVals[n, 1, :] .= real.(S)
      else
        S, U = eigen!(@view kern2[:, :, n]; sortby = (λ) -> -abs(λ))
        U = @view U[:,1:nmaps]
        U .*= transpose(exp.(-1im .* angle.(@view U[1, :])))
        @views eigenVals[n, 1, :] .= real.(S[1:nmaps])
      end

      @views mul!(eigenVecs[n, :, :], usv.V, U)
    end
  finally
    BLAS.set_num_threads(nblas)
  end

  return eigenVecs, eigenVals
end


"""
    power_iterations!(A;
    b = randn(eltype(A), size(A,2)),
    bᵒˡᵈ = similar(b),
    rtol=1e-6, maxiter=1000, verbose=false)

Power iterations to determine the maximum eigenvalue of a normal operator or square matrix.

# Arguments
* `A`                   - operator or matrix; has to be square

# Keyword Arguments
* `rtol=1e-6`           - relative tolerance; function terminates if the change of the max. eigenvalue is smaller than this values
* `b = randn(eltype(A), size(A,2))` - pre-allocated random vector; will be modified
* `bᵒˡᵈ = similar(b)`   - pre-allocated vector; will be modified
* `maxiter=1000`        - maximum number of power iterations
* `verbose=false`       - print maximum eigenvalue if `true`

# Output
maximum eigenvalue of the operator
"""
function power_iterations!(A;
  b = randn(eltype(A), size(A,2)),
  bᵒˡᵈ = similar(b),
  rtol=1e-6, maxiter=1000, verbose=false)

  b ./= norm(b)
  λ = Inf

  for i = 1:maxiter
    copy!(bᵒˡᵈ, b)
    mul!(b, A, bᵒˡᵈ)

    λᵒˡᵈ = λ
    λ = (bᵒˡᵈ' * b) / (bᵒˡᵈ' * bᵒˡᵈ)
    b ./= norm(b)

    verbose && println("iter = $i; λ = $λ")
    abs(λ/λᵒˡᵈ - 1) < rtol && return λ, b
  end

  return λ, b
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
function geometricCCMat(kdata::Matrix{T}, numVC::Int64 = size(kdata, 2)) where T <: Complex
  return svd(kdata).Vt[:, 1:numVC]
end

"""
perform SVD-based coil compression for the given `kdata` and `smaps`
"""
function geometricCC(kdata::Matrix{T}, smaps::Array{T,4}, numVC::Int64 = size(kdata, 2)) where T <: Complex
  nx, ny, nz, nc = size(smaps)
  usv = svd(kdata)
  kdataCC = kdata * usv.Vt[:, 1:numVC]
  smapsCC = zeros(T, nx, ny, nz, numVC)
  for j = 1:nz, i = 1:ny
    smapsCC[:, i, j, :] .= smaps[:, i, j, :] * usv.Vt[:, 1:numVC]
  end

  return kdataCC, smapsCC
end

"""
perform SVD-based coil compression for the given 2d-encoded `acqData` and `smaps`
"""
function geometricCC_2d(acqData::AcquisitionData{T}, smaps::Array{Complex{T},4}, numVC::Int64 = size(smaps, 4)) where T 
  nx, ny, nz, nc = size(smaps)
  acqDataCC = deepcopy(acqData)
  smapsCC = zeros(Complex{T}, nx, ny, nz, numVC)
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
