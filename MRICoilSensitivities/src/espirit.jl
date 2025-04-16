"""
    espirit(acqData::AcquisitionData, ksize::NTuple{D,Int} = (6,6[,6]), ncalib::Int = 24,
               imsize::NTuple{D,Int}=Tuple(acqData.encodingSize[1:D])
               ; eigThresh_1::Number=0.02, eigThresh_2::Number=0.95, nmaps = 1) where D

    espirit(calibData::Array{T}, imsize::NTuple{D,Int}, ksize::NTuple{S,Int} = (6,6[,6])
               ; eigThresh_1::Number = 0.02, eigThresh_2::Number = 0.95, nmaps = 1) where D

Obtains coil sensitivities from a calibration area using ESPIRiT. Works with 2D and 3D data.

The code is adapted from the MATLAB code by Uecker et al. (cf. [Uecker et al. "ESPIRiT—an eigenvalue approach to autocalibrating parallel MRI:
Where SENSE meets GRAPPA"](https://doi.org/10.1002/mrm.24751)). The matlab code can be found at:
[http://people.eecs.berkeley.edu/~mlustig/Software.html]

# Required Arguments
  * `acqData::AcquisitionData`  - AcquisitionData
  or
  * `calibData::Array{T}`       - Center of k-space in the format `kx × ky [× kz] × coils`, where the kz dimension is optional. The type `T` of the input data determines the type of the calculated maps. Reasonable choises are in the range of `Nx = Ny [= Nz] = 24`.
  * `imsize::NTuple{D,Int}`     - matrix size of the final maps (optional when passing acqData)

# Optional Arguments
  * `ksize::NTuple{D,Int}`      - number of calibration points in each dimension; the default is `(6,6)` for 2D and `(6,6,6)` for 3D.
  * `ncalib::Int`      - size of the central part corresponding to the calibration data.

# Keyword Arguments
  * `eigThresh_1::Number=0.02`  - threshold for the singular values of the calibration matrix (relative to the largest value); reduce for more accuracy, increase for saving memory and computation time.
  * `eigThresh_2::Number=0.95`  - threshold to mask the final maps: for each voxel, the map will be set to 0, if, for this voxel, no singular value > `eigThresh_2` exists.
  * `nmaps = 1`                 - Number of maps that are calculated. Set to 1 for regular SENSE; set to 2 for soft-SENSE (cf. [Uecker et al. "ESPIRiT—an eigenvalue approach to autocalibrating parallel MRI: Where SENSE meets GRAPPA"](https://doi.org/10.1002/mrm.24751)).
  * `use_poweriterations = true` - flag to determine if power iterations are used; power iterations are only used if `nmaps == 1`. They provide speed benefits over the full eigen decomposition, but are an approximation.
  * `echo::Int = 1` - select the k-space corresponding to the echo
  * `rep::Int = 1` - select the k-space corresponding to the rep
"""
function espirit(acqData::AcquisitionData{T,D}, ksize::NTuple{D,Int}=ntuple(d -> 6, D), ncalib::Int=24,
  imsize::NTuple{D,Int}=encodingSize(acqData);
  nmaps::Int=1,
  echo::Int=1, rep::Int=1, kargs...) where {T,D}

  if !isCartesian(trajectory(acqData, 1))
    @error "espirit does not yet support non-cartesian sampling"
  end

  acqsize = encodingSize(acqData)
  match_acq_size = all(imsize .== acqsize)

  #  Force maps to be at least same size as acqData.encodingSize.
  if any(acqsize .> imsize)
    imsize = acqsize
  end

  idx = match_acq_size ? acqData.subsampleIndices[1] : findIndices(imsize, acqsize)[acqData.subsampleIndices[1]]

  numChan, numSl = numChannels(acqData), numSlices(acqData)
  maps = zeros(Complex{T}, imsize..., numSl, numChan, nmaps)

  for slice = 1:numSl
    # form zeropadded array with kspace data
    kdata = zeros(Complex{T}, prod(imsize), numChan)
    for coil = 1:numChan
      kdata[idx, coil] .= kData(acqData, echo, coil, slice; rep=rep)
    end
    kdata = reshape(kdata, imsize..., numChan)

    calib = crop(kdata, ntuple(d -> ncalib, D))

    maps[CartesianIndices(imsize), slice, :, :] .= espirit(calib, imsize, ksize; kargs...)
  end

  if D == 3 # in case of 3D measurements, we have no slices
    maps = dropdims(maps, dims=4)
  end

  if nmaps == 1
    maps = maps[:, :, :, :, 1]
  end
  return maps
end

function espirit(calibData::Array{T}, imsize::NTuple{D,Int}, ksize::NTuple{D,Int}=ntuple(_ -> 6, D)
  ; eigThresh_1::Number=0.02, eigThresh_2::Number=0.95, nmaps=1, use_poweriterations=true) where {T,D}
  nc = size(calibData)[end]

  # compute calibration matrix, perform SVD and convert right singular vectors
  # into k-space kernels
  k, S = dat2Kernel(calibData, ksize)
  idx = findlast(x -> x >= S[1] * eigThresh_1, S) # last singular vector to keep

  # crop kernels and compute eigenvalue decomposition in images space
  # to get sensitivity maps
  maps, W = kernelEig(k[CartesianIndices((ksize..., nc)), 1:idx], imsize, nmaps; use_poweriterations=use_poweriterations)

  # sensitivity maps correspond to eigen vectors with a singular value of 1
  msk = zeros(Bool, size(W))
  msk[W.>eigThresh_2] .= true
  maps .*= msk
  maps_ = fftshift(maps, 1:length(imsize))

  return maps_
end

#
# computes ESPIRiT step I:
#   form calibration matrix from data, perform SVD and convert right singular vectors
#   into k-space kernels
#
function dat2Kernel(data::Array{T,M}, ksize::NTuple{D,Int64}) where {T,D,M}
  M != D + 1 && error("data must have 1 dimension more than the length of ksize")

  nc = size(data)[end]

  A = im2row(data, ksize)
  A = reshape(A, size(A, 1), :)

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
  sizePadded = 2 .* ksize
  nc = kern_size[end-1]
  nv = kern_size[end]
  flip_nc_nv = [1:length(ksize); length(ksize) + 2; length(ksize) + 1]

  nmaps = min(nc, nv, nmaps)

  # rotate kernel to order by maximum variance
  kernel = permutedims(kernel, flip_nc_nv)
  kernel = reshape(kernel, :, nc)

  usv = svd(kernel, full=(size(kernel, 1) < size(kernel, 2)))

  kernel = kernel * usv.V
  kernel = reshape(kernel, ksize..., nv, nc)

  kern2_ = zeros(T, nc, nc, sizePadded...)
  kernel_k = zeros(T, sizePadded..., nc)
  kernel_i = similar(kernel_k)

  fftplan = plan_fft(kernel_i, 1:ndims(kernel_i)-1; flags=FFTW.MEASURE, num_threads=Threads.nthreads())

  for iv ∈ axes(kernel, length(kern_size) - 1)
    @views kernel_k[CartesianIndices(ksize), :] .= kernel[CartesianIndices(ksize), iv, :]
    mul!(kernel_i, fftplan, kernel_k)

    @floop for j ∈ axes(kernel_i, ndims(kernel_i)), ix ∈ CartesianIndices(sizePadded)
      @inbounds @simd for i ∈ axes(kernel_i, ndims(kernel_i))
        kern2_[i, j, ix] += conj(kernel_i[ix, i]) * kernel_i[ix, j]
      end
    end
  end
  # resize is necessary for compatibility with MKL
  reshape(fft!(reshape(kern2_, nc^2, sizePadded...), 2:ndims(kern2_)-1), nc, nc, sizePadded...)

  kern2 = zeros(T, nc, nc, imsize...)

  # This code works but requires one additional allocation of kern2
  # idx_center = CartesianIndices(sizePadded) .- CartesianIndex(sizePadded .÷ 2) .+ CartesianIndex(imsize .÷ 2)
  # @views ifftshift!(kern2[:,:,idx_center], kern2_, 3:ndims(kern2_))
  # kern2 = ifftshift(kern2, 3:ndims(kern2_))

  idx_center = collect(CartesianIndices(sizePadded) .- CartesianIndex(sizePadded .÷ 2) .+ CartesianIndex(imsize .÷ 2))
  [idx_center[i] = CartesianIndex(mod1.(Tuple(idx_center[i]) .+ cld.(imsize, 2), imsize)) for i in eachindex(idx_center)]
  @views ifftshift!(kern2[:, :, idx_center], kern2_, 3:ndims(kern2_))

  reshape(ifft!(reshape(kern2, nc^2, imsize...), 2:ndims(kern2)-1), nc, nc, imsize...)
  kern2 .*= prod(imsize) / prod(sizePadded) / prod(ksize)


  eigenVecs = Array{T}(undef, imsize..., nc, nmaps)
  eigenVals = Array{real(T)}(undef, imsize..., 1, nmaps)

  nblas = BLAS.get_num_threads()
  try
    BLAS.set_num_threads(1)
    b = [randn(T, nc) for _ = 1:Threads.nthreads()]
    bᵒˡᵈ = [similar(b[1]) for _ = 1:Threads.nthreads()]

    @floop for n ∈ CartesianIndices(imsize)
      if use_poweriterations && nmaps == 1
        S, U = power_iterations!(view(kern2, :, :, n), b=b[Threads.threadid()], bᵒˡᵈ=bᵒˡᵈ[Threads.threadid()])
        # The following uses the method from IterativeSolvers.jl but is currently slower and allocating more
        # this is probably because we pre-allocate bᵒˡᵈ
        #S, U = RegularizedLeastSquares.IterativeSolvers.powm!(view(kern2, :, :, n), b[Threads.threadid()], maxiter=5)

        U .*= transpose(exp.(-1im .* angle.(U[1])))
        @views eigenVals[n, 1, :] .= real.(S)
      else
        S, U = eigen!(view(kern2, :, :, n); sortby=(λ) -> -abs(λ))
        U = @view U[:, 1:nmaps]
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
  b=randn(eltype(A), size(A, 2)),
  bᵒˡᵈ=similar(b),
  rtol=1e-6, maxiter=1000, verbose=false)

  b ./= norm(b)
  λ = Inf

  for i = 1:maxiter
    copy!(bᵒˡᵈ, b)
    mul!(b, A, bᵒˡᵈ)

    λᵒˡᵈ = λ
    λ = dot(bᵒˡᵈ, b) / dot(bᵒˡᵈ, bᵒˡᵈ)
    b ./= norm(b)

    #verbose && println("iter = $i; λ = $λ")
    abs(λ / λᵒˡᵈ - 1) < rtol && break
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
    res[:, cnt, :] .= reshape(img[m.+ms, :], :, ncoils)
  end
  return res
end

"""
    crop(A::Array{T,D}, s::NTuple{R,Int64}) where {T,D,R}

crop the central area (of size s) of the first dimensions of the array A
"""
function crop(A::Array{T,D}, s::NTuple{R,Int64}) where {T,D,R}

  idx_center = CartesianIndices(s) .- CartesianIndex(s .÷ 2) .+ CartesianIndex(size(A)[1:R] .÷ 2)

  return A[idx_center, CartesianIndices(size(A)[(R+1):end])]
end

"""
  findIndices()
  Calculates the indices which place the center of k-space sampled on oldEnc in the correct place w.r.t a grid of size newEnc
"""
function findIndices(newEnc::NTuple{D,Int64}, oldEnc::NTuple{D,Int64}) where {D}
  shift = (newEnc .- oldEnc) .÷ 2

  oldCart = CartesianIndices(oldEnc) .+ CartesianIndex(shift)
  newCart = CartesianIndices(newEnc)

  if any(newEnc .> oldEnc)
    return LinearIndices(newCart)[oldCart][:]
  else
    return LinearIndices(oldCart)[:]
  end
end
