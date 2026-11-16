export FieldmapNFFTOp, InhomogeneityData, createInhomogeneityData_

include("ExpApproximation.jl")

mutable struct InhomogeneityData{T, matT <: AbstractArray{Complex{T}, 2}, vecT <: AbstractArray{T, 1}}
  A_k::matT
  C_k::matT
  times::vecT
  Cmap::matT
  t_hat::T
  z_hat::Complex{T}
  method::String
end

Adapt.adapt_structure(::Type{arrT}, data::InhomogeneityData) where arrT = InhomogeneityData(adapt(arrT, data.A_k), adapt(arrT, data.C_k), adapt(arrT, data.times), adapt(arrT, data.Cmap), data.t_hat, data.z_hat, data.method)

mutable struct FieldmapNFFTOp{T, vecT <: AbstractVector{Complex{T}},F1,F2,D, vecI, matT, vecTR, S, M, B} <: AbstractMRIOperator{Complex{T}}
  const nrow :: Int
  const ncol :: Int
  const symmetric :: Bool
  const hermitian :: Bool
  const prod! :: Function
  const tprod! :: F1
  const ctprod! :: F2
  nprod :: Int
  ntprod :: Int
  nctprod :: Int
  Mv :: vecT
  Mtu :: vecT
  const plans
  const idx::Vector{vecI}
  circTraj::Bool
  const shape::NTuple{D,Int64}
  cparam::InhomogeneityData{T, matT, vecTR}
  const scheduler::S
  # inputs of the NFFT plans, kept to rebuild them in copy and Adapt
  const nodes::Vector{M}
  const backend::B
end

LinearOperators.storage_type(::FieldmapNFFTOp{T, vecT}) where {T, vecT} = vecT

"""
    FieldmapNFFTOp(shape::NTuple{D,Int64}, tr::Trajectory,
                        correctionmap::Array{Complex{T},D};
                        method::String="nfft",
                        echoImage::Bool=true,
                        alpha::Float64=1.25,
                        m::Float64=2.0,
                        K=20,
                        K_tol::Float64=1.e-3,
                        step::Int64=10,
                        S = Vector{Complex{T}},
                        scheduler = DynamicScheduler(),
                        kargs...) where {T,D}

generates a `FieldmapNFFTOp` which evaluates the MRI Fourier signal encoding operator,
including off-resonance and relaxation, as a weighted sum of NFFTs.

# Arguments:
* `shape::NTuple{D,Int64}`             - size of image to encode/reconstruct
* `tr::Trajectory`                     - Trajectory with the kspace nodes to sample
* `correctionmap::Array{Complex{T},D}` - complex correction map `z` with the relaxation rate R₂* in 1/s as its real
                                         part and the off-resonance frequency in rad/s as its imaginary part
* (`method::String="nfft"`)            - approximation of `exp(-z t)`, where `t` is the sampling time: `"nfft"`,
                                         `"leastsquare"`, `"hist"` or `"svd"`; for current failures, see
                                         [#360](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/360) and
                                         [#363](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/363)
* (`echoImage::Bool=true`)             - if true, the sampling times are counted from the echo time `TE` of the
                                         trajectory, so the image is the magnetization at that time, which with
                                         off-resonance is complex-valued even for a real-valued object
* (`alpha::Float64=1.25`)              - oversampling factor of the `"nfft"` approximation
* (`m::Float64=2.0`)                   - kernel size of the `"nfft"` approximation
* (`K=20`)                             - number of terms of `"leastsquare"` and `"hist"` and maximum number of terms
                                         of `"svd"`; `"nfft"` determines the number itself
* (`K_tol::Float64=1.e-3`)             - tolerance with which `"svd"` chooses the number of terms
* (`step::Int64=10`)                   - `"svd"` fits its coefficients to every `step`-th sample of a readout
* (`S = Vector{Complex{T}}`)           - storage vector type, for example `CuVector{ComplexF32}` on a GPU
* (`scheduler = DynamicScheduler()`)   - scheduler of the parallel loop over the terms
* (`kargs`)                            - further keyword arguments, which are ignored, so NFFT options such as
                                         `oversamplingFactor` have no effect

# Moving to another array backend

    Adapt.adapt(CuArray, op::FieldmapNFFTOp)

returns a new `FieldmapNFFTOp` with the approximation coefficients of `op` moved by `Adapt.adapt(CuArray, ...)`,
so they are not computed again. Its NFFT plans and work arrays are new: the plans are built from k-space nodes kept
on the CPU, with the NFFT backend that was active when `op` was constructed. `op` itself is not modified.
The backend may choose other plan defaults on the new device (NonuniformFFTs uses another window on CUDA),
so the results can differ from those of `op` within the NFFT accuracy.
With JLArrays the operator was tested with the NFFT.jl backend only; NonuniformFFTs does not create plans for them.
`SensitivityOp`, `SubspaceOp`, `MapSliceOp`, `ExplicitOp` and `FieldmapNFFTOp` are the MRI operators that support `Adapt.adapt` so far.
"""
function FieldmapNFFTOp(shape::NTuple{D,Int64}, tr::Trajectory,
                        correctionmap::Array{Complex{T},D};
                        method::String="nfft",
                        echoImage::Bool=true,
                        alpha::Float64=1.25,
                        m::Float64=2.0,
                        K=20,
                        K_tol::Float64=1.e-3,
                        step::Int64=10,
                        S = Vector{Complex{T}},
                        scheduler = DynamicScheduler(),
                        kargs...) where {T,D}
  nodes,times = kspaceNodes(tr), readoutTimes(tr)
  if echoImage
    times = times .- echoTime(tr)
  end

 # create and truncate low-rank expansion
  cparam = createInhomogeneityData_(vec(times), correctionmap; K=K, alpha=alpha, m=m, method=method, K_tol=K_tol, numSamp=numSamplingPerProfile(tr),step=step)
  K = size(cparam.A_k,2)

  @debug "K = $K"

  idx = Vector{Int}[]
  segmentNodes = Matrix{eltype(nodes)}[]
  for κ=1:K
    posIndices = findall(x->x!=0.0, cparam.A_k[:,κ])
    if !isempty(posIndices)
      push!(idx, posIndices)
      push!(segmentNodes, nodes[:,idx[κ]])
    end
  end

  baseArrayType = stripParameters(S)
  return FieldmapNFFTOp(S, shape, segmentNodes, adapt.(baseArrayType, idx), adapt(baseArrayType, cparam), isCircular(tr), scheduler, AbstractNFFTs.active_backend())
end

# Builds the operator from the stored nodes and time-segmentation data; used by the constructor above, copy and Adapt.
function FieldmapNFFTOp(S, shape::NTuple{D,Int64}, nodes, idx, cparam::InhomogeneityData{T}, circTraj, scheduler, backend) where {T,D}
  plans = [plan_nfft(backend, stripParameters(S), k, shape, m=3, σ=1.25, precompute = AbstractNFFTs.POLYNOMIAL) for k in nodes]
  idx = map(i -> isa(i, Array) ? convert(Vector{Int}, i) : Int32.(i), idx)
  nrow, ncol = size(cparam.A_k, 1), prod(shape)
  K = length(plans)

  tmp = S(undef, 0)
  x_tmp = fill!(similar(tmp, Complex{T}, ncol), zero(Complex{T}))
  y_tmp = fill!(similar(tmp, Complex{T}, nrow), zero(Complex{T}))
  
  d = [fill!(similar(tmp, Complex{T}, length(idx[κ])), zero(Complex{T})) for κ=1:K ]
  p = [fill!(similar(tmp, Complex{T}, shape), zero(Complex{T})) for κ=1:K]

  return FieldmapNFFTOp{T, typeof(tmp), Nothing, Function, D, eltype(idx), typeof(cparam.A_k), typeof(cparam.times), typeof(scheduler), eltype(nodes), typeof(backend)}(nrow, ncol, false, false
            , (res,x) -> produ!(res,x,x_tmp,shape,plans,idx,cparam,circTraj,d,p,scheduler)
            , nothing
            , (res,y) -> ctprodu!(res,y,y_tmp,shape,plans,idx,cparam,circTraj,d,p,scheduler), 0, 0, 0, tmp, tmp 
            , plans, idx, circTraj, shape, cparam, scheduler, nodes, backend)
end

Adapt.adapt_structure(to, op::FieldmapNFFTOp) =
  FieldmapNFFTOp(typeof(adapt(to, LinearOperators.storage_type(op)(undef, 0))), op.shape, op.nodes, map(i -> adapt(to, i), op.idx),
                 adapt(to, op.cparam), op.circTraj, op.scheduler, op.backend)

function Base.copy(cparam::InhomogeneityData{T}) where T

  return cparam

end

Base.copy(S::FieldmapNFFTOp) = FieldmapNFFTOp(LinearOperators.storage_type(S), S.shape, S.nodes, S.idx, copy(S.cparam),
                                              S.circTraj, S.scheduler, S.backend)

function produ!(s::AbstractVector{T}, x::AbstractVector{T}, x_tmp::AbstractVector{T},shape::Tuple, plan,
               idx, cparam::InhomogeneityData,
               shutter::Bool, d, p, scheduler::Scheduler = DynamicScheduler()) where {T}

  s .= zero(T)
  K=length(plan)

  if shutter
    circularShutter!(reshape(x, shape), 1.0)
  end

  # Preprocessing step when time and correctionMap are centered
  x_tmp .= x
  if cparam.method == "nfft"
    x_tmp .*= exp.(-vec(cparam.Cmap) * cparam.t_hat )
  end

  sp = Threads.SpinLock()
  produ_inner!(K,cparam.C_k, cparam.A_k, shape, d, s, sp, plan, idx, x_tmp, p, scheduler)

  # Postprocessing step when time and correctionMap are centered
  if cparam.method == "nfft"
      s .*= exp.(-cparam.z_hat*(cparam.times .- cparam.t_hat) )
  end
end

function produ_inner!(K, C, A, shape, d, s, sp, plan, idx, x_, p, scheduler::Scheduler = DynamicScheduler())
  @tasks for κ=1:K
    @set scheduler = scheduler
    p[κ][:] .= C[κ,:] .* x_
    mul!(d[κ], plan[κ], p[κ])
    
    lock(sp)
    for k=1:length(idx[κ])
      l=idx[κ][k]
      s[l] += d[κ][k]*A[l,κ]
    end
    unlock(sp)
  end
  
  return
end


function ctprodu!(y::AbstractVector{Complex{T}}, x::AbstractVector{Complex{T}}, x_tmp::AbstractVector{Complex{T}}, shape::Tuple, plan, idx,
                 cparam::InhomogeneityData{T}, shutter::Bool, d, p, scheduler::Scheduler = DynamicScheduler()) where {T}

  y .= zero(Complex{T})
  K=length(plan)

  x_tmp .= x

  # Preprocessing step when time and correctionMap are centered
  if cparam.method == "nfft"
      x_tmp .*= conj.(exp.(-cparam.z_hat*(cparam.times .- cparam.t_hat)))
  end

  sp = Threads.SpinLock()
  ctprodu_inner!(K,cparam.C_k, cparam.A_k, shape, d, y, sp, plan, idx, x_tmp, p, scheduler)

  if cparam.method == "nfft"
    y .*=  conj(exp.(-vec(cparam.Cmap) * cparam.t_hat))
  end

  if shutter
    circularShutter!(reshape(y, shape), 1.0)
  end

end

function ctprodu_inner!(K, C, A, shape, d, y, sp, plan, idx, x_, p, scheduler::Scheduler = DynamicScheduler())

  @tasks for κ=1:K
    @set scheduler = scheduler
     for k=1:length(idx[κ])
       d[κ][k] = conj.(A[idx[κ][k],κ]) * x_[idx[κ][k]]
     end
     
     mul!(p[κ], adjoint(plan[κ]), d[κ])
     
     lock(sp)
     for k=1:length(p[κ])
       y[k] += conj(C[κ,k]) * p[κ][k]
     end
     unlock(sp)
  end
    
  return
end


function createInhomogeneityData_(times::Vector{T},
                                  correctionmap::Array{Complex{T},D};
                                  K::Int64=20,
                                  K_tol::Float64=1.e-3,
                                  alpha::Float64=1.75,
                                  m::Float64 = 4.0,
                                  method="nfft",
                                  numSamp::Int64=length(times),
                                  step::Int64=10) where {T,D}

    if method == "nfft"
      A,C = get_AC_coefficients_nfft(K,vec(times),vec(correctionmap), alpha, m)
    elseif method == "leastsquare"
      A,C = get_AC_coefficients_lsqr(K,vec(times),vec(correctionmap))
    elseif method == "hist"
      A = get_A_coefficients_hist_lsqr(K,vec(times),vec(correctionmap))
      C = get_C_coefficients_hist_lsqr(K,vec(times),vec(correctionmap))
    elseif method == "leastsquareFS"
      A,C = get_AC_coefficients_lsqr_fs(K,vec(times),vec(correctionmap))
    elseif method == "histFS"
      A = get_A_coefficients_hist_lsqr_fs(K,vec(times),vec(correctionmap))
      # C = get_C_coefficients_hist_lsqr_fs(K,vec(times),vec(correctionmap))
      C = get_C_coefficients_lsqr(vec(times),vec(correctionmap),A,numSamp)
    elseif method == "svd"
      A = get_A_coefficients_svd(K,vec(times)[1:numSamp],vec(correctionmap); K_tol=K_tol)
      C = get_C_coefficients_lsqr(vec(times),vec(correctionmap),A,numSamp,step)
      rep = div(length(times),numSamp)
      A = repeat(A,outer=(rep,1))
    else
      error("approximation scheme $(interp) is not yet implemented")
    end
    if size(A,2) != size(C,1)
      @info "method: $(method), A: $(size(A)), C: $(size(C))"
      error("Consistency check failed! A and C are not compatible")
    end

    t_hat = (times[1] + times[end])/2
    z_hat = minimum(real(correctionmap)) + maximum(real(correctionmap))
    z_hat += 1im*(minimum(imag(correctionmap)) + maximum(imag(correctionmap)))
    z_hat *= T(0.5)

    return InhomogeneityData(A,C,vec(times),correctionmap,t_hat,z_hat,method)
end
