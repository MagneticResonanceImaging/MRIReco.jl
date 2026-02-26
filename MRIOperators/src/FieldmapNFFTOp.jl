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

mutable struct FieldmapNFFTOp{T, vecT <: AbstractVector{Complex{T}},F1,F2,D, vecI, matT, vecTR, S} <: AbstractMRIOperator{Complex{T}}
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
end

LinearOperators.storage_type(::FieldmapNFFTOp{T, vecT}) where {T, vecT} = vecT

"""
    FieldmapNFFTOp(shape::NTuple{D,Int64}, tr::Trajectory,
                        correctionmap::Array{ComplexF64,D};
                        method::String="nfft",
                        echoImage::Bool=true,
                        alpha::Float64=1.75,
                        m::Float64=3.0,
                        K=20,
                        kargs...) where D

generates a `FieldmapNFFTOp` which evaluates the MRI Fourier signal encoding operator,
including B0-inhomogeneities using time-segmented NFFTs.

# Arguments:
* `shape::NTuple{D,Int64}`             - size of image to encode/reconstruct
* `tr::Trajectory`                     - Trajectory with the kspace nodes to sample
* `correctionmap::Array{ComplexF64,D}` - fieldmap for the correction of off-resonance effects
* (`method::String="nfft"`)            - method to use for time-segmentation when correctio field inhomogeneities
* (`echoImage::Bool=false`)            - if true sampling times will only be considered relative to the echo time
                                         this results in complex valued image even for real-valued input.
* (`alpha::Float64=1.75`)              - oversampling factor for interpolation
* (`m::Float64=3.0`)                   - truncation size of interpolation kernel
* (`K=20`)                             - number of translates for LeastSquares approaches
                                         (not NFFT-approach) to time-segmentation
* (`kargs`)                            - additional keyword arguments
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
  nrow = size(nodes,2)
  ncol = prod(shape)

 # create and truncate low-rank expansion
  cparam = createInhomogeneityData_(vec(times), correctionmap; K=K, alpha=alpha, m=m, method=method, K_tol=K_tol, numSamp=numSamplingPerProfile(tr),step=step)
  K = size(cparam.A_k,2)

  @debug "K = $K"

  plans = [] # Dont fully specify type yet
  idx = []
  baseArrayType = stripParameters(S)
  for κ=1:K
    posIndices = findall(x->x!=0.0, cparam.A_k[:,κ])
    if !isempty(posIndices)
      push!(idx, posIndices)
      push!(plans, plan_nfft(baseArrayType, nodes[:,idx[κ]], shape, m=3, σ=1.25, precompute = AbstractNFFTs.POLYNOMIAL))
    end
  end
  plans = identity.(plans) # This gives the vectors (more) concrete types
  idx = identity.(idx)
  K=length(plans)
  
  cparam = adapt(baseArrayType, cparam)
  idx = adapt.(baseArrayType, idx)

  if !isa(first(idx), Array)
    idx = map(i -> Int32.(i), idx)
  end

  tmp = S(undef, 0)
  x_tmp = fill!(similar(tmp, Complex{T}, ncol), zero(Complex{T}))
  y_tmp = fill!(similar(tmp, Complex{T}, nrow), zero(Complex{T}))
  
  d = [fill!(similar(tmp, Complex{T}, length(idx[κ])), zero(Complex{T})) for κ=1:K ]
  p = [fill!(similar(tmp, Complex{T}, shape), zero(Complex{T})) for κ=1:K]

  circTraj = isCircular(tr)

  return FieldmapNFFTOp{T, typeof(tmp), Nothing, Function, D, eltype(idx), typeof(cparam.A_k), typeof(cparam.times), typeof(scheduler)}(nrow, ncol, false, false
            , (res,x) -> produ!(res,x,x_tmp,shape,plans,idx,cparam,circTraj,d,p,scheduler)
            , nothing
            , (res,y) -> ctprodu!(res,y,y_tmp,shape,plans,idx,cparam,circTraj,d,p,scheduler), 0, 0, 0, tmp, tmp 
            , plans, idx, circTraj, shape, cparam, scheduler)
end

function Base.copy(cparam::InhomogeneityData{T}) where T

  return cparam

end

function Base.copy(S::FieldmapNFFTOp{T}) where {T}

  shape = S.shape

  K=length(S.plans)
  plans = [copy(S.plans[i]) for i=1:K]
  idx = copy(S.idx)

  x_tmp = fill!(similar(S.Mv, Complex{T}, S.ncol), zero(Complex{T}))
  y_tmp = fill!(similar(S.Mv, Complex{T}, S.nrow), zero(Complex{T}))

  cparam = copy(S.cparam)
  d = [fill!(similar(S.Mv, Complex{T}, length(idx[κ])), zero(Complex{T})) for κ=1:K ]
  p = [fill!(similar(S.Mv, Complex{T}, shape), zero(Complex{T})) for κ=1:K]

  D_ = length(shape)
  circTraj = S.circTraj

  scheduler = S.scheduler

  mul! = (res,x) -> produ!(res,x,x_tmp,shape,plans,idx,cparam,circTraj,d,p, scheduler)
  ctmul! = (res,y) -> ctprodu!(res,y,y_tmp,shape,plans,idx,cparam,circTraj,d,p, scheduler)

  return FieldmapNFFTOp{T, typeof(x_tmp), Nothing, Function, D_, eltype(idx), typeof(cparam.A_k), typeof(cparam.times), typeof(scheduler)}(S.nrow, S.ncol, false, false
            , mul!
            , nothing
            , ctmul!, 0, 0, 0, similar(x_tmp, 0), similar(x_tmp, 0)
            , plans, idx, circTraj, shape, cparam, scheduler)
end

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
