export FieldmapNFFTOp, InhomogeneityData, createInhomogeneityData_

mutable struct InhomogeneityData{T}
  A_k::Matrix{Complex{T}}
  C_k::Matrix{Complex{T}}
  times::Vector{T}
  Cmap::Matrix{Complex{T}}
  t_hat::T
  z_hat::Complex{T}
  method::String
end

mutable struct FieldmapNFFTOp{T,F1,F2,D} <:AbstractLinearOperator{Complex{T}}
  nrow :: Int
  ncol :: Int
  symmetric :: Bool
  hermitian :: Bool
  prod! :: Function
  tprod! :: F1
  ctprod! :: F2
  nprod :: Int
  ntprod :: Int
  nctprod :: Int
  args5 :: Bool
  use_prod5! :: Bool
  allocated5 :: Bool
  Mv5 :: Vector{Complex{T}}
  Mtu5 :: Vector{Complex{T}}
  plans
  idx::Vector{Vector{Int64}}
  circTraj::Bool
  shape::NTuple{D,Int64}
  cparam::InhomogeneityData{T}
end

# LinearOperators.storage_type(op::FieldmapNFFTOp) = typeof(op.Mv5)

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

  plans = Vector{NFFT.NFFTPlan{T,D,1}}(undef,K)
  idx = Vector{Vector{Int64}}(undef,K)
  for κ=1:K
    idx[κ] = findall(x->x!=0.0, cparam.A_k[:,κ])
    plans[κ] = plan_nfft(nodes[:,idx[κ]], shape, m=3, σ=1.25, precompute = NFFT.FULL)
  end

  d = [zeros(ComplexF64, length(idx[κ])) for κ=1:K ]
  
  circTraj = isCircular(tr)

  mul!(res, x::Vector{T}) where T =
     (res .= produ(x,nrow,ncol,shape,plans,idx,cparam,circTraj,d))
  ctmul!(res, y::Vector{T}) where T =
     (res .= ctprodu(y,shape,plans,idx,cparam,circTraj,d))

  return FieldmapNFFTOp{T,Nothing,Function,D}(nrow, ncol, false, false
            , mul!
            , nothing
            , ctmul!, 0, 0, 0, false, false, false, ComplexF64[], ComplexF64[]
            , plans, idx, circTraj, shape, cparam)
end

function Base.copy(S::FieldmapNFFTOp{T,Nothing,Function,D}) where {T,D}
  K=length(S.plans)
  plans = [copy(S.plans[i]) for i=1:K]
  idx = deepcopy(S.idx)

  d = [zeros(Complex{T}, length(idx[κ])) for κ=1:K ]

  cparam = deepcopy(S.cparam)

  mul!(res, x::Vector{T}) where T =
     (res .= produ(x,S.nrow,S.ncol,S.shape,plans,idx,cparam,S.circTraj,d))
  ctmul!(res, y::Vector{T}) where T =
     (res .= ctprodu(y,S.shape,plans,idx,cparam,S.circTraj,d))

  D_ = length(S.shape)
  return FieldmapNFFTOp{T,Nothing,Function,D_}(S.nrow, S.ncol, false, false
            , mul!
            , nothing
            , ctmul!, 0, 0, 0, false, false, false, Complex{T}[], Complex{T}[]
            , plans, idx, S.circTraj, S.shape, cparam)
end

# function produ{T<:Float64}(x::Vector{T}, numOfNodes::Int, numOfPixel::Int, shape::Tuple, plan::Vector{NFFTPlan{Float64,2,1}}, cparam::InhomogeneityData, density::Vector{Float64}, symmetrize::Bool)
function produ(x::Vector{T}, numOfNodes::Int, numOfPixel::Int, shape::Tuple, plan,
               idx::Vector{Vector{Int64}}, cparam::InhomogeneityData,
               shutter::Bool, d) where T
  K = size(cparam.A_k,2)
  s = zeros(T,numOfNodes)

  if shutter
    circularShutter!(reshape(x, shape), 1.0)
  end

  # Preprocessing step when time and correctionMap are centered
  if cparam.method == "nfft"
    x_ = x .* exp.(-vec(cparam.Cmap) * cparam.t_hat )
  else
    x_ = copy(x)
  end

  sp = Threads.SpinLock()
  produ_inner(K,cparam.C_k, cparam.A_k, shape, d, s, sp, plan, idx, x_)

  # Postprocessing step when time and correctionMap are centered
  if cparam.method == "nfft"
      s .*= exp.(-cparam.z_hat*(cparam.times .- cparam.t_hat) )
  end
  return s
end

function produ_inner(K, C, A, shape, d, s, sp, plan, idx, x_)
  @sync for κ=1:K
    Threads.@spawn begin
      p_ = C[κ,:] .* x_
      mul!(d[κ], plan[κ], reshape(p_, shape))
    
      lock(sp)
      for k=1:length(idx[κ])
        l=idx[κ][k]
        s[l] += d[κ][k]*A[l,κ]
      end
      unlock(sp)
    end
  end
  
  return
end


# function ctprodu{T<:ComplexF64}(x::Vector{T}, shape::Tuple, plan::Vector{NFFTPlan{Float64,2,1}, cparam::InhomogeneityData, density::Vector{Float64}, symmetrize::Bool)
function ctprodu(x::Vector{Complex{T}}, shape::Tuple, plan, idx::Vector{Vector{Int64}},
                 cparam::InhomogeneityData{T}, shutter::Bool, d) where T

  y = zeros(Complex{T},prod(shape))
  K = size(cparam.A_k,2)

  x_ = copy(x)

  # Preprocessing step when time and correctionMap are centered
  if cparam.method == "nfft"
      x_ .*= conj.(exp.(-cparam.z_hat*(cparam.times .- cparam.t_hat)))
  end

  sp = Threads.SpinLock()
  ctprodu_inner(K,cparam.C_k, cparam.A_k, shape, d, y, sp, plan, idx, x_)

  if cparam.method == "nfft"
    y .*=  conj(exp.(-vec(cparam.Cmap) * cparam.t_hat))
  end

  if shutter
    circularShutter!(reshape(y, shape), 1.0)
  end

  return y
end

function ctprodu_inner(K, C, A, shape, d, y, sp, plan, idx, x_)

  @sync for κ=1:K
    Threads.@spawn begin
       for k=1:length(idx[κ])
         d[κ][k] = conj.(A[idx[κ][k],κ]) * x_[idx[κ][k]]
       end
       p_ = zeros(ComplexF64, shape) 
     
       mul!(p_, adjoint(plan[κ]), d[κ])
     
       lock(sp)
       for k=1:length(p_)
         y[k] += conj(C[κ,k]) * p_[k]
       end
       unlock(sp)
    end   
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
