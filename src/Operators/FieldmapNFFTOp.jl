export FieldmapNFFTOp, createInhomogeneityData_

mutable struct FieldmapNFFTOp{T,F1<:FuncOrNothing,F2<:FuncOrNothing} <:
                      AbstractLinearOperator{T,Function,F1,F2}
  nrow :: Int
  ncol :: Int
  symmetric :: Bool
  hermitian :: Bool
  prod :: Function
  tprod :: F1
  ctprod :: F2
end

mutable struct InhomogeneityData
  A_k::Matrix{ComplexF64}
  C_k::Matrix{ComplexF64}
  times::Vector{Float64}
  Cmap::Matrix{ComplexF64}
  t_hat::Float64
  z_hat::ComplexF64
  method::String
end

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
                        correctionmap::Array{ComplexF64,D};
                        method::String="nfft",
                        echoImage::Bool=true,
                        alpha::Float64=1.75,
                        m::Float64=3.0,
                        K=20,
                        kargs...) where D

  nodes,times = kspaceNodes(tr), readoutTimes(tr)
  if echoImage
    times = times .- echoTime(tr)
  end
  nrow = size(nodes,2)
  ncol = prod(shape)

 # create and truncate low-rank expansion
  cparam = createInhomogeneityData_(vec(times), correctionmap; K=K, alpha=alpha, m=m, method=method)
  K = size(cparam.A_k,2)

  @debug "K = $K"

  plan = Vector{NFFTPlan{D,0,Float64}}(undef,K)
  idx = Vector{Vector{Int64}}(undef,K)
  for κ=1:K
    idx[κ] = findall(x->x!=0.0, cparam.A_k[:,κ])
    plan[κ] = NFFTPlan(nodes[:,idx[κ]], shape, 3, 1.25, precompute = NFFT.FULL)
  end

  p = [zeros(ComplexF64, ncol) for t=1:Threads.nthreads() ]
  y = [zeros(ComplexF64, nrow) for t=1:Threads.nthreads() ]
  d = [zeros(ComplexF64, length(idx[κ])) for κ=1:K ]

  mul(x::Vector{T}) where T<:ComplexF64 =
     produ(x,nrow,ncol,shape,plan,idx,cparam,isCircular(tr),p,y,d)
  ctmul(y::Vector{T}) where T<:ComplexF64 =
     ctprodu(y,shape,plan,idx,cparam,isCircular(tr),p,y,d)
  inverse(y::Vector{T}) where T<:ComplexF64 =
     inv(y,shape,plan,idx,cparam,isCircular(tr),p,y,d)

  return FieldmapNFFTOp{ComplexF64,Nothing,Function}(nrow, ncol, false, false
            , mul
            , nothing
            , ctmul)
end

# function produ{T<:ComplexF64}(x::Vector{T}, numOfNodes::Int, numOfPixel::Int, shape::Tuple, plan::Vector{NFFTPlan{2,0,ComplexF64}}, cparam::InhomogeneityData, density::Vector{Float64}, symmetrize::Bool)
function produ(x::Vector{T}, numOfNodes::Int, numOfPixel::Int, shape::Tuple, plan,
               idx::Vector{Vector{Int64}}, cparam::InhomogeneityData,
               shutter::Bool, p, y, d) where T<:ComplexF64
  K = size(cparam.A_k,2)
  s = zeros(ComplexF64,numOfNodes)

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
  produ_inner(K,cparam.C_k, cparam.A_k, shape, p, d, y, s, sp, plan, idx, x_)

  # Postprocessing step when time and correctionMap are centered
  if cparam.method == "nfft"
      s .*= exp.(-cparam.z_hat*(cparam.times .- cparam.t_hat) )
  end
  return s
end

function produ_inner(K, C, A, shape, p, d, y, s, sp, plan, idx, x_)
  Threads.@threads for κ=1:K
    t = Threads.threadid()
    for l=1:length(x_)
      p[t][l] = C[κ,l] * x_[l]
    end
    NFFT.nfft!(plan[κ], reshape(p[t], shape), d[κ])
    fill!(y[t], 0.0)
    for k=1:length(idx[κ])
      y[t][idx[κ][k]] = d[κ][k]
    end
    for l=1:size(A,1)
      y[t][l] *= A[l,κ]
    end
    lock(sp)
    s .+= y[t]
    unlock(sp)
  end
  return
end


# function ctprodu{T<:ComplexF64}(x::Vector{T}, shape::Tuple, plan::Vector{NFFTPlan{2,0,ComplexF64}}, cparam::InhomogeneityData, density::Vector{Float64}, symmetrize::Bool)
function ctprodu(x::Vector{T}, shape::Tuple, plan, idx::Vector{Vector{Int64}},
                 cparam::InhomogeneityData, shutter::Bool, p, y_, d) where T<:ComplexF64

  y = zeros(ComplexF64,prod(shape))
  K = size(cparam.A_k,2)

  x_ = copy(x)

  # Preprocessing step when time and correctionMap are centered
  if cparam.method == "nfft"
      x_ .*= conj.(exp.(-cparam.z_hat*(cparam.times .- cparam.t_hat)))
  end

  sp = Threads.SpinLock()
  ctprodu_inner(K,cparam.C_k, cparam.A_k, shape, p, d, y, sp, plan, idx, x_)

  if cparam.method == "nfft"
    y .*=  conj(exp.(-vec(cparam.Cmap) * cparam.t_hat))
  end

  if shutter
    circularShutter!(reshape(y, shape), 1.0)
  end

  return y
end

function ctprodu_inner(K, C, A, shape, p, d, y, sp, plan, idx, x_)
  Threads.@threads for κ=1:K
    t = Threads.threadid()
    for k=1:length(idx[κ])
      d[κ][k] = conj.(A[idx[κ][k],κ]) * x_[idx[κ][k]]
    end
    NFFT.nfft_adjoint!(plan[κ], d[κ], reshape(p[t], shape))
    for k=1:length(p[t])
      p[t][k] = conj(C[κ,k]) * p[t][k]
    end
    lock(sp)
    y .+= p[t]
    unlock(sp)
  end
  return
end


function createInhomogeneityData_(times::Vector,
                                  correctionmap::Array{ComplexF64,D};
                                  K::Int64=20,
                                  alpha::Float64=1.75,
                                  m::Float64 = 4.0,
                                  method="nfft") where D

    if method == "const"
      A = get_AC_coefficients_one_term(K,vec(times),vec(correctionmap))
      C = get_C_coefficients_hist_lsqr(K,vec(times),vec(correctionmap))
    elseif method== "linear"
      A = get_AC_coefficients_two_terms(K,vec(times),vec(correctionmap))
      C = get_C_coefficients_hist_lsqr(K,vec(times),vec(correctionmap))
    elseif method == "nfft"
      A,C = get_AC_coefficients_nfft(K,vec(times),vec(correctionmap), alpha, m)
    elseif method == "leastsquare"
      A,C = get_AC_coefficients_lsqr(K,vec(times),vec(correctionmap))
    elseif method == "hist"
      A = get_A_coefficients_hist_lsqr(K,vec(times),vec(correctionmap))
      C = get_C_coefficients_hist_lsqr(K,vec(times),vec(correctionmap))
    else
      error("approximation scheme $(interp) is not yet implemented")
    end
    if size(A,2) != size(C,1)
      error("Consistency check failed! A and C are not compatible")
    end

    t_hat = (times[1] + times[end])/2
    z_hat = minimum(real(correctionmap)) + maximum(real(correctionmap))
    z_hat += 1im*(minimum(imag(correctionmap)) + maximum(imag(correctionmap)))
    z_hat *= 0.5

    return InhomogeneityData(A,C,vec(times),correctionmap,t_hat,z_hat,method)
end

function adjoint(op::FieldmapNFFTOp{T}) where T
  return LinearOperator{T,Function,Nothing,Function}(op.ncol, op.nrow, op.symmetric, op.hermitian,
                        op.ctprod, nothing, op.prod)
end
