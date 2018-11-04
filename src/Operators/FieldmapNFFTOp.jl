export FieldmapNFFTOp, createInhomogeneityData_

mutable struct FieldmapNFFTOp{T,F1<:FuncOrNothing,F2<:FuncOrNothing,F3<:FuncOrNothing} <:
                      AbstractLinearOperator{T,F1,F2,F3}
  nrow :: Int
  ncol :: Int
  symmetric :: Bool
  hermitian :: Bool
  prod :: Function
  tprod :: F1
  ctprod :: F2
  inv :: F3
  density::Vector{Float64}
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

#
# Linear Operator to perform NFFT
#
function FieldmapNFFTOp(shape::Tuple, tr::AbstractTrajectory, correctionmap::Matrix
                        ; method="nfft"
                        , symmetrize=true
                        , echoImage=true
                        , alpha::Float64=1.75
                        , m::Float64=4.0
                        , K=20)

  nodes,times = kspaceNodes(tr), readoutTimes(tr)
  if echoImage
    times = times .- echoTime(tr)
  end
  nrow = size(nodes,2)
  ncol = prod(shape)

 # create and truncate low-rank expansion
  cparam = createInhomogeneityData_(nrow,ncol,vec(times),correctionmap; K=K, alpha=alpha, m=m, method=method)
  K = size(cparam.A_k,2)

  @info "K = $K"

  plan = Vector{NFFTPlan}(undef,K)
  idx = Vector{Vector{Int64}}(undef,K)
  for κ=1:K
    idx[κ] = findall(x->x!=0.0, cparam.A_k[:,κ])
    plan[κ] = NFFTPlan(nodes[:,idx[κ]], shape, 3, 1.25)
  end

  planTmp = NFFTPlan(nodes, shape, 3, 1.25)
  density = convert(Vector{Float64}, sdc(planTmp))

  mul(x::Vector{T}) where T<:ComplexF64 =
     produ(x,nrow,ncol,shape,plan,idx,cparam,density,symmetrize)
  ctmul(y::Vector{T}) where T<:ComplexF64 =
     ctprodu(y,shape,plan,idx,cparam,density,symmetrize,isCircular(tr))
  inverse(y::Vector{T}) where T<:ComplexF64 =
     inv(y,shape,plan,idx,cparam,density,symmetrize,isCircular(tr))

  return FieldmapNFFTOp{ComplexF64,Nothing,Function,Function}(nrow, ncol, false, false
            , mul
            , nothing
            , ctmul
            , inverse
            , density )
end

# function produ{T<:ComplexF64}(x::Vector{T}, numOfNodes::Int, numOfPixel::Int, shape::Tuple, plan::Vector{NFFTPlan{2,0,ComplexF64}}, cparam::InhomogeneityData, density::Vector{Float64}, symmetrize::Bool)
function produ(x::Vector{T}, numOfNodes::Int, numOfPixel::Int, shape::Tuple, plan,
               idx::Vector{Vector{Int64}}, cparam::InhomogeneityData,
                density, symmetrize::Bool) where T<:ComplexF64
  K = size(cparam.A_k,2)
  s = zeros(ComplexF64,numOfNodes)
  # Preprocessing step when time and correctionMap are centered
  if cparam.method == "nfft"
      x = x .* exp.(-vec(cparam.Cmap) * cparam.t_hat )
  end

  sp = Threads.SpinLock()
  @time Threads.@threads for κ=1:K
    p_tild = cparam.C_k[κ,:] .* x;
    s_tild = zeros(ComplexF64, numOfNodes);
    s_tild[idx[κ]] = reshape(nfft(plan[κ], reshape(p_tild, shape)), length(idx[κ]) );
    s_tild = cparam.A_k[:,κ] .* s_tild
    lock(sp)
    s[:] += s_tild
    unlock(sp)
  end

  # Postprocessing step when time and correctionMap are centered
  if cparam.method == "nfft"
      s[:] .*= exp.(-cparam.z_hat*(cparam.times .- cparam.t_hat) )
  end
  if symmetrize
      s[:] .*= sqrt.(density) # <- use for FISTA
  end
  return s
end

# function inv{T<:ComplexF64}(x::Vector{T}, shape::Tuple, plan::Vector{NFFTPlan{2,0,ComplexF64}}, cparam::InhomogeneityData, density::Vector{Float64}, symmetrize::Bool)
function inv(x::Vector{T}, shape::Tuple, plan, idx::Vector{Vector{Int64}},
             cparam::InhomogeneityData, density, symmetrize::Bool,
            shutter::Bool) where T<:ComplexF64
  if symmetrize
    x = x .* sqrt.(density)
  else
    x = x .* density
  end

  y = ctprodu(x,shape,plan,cparam,density,false,shutter)
end

# function ctprodu{T<:ComplexF64}(x::Vector{T}, shape::Tuple, plan::Vector{NFFTPlan{2,0,ComplexF64}}, cparam::InhomogeneityData, density::Vector{Float64}, symmetrize::Bool)
function ctprodu(x::Vector{T}, shape::Tuple, plan, idx::Vector{Vector{Int64}},
                 cparam::InhomogeneityData, density, symmetrize::Bool,
                 shutter::Bool) where T<:ComplexF64

  y = zeros(ComplexF64,prod(shape))
  K = size(cparam.A_k,2)

  if symmetrize
      x = x .* sqrt.(density) # <- using for FISTA
  end

  # Preprocessing step when time and correctionMap are centered
  if cparam.method == "nfft"
      x = x .* conj.(exp.(-cparam.z_hat*(cparam.times .- cparam.t_hat)))
  end
  # Algorithm for fast adjoint Transformation with correctionterm
  #@time y = @distributed (+)
  sp = Threads.SpinLock()
  @time Threads.@threads for κ=1:K
      s_tild = conj.(cparam.A_k[:,κ]) .* x
      p_tild = reshape(nfft_adjoint(plan[κ], s_tild[idx[κ]]), prod(shape) )
      p_tild = conj.(cparam.C_k[κ,:]) .* p_tild
      lock(sp)
      y[:] += p_tild
      unlock(sp)
  end

  if cparam.method == "nfft"
    y[:] .*=  conj(exp.(-vec(cparam.Cmap) * cparam.t_hat))
  end

  if shutter
    circularShutter!(reshape(y, shape), 1.0)
  end

  return y
end

####################### Helper Function ########################################
function createInhomogeneityData_( numOfNodes::Int64
                                , numOfPixel::Int64
                                , times::Vector
                                , correctionmap::Matrix
                                ; K::Int64=20
                                , alpha::Float64=1.75
                                , m::Float64 = 4.0
                                , method="nfft")

    C = getC_Coefficients_hist_lsqr(K,numOfNodes,numOfPixel,vec(times),vec(correctionmap))
    if method == "const"
      A = getA_Coefficients_one_term(K,numOfNodes,numOfPixel,vec(times),vec(correctionmap))
    elseif method== "linear"
      A = getA_Coefficients_two_terms(K,numOfNodes,numOfPixel,vec(times),vec(correctionmap))
    elseif method == "nfft"
      A,C = getA_Ccoefficients_nfft(K,numOfNodes,numOfPixel,vec(times),vec(correctionmap), alpha, m)
    elseif method == "leastsquare"
        A,C = getA_Coefficients_least_Squares(K,numOfNodes,numOfPixel,vec(times),vec(correctionmap))
    elseif method == "hist"
      A = getA_Coefficients_hist_lsqr(K,numOfNodes,numOfPixel,vec(times),vec(correctionmap))
      C = getC_Coefficients_hist_lsqr(K,numOfNodes,numOfPixel,vec(times),vec(correctionmap))
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
