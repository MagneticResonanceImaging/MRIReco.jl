export FieldmapToeplitzOp, getMatrixElement

#
# Linear Operator to perform NFFT
#
# function FieldmapToeplitzOp(shape::Tuple, tr::AbstractTrajectory, correctionmap::Matrix
#                         ; symmetrize=true
#                         , K::Int = 30)
function FieldmapToeplitzOp(shape::Tuple, tr::AbstractTrajectory, correctionmap::Matrix
                        ; symmetrize=true, m::Float64 = 4.0, alpha::Float64=1.75)

  nodes,times = kspaceNodes(tr), readoutTimes(tr)
  nrow = prod(shape)
  ncol = prod(shape)

  plan = NFFTPlan(nodes, shape, 4, 1.75)
  density = sdc(plan)

  # cparam = lowRankParams(nrow,ncol,vec(times),correctionmap,K)
  cparam = lowRankParams(nrow,ncol,vec(times),correctionmap,m,alpha)
  K = size(cparam.A_k,2)

  ft = plan_fft( zeros(ComplexF64, Tuple(2*collect(shape)));flags=FFTW.MEASURE )
  ift = plan_ifft(zeros(ComplexF64, Tuple(2*collect(shape)));flags=FFTW.MEASURE)
  λ = zeros(ComplexF64,4*prod(shape),K)

  for κ=1:K
    λ[:,κ] = toeplitzEigvals(κ, shape, tr, ft, cparam; nodes=nodes, symmetrize=symmetrize)
  end

  mul(x::Vector{T}) where T<:ComplexF64 = produ(x,nrow,shape,λ,cparam,ft,ift)

  return LinearOperator{ComplexF64}(nrow, ncol, false, true
            , mul
            , nothing
            , nothing )
end

function produ(x::Vector{T}, numOfNodes::Int, shape::Tuple, λ::Matrix, cparam::InhomogeneityData, ft, ift) where T<:ComplexF64
  K = size(cparam.A_k,2)
  s = zeros(ComplexF64,numOfNodes)
  # Preprocessing step when time and correctionMap are centered
  # x = x .* exp(-vec(cparam.Cmap) * cparam.t_hat )
  s = @distributed (+) for κ=1:K
                        p_tild = cparam.C_k[κ,:] .* x;
                        s_tild = toeplitzProd(p_tild, shape, λ[:,κ], ft, ift)
                        s_tild = conj.( cparam.C_k[κ,:] ) .* s_tild
                    end
  # Postprocessing step when time and correctionMap are centered
  # s = s .* exp.(-conj.(vec(cparam.Cmap)) * cparam.t_hat)

  return s
end

#
# TODO: optimize allocation of xˡᵃʳᵍᵉ
#
function toeplitzProd(x::Vector{T}, shape::Tuple, λ::Vector, fftplan, ifftplan) where T
  xˡᵃʳᵍᵉ = zeros(T,Tuple(2*collect(shape)))
  xˡᵃʳᵍᵉ[1:shape[1],1:shape[2]] = x
  λ = reshape(λ,Tuple(2*collect(shape)))

  return vec( (ifftplan*( λ.*(fftplan*xˡᵃʳᵍᵉ)))[1:shape[1],1:shape[2]] )

end

#
# calculate the matrix element A_{j,k} explicitely
#
function getMatrixElement(j::Int, k::Int, κ::Int, shape::Tuple, cparam::InhomogeneityData, nodes::Matrix; weights=nothing)
    elem=0.
    x = mod(k-1,shape[1])-mod(j-1,shape[1])
    y = div(k-1,shape[1])-div(j-1,shape[1])
    if weights != nothing
      for i=1:size(nodes,2)
        elem += cparam.A_k[i,κ]*exp( -2*pi*1im*(nodes[1,i]*x + nodes[2,i]*y) )*weights[i]
      end
    else
      for i=1:size(nodes,2)
        elem += cparam.A_k[i,κ]*exp( -2*pi*1im*(nodes[1,i]*x + nodes[2,i]*y) )
      end
    end

    return elem
end

function toeplitzEigvals(κ::Int, shape::Tuple, tr::AbstractTrajectory, ft, cparam::InhomogeneityData; nodes=nothing, symmetrize=true,m=2,σ=1.25)

  # calculate first column of block Toeplitz matrix
  nodes==nothing ? nodes=kspaceNodes(tr) : nothing
  p = NFFTPlan(nodes, shape, m, σ)
  density = sdc(p)

  if symmetrize
    firstCol = [getMatrixElement(j,1,κ,shape,cparam,nodes,weights=density) for j=1:prod(shape)]
  else
    firstCol = [getMatrixElement(j,1,κ,shape,cparam,nodes) for j=1:prod(shape)]
  end
  firstCol = reshape(firstCol, shape)

  # calculate first rows of the leftmost Toeplitz blocks
  firstRow = zeros(ComplexF64,shape[2],shape[1])
  if symmetrize
     for i=1:shape[2]
       firstRow[i,:] = [ getMatrixElement((i-1)*shape[1]+1,k,κ,shape,cparam,nodes,weights=density) for k=1:shape[1] ] # first row of the relevant blocks
     end
   else
     for i=1:shape[2]
       firstRow[i,:] = [ getMatrixElement((i-1)*shape[1]+1,k,κ,shape,cparam,nodes) for k=1:shape[1] ] # first row of the relevant blocks
     end
   end

  # construct FT matrix of the eigentvalues
  eigMat = zeros(ComplexF64, Tuple(2*collect(shape)))
  eigMat[1:shape[1], 1:shape[2]] = firstCol
  eigMat[shape[1]+2:end, 1:shape[2]] = copy(transpose(firstRow[:, end:-1:2]))
  eigMat[1:shape[1], shape[2]+2:end] = firstRow[end:-1:2,:]'
  eigMat[shape[1]+2:end, shape[2]+2:end] =  conj.(firstCol[end:-1:2,end:-1:2])

  return vec(ft*eigMat)
end


function lowRankParams(numOfNodes::Int64
                                , numOfPixel::Int64
                                , times::Vector
                                , cmap::Matrix
                                , m::Float64 = 4.0
                                , alpha::Float64=1.75)
                                # , K::Int64=20)

  cmap_sum = vec( [conj(cmap[i])+cmap[j] for i=1:length(cmap), j=1:length(cmap)] )

  N = Int64(8*ceil( maximum(abs.(cmap)) * maximum(abs.(times)) / (2*pi) ))
  K = ceil(Int64, alpha*N + 2*m)

  A = getA_Coefficients_hist_lsqr(K,numOfNodes,numOfPixel,vec(times),vec(cmap_sum))
  C = getC_Coefficients_hist_lsqr(K,numOfNodes,numOfPixel,vec(times),vec(cmap))
  if size(A,2) != size(C,1)
    error("Consistency check failed! A and C are not compatible")
  end

  t_hat = (times[1] + times[end])/2
  z_hat = minimum(real.(cmap)) + maximum(real.(cmap))
  z_hat += 1im*(minimum(imag.(cmap)) + maximum(imag.(cmap)))
  z_hat *= 0.5

  return InhomogeneityData(A,C,vec(times),cmap,t_hat,z_hat,"hist")
end
