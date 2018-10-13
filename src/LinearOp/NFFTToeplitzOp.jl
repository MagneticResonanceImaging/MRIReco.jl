export NFFTToeplitzOp

mutable struct NFFTToeplitzOp{T,F1<:FuncOrNothing,F2<:FuncOrNothing,F3<:FuncOrNothing} <:
                      AbstractLinearOperator{T,F1,F2,F3}
  nrow :: Int
  ncol :: Int
  symmetric :: Bool
  hermitian :: Bool
  prod :: F1
  tprod :: F2
  ctprod :: F3
  density
end

#
# Linear Operator to perform NFFT
#
function NFFTToeplitzOp(shape::Tuple, tr::AbstractTrajectory; nodes=nothing, symmetrize=true,m=2,σ=1.25)

  # precalculate eigenvalues and the corresponding FFTs
  λ, ft, ift = diagonalizeOp(shape, tr; nodes=nodes, symmetrize=symmetrize, m=m, σ=σ)

  return NFFTToeplitzOp{ComplexF64}(prod(shape), prod(shape), false, true
          , x->produ(x, shape, λ, ft, ift)
          , nothing
          , nothing
          , density )
end

function produ(x::Vector{T}, shape::Tuple, λ::Vector, fftplan, ifftplan) where T
  xˡᵃʳᵍᵉ = zeros(T,Tuple(2*collect(shape)))
  xˡᵃʳᵍᵉ[1:shape[1],1:shape[2]] = x
  λ = reshape(λ,Tuple(2*collect(shape)))

  return vec( (ifftplan*( λ.*(fftplan*xˡᵃʳᵍᵉ)))[1:shape[1],1:shape[2]] )

end

#
# calculate the matrix element A_{j,k} explicitely
#
function getMatrixElement(j::Int, k::Int, shape::Tuple, nodes::Matrix; weights=nothing)
    elem=0.
    x = mod(k-1,shape[1])-mod(j-1,shape[1])
    y = div(k-1,shape[1])-div(j-1,shape[1])
    if weights != nothing
      for i=1:size(nodes,2)
        elem += exp( -2*pi*1im*(nodes[1,i]*x + nodes[2,i]*y) )*weights[i]
      end
    else
      for i=1:size(nodes,2)
        elem += exp( -2*pi*1im*(nodes[1,i]*x + nodes[2,i]*y) )
      end
    end

    return elem
end

function diagonalizeOp(shape::Tuple, tr::AbstractTrajectory; nodes=nothing, symmetrize=true,m=2,σ=1.25)

  # plan the FFTs
  fftplan = plan_fft( zeros(ComplexF64, Tuple(2*collect(shape)));flags=FFTW.MEASURE )
  ifftplan = plan_ifft(zeros(ComplexF64, Tuple(2*collect(shape)));flags=FFTW.MEASURE)

  # calculate first column of block Toeplitz matrix
  nodes==nothing ? nodes=kspaceNodes(tr) : nothing
  p = NFFTPlan(nodes, shape, m, σ)
  density = convert(Vector{Float64}, sdc(p))

  if symmetrize
    firstCol = [getMatrixElement(j,1,shape,nodes,weights=density) for j=1:prod(shape)]
  else
    firstCol = [getMatrixElement(j,1,shape,nodes) for j=1:prod(shape)]
  end
  firstCol = reshape(firstCol, shape)

  # TODO: more efficient implementation using NFFTs
  # e1 = zeros(ComplexF64,shape)
  # e1[1] = 1.
  # if symmetrize
  #   firstCol = ndft_adjoint(p, ndft(p,e1)[:].*density)
  # else
  #   firstCol = ndft_adjoint(p, ndft(p,e1))
  # end
  # firstCol = reshape(firstCol, shape)

  # calculate first rows of the leftmost Toeplitz blocks
  firstRow = zeros(ComplexF64,shape[2],shape[1])
   if symmetrize
     for i=1:shape[2]
       firstRow[i,:] = [ getMatrixElement((i-1)*shape[1]+1,k,shape,nodes,weights=density) for k=1:shape[1] ] # first row of the relevant blocks
     end
   else
     for i=1:shape[2]
       firstRow[i,:] = [ getMatrixElement((i-1)*shape[1]+1,k,shape,nodes) for k=1:shape[1] ] # first row of the relevant blocks
     end
   end

  # firstRow = zeros(ComplexF64,shape[2],shape[1])
  # if symmetrize
  #   for i=1:shape[2]
  #     e1 = zeros(ComplexF64,shape)
  #     e1[1,i] = 1.
  #     firstRow[i,:] = conj( ndft_adjoint(p, ndft(p,e1)[:].*density)[1:shape[1]] )
  #   end
  # else
  #   for i=1:shape[2]
  #     e1 = zeros(ComplexF64,shape)
  #     e1[1,i] = 1.
  #     firstRow[i,:] = conj( ndft_adjoint(p, ndft(p,e1)[1:shape[1]]) )
  #   end
  # end

  # construct FT matrix of the eigentvalues
  eigMat = zeros(ComplexF64, Tuple(2*collect(shape)))
  eigMat[1:shape[1], 1:shape[2]] = firstCol
  eigMat[shape[1]+2:end, 1:shape[2]] = copy(transpose(firstRow[:, end:-1:2]))
  eigMat[1:shape[1], shape[2]+2:end] = firstRow[end:-1:2,:]'
  eigMat[shape[1]+2:end, shape[2]+2:end] =  conj(firstCol[end:-1:2,end:-1:2])

  return vec(fftplan*eigMat), fftplan, ifftplan
end
