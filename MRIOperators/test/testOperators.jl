
# test FourierOperators
function testFT(N=16)
  # random image
  x = zeros(ComplexF64,N,N)
  for i=1:N,j=1:N
    x[i,j] = rand()
  end

  # FourierMatrix
  idx = CartesianIndices((N,N))[collect(1:N^2)]
  F = [ exp(-2*pi*im*((idx[j][1]-1)*(idx[k][1]-1)+(idx[j][2]-1)*(idx[k][2]-1))/N) for j=1:N^2, k=1:N^2 ]
  F_adj = F'

  # Operators
  tr = CartesianTrajectory(Float64,N,N)
  F_nfft = NFFTOp((N,N),tr,symmetrize=false)
  F_exp = ExplicitOp((N,N),tr,zeros(ComplexF64,N,N),symmetrize=false)

  # test agains FourierOperators
  y = vec( ifftshift(reshape(F*vec(fftshift(x)),N,N)) )
  y_adj = vec( ifftshift(reshape(F_adj*vec(fftshift(x)),N,N)) )

  y_nfft = F_nfft*vec(x)
  y_adj_nfft = adjoint(F_nfft) * vec(x)

  y_exp = F_exp*vec(x)
  y_adj_exp = adjoint(F_exp) * vec(x)

  @test y     ≈ y_nfft      rtol = 1e-2
  @test y     ≈ y_exp       rtol = 1e-2
  @test y_adj ≈ y_adj_nfft  rtol = 1e-2
  @test y_adj ≈ y_adj_exp   rtol = 1e-2

  # test AHA w/o Toeplitz
  F_nfft.toeplitz = false
  AHA = SparsityOperators.normalOperator(F_nfft)
  y_AHA_nfft = AHA * vec(x)
  y_AHA = F' * F * vec(x)
  @test y_AHA ≈ y_AHA_nfft   rtol = 1e-2

  # test AHA with Toeplitz
  F_nfft.toeplitz = true
  AHA = SparsityOperators.normalOperator(F_nfft)
  y_AHA_nfft = AHA * vec(x)
  y_AHA_nfft = adjoint(F_nfft) * F_nfft * vec(x)
  y_AHA = F' * F * vec(x)
  @test y_AHA ≈ y_AHA_nfft   rtol = 1e-2

  # test type stability;
  # TODO: Ensure type stability for Trajectory objects and test here
  nodes = Float32.(tr.nodes)
  F_nfft = NFFTOp((N,N),nodes,symmetrize=false)

  y_nfft = F_nfft * vec(ComplexF32.(x))
  y_adj_nfft = adjoint(F_nfft) * vec(ComplexF32.(x))

  @test Complex{eltype(nodes)} === eltype(y_nfft)
  @test Complex{eltype(nodes)} === eltype(y_adj_nfft)
end

function testFT3d(N=12)
  # random image
  x = zeros(ComplexF64,N,N,N)
  for i=1:N,j=1:N,k=1:N
    x[i,j,k] = rand()
  end

  # FourierMatrix
  idx = CartesianIndices((N,N,N))[collect(1:N^3)]
  F = [ exp(-2*pi*im*((idx[j][1]-1)*(idx[k][1]-1)+(idx[j][2]-1)*(idx[k][2]-1)+(idx[j][3]-1)*(idx[k][3]-1))/N) for j=1:N^3, k=1:N^3 ]
  F_adj = F'

  # Operators
  tr = CartesianTrajectory3D(Float64,N,N,numSlices=N)
  F_nfft = NFFTOp((N,N,N),tr,symmetrize=false)
  F_exp = ExplicitOp((N,N,N),tr,zeros(ComplexF64,N,N,N),symmetrize=false)

  # test agains FourierOperators
  y = vec( ifftshift(reshape(F*vec(fftshift(x)),N,N,N)) )
  y_adj = vec( ifftshift(reshape(F_adj*vec(fftshift(x)),N,N,N)) )

  y_nfft = F_nfft*vec(x)
  y_adj_nfft = adjoint(F_nfft) * vec(x)

  y_exp = F_exp*vec(x)
  y_adj_exp = adjoint(F_exp) * vec(x)

  @test  y     ≈ y_exp      rtol = 1e-2
  @test  y     ≈ y_nfft     rtol = 1e-2
  @test  y_adj ≈ y_adj_exp  rtol = 1e-2
  @test  y_adj ≈ y_adj_nfft rtol = 1e-2

  # test AHA w/o Toeplitz
  F_nfft.toeplitz = false
  AHA = SparsityOperators.normalOperator(F_nfft)
  y_AHA_nfft = AHA * vec(x)
  y_AHA = F' * F * vec(x)
  @test y_AHA ≈ y_AHA_nfft   rtol = 1e-2

  # test AHA with Toeplitz
  F_nfft.toeplitz = true
  AHA = SparsityOperators.normalOperator(F_nfft)
  y_AHA_nfft = AHA * vec(x)
  y_AHA_nfft = adjoint(F_nfft) * F_nfft * vec(x)
  y_AHA = F' * F * vec(x)
  @test y_AHA ≈ y_AHA_nfft   rtol = 1e-2
end

function testUndersampledFourierOp(N=16)
  x = rand(ComplexF64,N,N) 
  tr = CartesianTrajectory(Float64, N÷2, N)
  
  # FourierOperator
  F_ft = fourierEncodingOp((N,N), tr, "fast")
  # Explicit operator
  F = ExplicitOp((N,N),tr,zeros(ComplexF64,N,N))
  
  y1 = F_ft*vec(x)
  y2 = F*vec(x)
  x1 = adjoint(F_ft)*y1
  x2 = adjoint(F)*y2
  @test y1 ≈ y2   rtol = 1e-2
  @test x1 ≈ x2   rtol = 1e-2
end

## test FieldmapNFFTOp
function testFieldmapFT(N=16)
  # random image
  x = zeros(ComplexF64,N,N)
  for i=1:N,j=1:N
    x[i,j] = rand()
  end

  tr = CartesianTrajectory(Float64,N,N;TE=0.0,AQ=0.01)
  times = readoutTimes(tr)
  nodes = kspaceNodes(tr)
  cmap = im*quadraticFieldmap(N,N)[:,:,1]

  # FourierMatrix
  idx = CartesianIndices((N,N))[collect(1:N^2)]
  F = [exp(-2*1im*pi*(nodes[1,k]*(idx[l][1]-size(x,1)/2-1)+nodes[2,k]*(idx[l][2]-size(x,2)/2-1))-cmap[idx[l][1],idx[l][2]]*times[k]) for k=1:size(nodes,2), l=1:length(x)]
  F_adj = F'

  # Operators
  F_nfft = FieldmapNFFTOp((N,N),tr,cmap,symmetrize=false)

  # test agains FourierOperators
  y = F*vec(x)
  y_adj = F_adj*vec(x)

  y_nfft = F_nfft*vec(x)
  y_adj_nfft = adjoint(F_nfft) * vec(x)

  @test (norm(y-y_nfft)/norm(y)) < 1e-2
  @test (norm(y_adj-y_adj_nfft)/norm(y_adj)) < 1e-2
end

function testSparseOp(T::Type,shape)
    W = SparseOp(Complex{T},"Wavelet", shape)

    x = zeros(Complex{T},shape)
    for i=1:shape[1],j=1:shape[2],k=1:shape[3]
      x[i,j,k] = rand()
    end
    x=vec(x)

    xapprox = W' * W * x
    @test (norm(xapprox-x)/norm(x)) < 1e-3
end

## test FieldmapNFFTOp
function testCopySizes(N=16)
  # random image
  x = zeros(ComplexF64,N,N)
  for i=1:N,j=1:N
    x[i,j] = rand()
  end

  tr = CartesianTrajectory(Float64,N,N;TE=0.0,AQ=0.01)
  times = readoutTimes(tr)
  nodes = kspaceNodes(tr)
  cmap = im*quadraticFieldmap(N,N)[:,:,1]

  # FourierMatrix
  idx = CartesianIndices((N,N))[collect(1:N^2)]

  # Operators

  F_nfft = NFFTOp((N,N),tr,symmetrize=false)
  F_fmap_nfft = FieldmapNFFTOp((N,N),tr,cmap,symmetrize=false)

  # Copy the FieldmapNFFTOp operator and change the plans field of the new operator to empty 
  F_fmap_nfft_copy = copy(F_fmap_nfft)
  F_nfft_copy = copy(F_nfft)

  @test sizeof(F_fmap_nfft_copy) == sizeof(F_fmap_nfft)

end

function testOperators()
  @testset "Linear Operator" begin
    testFT()
    testFT3d()
    testFieldmapFT()
    testUndersampledFourierOp()
    testSparseOp(Float32,(128,80,1))
    testSparseOp(Float64,(128,80,1))
    testSparseOp(Float32,(128,128,80))
    testSparseOp(Float64,(128,128,80))
  end
end

testOperators()
