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
  tr = CartesianTrajectory(N,N)
  F_nfft = NFFTOp((N,N),tr,symmetrize=false)
  F_exp = ExplicitOp((N,N),tr,zeros(ComplexF64,N,N),symmetrize=false)

  # test agains FourierOperators
  y = vec( ifftshift(reshape(F*vec(fftshift(x)),N,N)) )
  y_adj = vec( ifftshift(reshape(F_adj*vec(fftshift(x)),N,N)) )

  y_nfft = F_nfft*vec(x)
  y_adj_nfft = adjoint(F_nfft) * vec(x)

  y_exp = F_exp*vec(x)
  y_adj_exp = adjoint(F_exp) * vec(x)

  @test (norm(y-y_nfft)/norm(y)) < 1e-2
  @test (norm(y_adj-y_adj_nfft)/norm(y_adj)) < 1e-2
  @test (norm(y-y_exp)/norm(y)) < 1e-2
  @test (norm(y_adj-y_adj_exp)/norm(y_adj)) < 1e-2
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
  tr = CartesianTrajectory3D(N,N,numSlices=N)
  F_nfft = NFFTOp((N,N,N),tr,symmetrize=false)
  F_exp = ExplicitOp((N,N,N),tr,zeros(ComplexF64,N,N,N),symmetrize=false)

  # test agains FourierOperators
  y = vec( ifftshift(reshape(F*vec(fftshift(x)),N,N,N)) )
  y_adj = vec( ifftshift(reshape(F_adj*vec(fftshift(x)),N,N,N)) )

  y_nfft = F_nfft*vec(x)
  y_adj_nfft = adjoint(F_nfft) * vec(x)

  y_exp = F_exp*vec(x)
  y_adj_exp = adjoint(F_exp) * vec(x)

  relError = (norm(y-y_exp)/norm(y))
  println("Relative error ExplicitOp: ", relError)
  @test  relError< 1e-2
  relError = (norm(y_adj-y_adj_exp)/norm(y_adj))
  println("Relative error adj. ExplicitOp: ", relError)
  @test  relError< 1e-2
  relError = (norm(y-y_nfft)/norm(y))
  println("Relative error NFFT: ", relError)
  @test relError < 1e-2
  relError = (norm(y_adj-y_adj_nfft)/norm(y_adj))
  println("Relative error adj. NFFT: ", relError)
  @test  relError< 1e-2
end

# test FieldmapNFFTOp
function testFieldmapFT(N=16)
  # random image
  x = zeros(ComplexF64,N,N)
  for i=1:N,j=1:N
    x[i,j] = rand()
  end

  tr = CartesianTrajectory(N,N;TE=0.0,AQ=0.01)
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

function testOperators()
  @testset "Linear Operator" begin
    testFT()
    testFT3d()
    testFieldmapFT()
  end
end

testOperators()
