using MRIOperators.Adapt

# CPU runs move between Array and JLArray; GPU runs (arrayTypes set by test/gpu/*.jl) between Array and the device array.
adaptTargets = arrayTypes == [Array] ? [Array, JLArray] : arrayTypes
# Element types of the tests below; test/gpu/metal.jl sets ComplexF32 only, as Metal has no Float64.
adaptEltypes = @isdefined(adaptEltypes) ? adaptEltypes : (ComplexF32, ComplexF64)
adaptDirections = vcat([(Array, to) for to in adaptTargets], [(to, Array) for to in adaptTargets if to !== Array])
adaptRtol(T) = real(T) == Float32 ? 5e-5 : 1e-11

storageType(op) = MRIOperators.LinearOperators.storage_type(op)
testData(T, dims...) = reshape(T[T <: Real ? cos(0.3k + 0.1) : complex(cos(0.3k + 0.1), sin(0.7k) + 0.5) for k in 1:prod(dims)], dims)
onTarget(x, to) = x isa to && x isa DenseArray

# Forward and adjoint results with three- and five-argument `mul!`, collected on the CPU.
function adaptActions(op, conv, x, y)
  T = eltype(op)
  α, β = T(0.5), T(0.75)
  f3, f5, a3, a5 = conv(zeros(T, size(op, 1))), conv(copy(y)), conv(zeros(T, size(op, 2))), conv(copy(x))
  mul!(f3, op, conv(x)); mul!(f5, op, conv(x), α, β)
  mul!(a3, adjoint(op), conv(y)); mul!(a5, adjoint(op), conv(y), α, β)
  return map(Array, (f3, f5, a3, a5))
end
adaptActions(op, conv) = adaptActions(op, conv, testData(eltype(op), size(op, 2)), testData(eltype(op), size(op, 1)) .+ 1)

function testAdaptActions(moved, to, ref, refConv, rtol)
  for (m, r) in zip(adaptActions(moved, to), adaptActions(ref, refConv))
    @test isapprox(m, r; rtol)
  end
end

function testMovedWrapper(moved, src, to)
  @test typeof(moved).name.wrapper === typeof(src).name.wrapper
  @test moved !== src && parent(moved) !== parent(src)
  @test eltype(moved) == eltype(src) && size(moved) == size(src)
  @test storageType(parent(moved)) <: to
  @test moved.nrow == src.nrow && moved.ncol == src.ncol # property forwarding to the parent
end
function testSensitivityLeaves(op, to)
  @test fieldnames(typeof(op)) == (:op, :sensMaps, :sensMapsC, :numContr)
  @test all(isconcretetype, fieldtypes(typeof(op)))
  @test onTarget(getfield(op, :sensMaps), to)
  @test onTarget(getfield(op, :sensMapsC), to)
end
function testSubspaceLeaf(op, to)
  @test fieldnames(typeof(op)) == (:op, :basis, :shape, :numContr)
  @test all(isconcretetype, fieldtypes(typeof(op)))
  @test onTarget(getfield(op, :basis), to)
end

# Contents and counters of an operator's workspaces (`tmp` for LOC products).
workspaceState(op) = (Array(op.Mv), Array(op.Mtu), op.nprod, op.ntprod, op.nctprod,
                      hasfield(typeof(op), :tmp) ? Array(op.tmp) : nothing)

@testset "Adapt $(nameof(from)) -> $(nameof(to))" for (from, to) in adaptDirections
  @testset "Wrappers $T $shape numContr=$numContr" for T in adaptEltypes, shape in ((3, 2), (2, 2, 2)), numContr in (1, 3)
    rtol = adaptRtol(T)
    numVox, numChan = prod(shape), 2
    maps4 = testData(T, shape..., ntuple(_ -> 1, 3 - length(shape))..., numChan)
    maps = reshape(maps4, numVox, numChan)
    for src in (SensitivityOp(from(maps), numContr), SensitivityOp(from(maps4), numContr))
      moved = adapt(to, src)
      testMovedWrapper(moved, src, to)
      testSensitivityLeaves(moved, to)
      testAdaptActions(moved, to, src, from, rtol)
    end
    src = SubspaceOp(from(testData(T, numContr, 2)), shape, numContr)
    moved = adapt(to, src)
    testMovedWrapper(moved, src, to)
    testSubspaceLeaf(moved, to)
    @test getfield(moved, :shape) == shape && getfield(moved, :numContr) == numContr
    testAdaptActions(moved, to, src, from, rtol)
  end

  @testset "Nested product and source workspaces $T" for T in adaptEltypes
    shape, numContr = (2, 2, 2), 3
    maps, basis = from(testData(T, prod(shape), 2)), from(testData(T, numContr, 2))
    S = SensitivityOp(maps, numContr)
    B = SubspaceOp(basis, shape, numContr)
    W = WeightingOp(from(real(T)[0.5 + k / 10 for k in 1:size(S, 1)]))
    inner = ProdOp(S, B)
    P = ProdOp(W, inner)
    sources = (S, B, W, inner, P)
    # Exercising the sources allocates and writes their workspaces; the results are the references.
    refs = map(op -> adaptActions(op, from), sources)
    holders = (parent(S), parent(B), getfield(W, :op), inner, P)
    before = map(workspaceState, holders)
    inputs = map(Array, (maps, getfield(S, :sensMapsC), basis, W.weights))

    moved = adapt(to, P)
    @test moved.B.A isa SensitivityOp && moved.B.B isa SubspaceOp
    @test onTarget(moved.A.weights, to)
    testSensitivityLeaves(moved.B.A, to)
    testSubspaceLeaf(moved.B.B, to)
    for op in (moved, moved.A, moved.B, parent(moved.B.A), parent(moved.B.B))
      @test storageType(op) <: to
    end
    for (op, ref) in zip((moved.B.A, moved.B.B, moved.A, moved.B, moved), refs)
      for (m, r) in zip(adaptActions(op, to), ref)
        @test isapprox(m, r; rtol = adaptRtol(T))
      end
    end
    @test map(workspaceState, holders) == before
    @test map(Array, (maps, getfield(S, :sensMapsC), basis, W.weights)) == inputs
  end

  @testset "Sensitivity maps changed after construction" begin
    T = ComplexF32
    maps = testData(T, 4, 2)
    x, y = testData(T, 4), testData(T, 8)
    m = from(copy(maps))
    src = SensitivityOp(m, 1)
    m .*= 2 # the forward operation reads the live maps; the adjoint keeps the conjugate taken at construction
    srcFwd, srcAdj = Array(src * from(x)), Array(adjoint(src) * from(y))
    @test srcFwd ≈ 2 .* (SensitivityOp(maps, 1) * x)
    @test srcAdj ≈ adjoint(SensitivityOp(maps, 1)) * y
    @test !isapprox(Array(adjoint(SensitivityOp(m, 1)) * from(y)), srcAdj; rtol = 0.1) # rebuilding from the maps alone differs
    moved = adapt(to, src)
    @test isapprox(Array(moved * to(x)), srcFwd; rtol = adaptRtol(T))
    @test isapprox(Array(adjoint(moved) * to(y)), srcAdj; rtol = adaptRtol(T))
  end

  @testset "Fourier encoding with FieldmapNFFTOp" begin
    N = 16
    tr = CartesianTrajectory(Float64, N, N; TE = 0.002, AQ = 0.01)
    # times relative to the readout start and least-squares time segmentation, both not the defaults
    src = with(nfft_backend => NFFT.backend()) do
      fourierEncodingOp((N, N), tr, "fast"; correctionMap = im * quadraticFieldmap(N, N), echoImage = false,
                        method = "leastsquare", K = 4, subsampleIdx = collect(1:2:N^2), S = typeof(from(ComplexF64[])))
    end
    F = src.B
    @test F isa FieldmapNFFTOp && length(F.plans) > 1
    F.cparam.C_k .*= 2 # the stored coefficients must be used, not recomputed ones
    ref = adaptActions(src, from)
    moved = with(nfft_backend => NonuniformFFTs.backend()) do
      adapt(to, src)
    end
    G = moved.B
    @test storageType(moved) <: to && storageType(G) <: to
    @test all(a -> onTarget(a, to), (G.cparam.A_k, G.cparam.C_k, G.cparam.times, G.cparam.Cmap))
    @test G.nodes == F.nodes && all(k -> k isa Matrix, G.nodes) # NFFT plans take nodes on the CPU
    @test all(i -> i isa (to === Array ? Vector{Int} : to{Int32}), G.idx)
    @test (G.shape, G.circTraj, G.scheduler, G.backend) === (F.shape, F.circTraj, F.scheduler, F.backend)
    @test parentmodule(typeof(first(G.plans))) !== NonuniformFFTs
    testAdaptActions(moved, to, src, from, 1e-4)
    testAdaptActions(adapt(from, moved), from, src, from, 1e-4)
    for (s, r) in zip(adaptActions(src, from), ref)
      @test isapprox(s, r; rtol = 1e-10)
    end
  end

  @testset "MapSliceOp $T dim=$dim" for T in (Float32, ComplexF32), dim in (1, 2)
    size1, size2 = dim == 1 ? ((2, 4), (3, 4)) : ((4, 2), (4, 3))
    for trafo in (from(testData(T, 3, 2)), WeightingOp(from(testData(T, 2))))
      src = MapSliceOp(trafo, dim, size1, trafo isa AbstractMatrix ? size2 : size1; T)
      moved = adapt(to, src)
      testMovedWrapper(moved, src, to)
      @test storageType(getfield(moved, :trafo)) <: to
      testAdaptActions(moved, to, src, from, adaptRtol(T))
    end
  end

  @testset "ExplicitOp $(D)D echoImage=$echoImage" for D in (2, 3), echoImage in (false, true)
    T = ComplexF32
    tr = D == 2 ? CartesianTrajectory(Float32, 3, 3; TE = 0.001f0, AQ = 0.004f0) :
                  CartesianTrajectory3D(Float32, 3, 3; numSlices = 3, TE = 0.001f0, AQ = 0.004f0)
    src = ExplicitOp(ntuple(_ -> 3, D), tr, from(testData(T, ntuple(_ -> 3, D)...)); echoImage, S = typeof(from(T[])))
    moved = adapt(to, src)
    @test moved isa ExplicitOp && moved !== src && storageType(moved) <: to
    @test getfield(moved, :shape) == getfield(src, :shape) && getfield(moved, :echoOffset) == getfield(src, :echoOffset)
    @test all(f -> onTarget(getfield(moved, f), to), (:nodes, :times, :disturbanceTerm))
    @test eltype(getfield(moved, :nodes)) == eltype(getfield(moved, :times)) == Float32
    testAdaptActions(moved, to, src, from, adaptRtol(T))
    # an adjoint built before adapting, alone and inside a product
    adj = adjoint(src)
    P = ProdOp(adj, WeightingOp(from(testData(T, size(src, 1)))))
    movedAdj, movedP = adapt(to, adj), adapt(to, P)
    for E in (movedAdj.parent, movedP.A.parent)
      @test E isa ExplicitOp && storageType(E) <: to && onTarget(getfield(E, :nodes), to)
    end
    testAdaptActions(movedAdj, to, adj, from, adaptRtol(T))
    testAdaptActions(movedP, to, P, from, adaptRtol(T))
  end

  @testset "MapSliceOp in a product" begin
    T = ComplexF32
    src = ProdOp(MapSliceOp(from(testData(T, 3, 2)), 2, (4, 2), (4, 3); T), SensitivityOp(from(testData(T, 4, 2)), 1))
    moved = adapt(to, src)
    @test moved.A isa MapSliceOp && moved.B isa SensitivityOp
    @test all(op -> storageType(op) <: to, (moved, parent(moved.A), parent(moved.B)))
    testAdaptActions(moved, to, src, from, adaptRtol(T))
  end
end

if arrayTypes == [Array]
  @testset "Adapt SubspaceOp with explicit storage type $T" for T in (ComplexF32, ComplexF64)
    basis = testData(T, 3, 2)
    src = SubspaceOp(basis, (2, 2), 3, JLArray{T, 1})
    ref = SubspaceOp(basis, (2, 2), 3)
    @test storageType(parent(src)) == JLArray{T, 1}
    testAdaptActions(src, JLArray, ref, Array, adaptRtol(T))
    for to in (Array, JLArray)
      moved = adapt(to, src) # the storage type follows the adapted basis
      testMovedWrapper(moved, src, to)
      testSubspaceLeaf(moved, to)
      testAdaptActions(moved, to, ref, Array, adaptRtol(T))
    end
  end
end
