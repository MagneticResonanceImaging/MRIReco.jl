function testSparsityParameters(arrayType=Array, T=Float64)
  @testset "Sparsity Parameters $arrayType, $T" begin

    @testset "SimpleSparsityParameters" begin
      # Default construction
      sp = SimpleSparsityParameters()
      @test sp.sparsity == ""

      # Single string
      sp = SimpleSparsityParameters("Wavelet")
      @test sp.sparsity == "Wavelet"

      # Vector of strings
      sp = SimpleSparsityParameters(["Wavelet", "Wavelet"])
      @test sp.sparsity == ["Wavelet", "Wavelet"]

      # Vector of strings and nothing
      sp = SimpleSparsityParameters(["Wavelet", nothing])
      @test sp.sparsity == ["Wavelet", nothing]

      # Nothing
      sp = SimpleSparsityParameters(nothing)
      @test isnothing(sp.sparsity)
    end

    @testset "CustomSparsityParameters" begin
      # Default construction
      cp = CustomSparsityParameters()
      @test cp.trafo === nothing

      # With operator
      op = opEye(T, 16)
      cp = CustomSparsityParameters(op)
      @test cp.trafo === op

      # With vector of operators
      op1 = opEye(T, 16)
      op2 = opEye(T, 16)
      cp = CustomSparsityParameters([op1, op2])
      @test length(cp.trafo) == 2

      # With vector of operators and nothing
      cp = CustomSparsityParameters([op1, nothing, op2])
      @test length(cp.trafo) == 3
    end

    @testset "Operator construction" begin
      shape = (8, 8, 8)
      ctx = MRIRecoContext(shape, T, arrayType)
      algoT = AbstractIterativeMRIRecoAlgorithm

      with(MRIRECO_CONTEXT => ctx) do

        @testset "SimpleSparsityParameters" begin
          # SimpleSparsityParameters: Empty string -> nothing (to be filled later)
          sp = SimpleSparsityParameters("")
          trafos = sp(algoT, 1)
          @test length(trafos) == 1
          @test isnothing(trafos[1])

          # SimpleSparsityParameters: Single string "Wavelet" -> SparseOp
          sp = SimpleSparsityParameters("Wavelet")
          trafos = sp(algoT, 1)
          @test length(trafos) == 1
          @test !isnothing(trafos[1])
          @test trafos[1] isa LinearOperator
          @test LinearOperators.storage_type(trafos[1]) == ctx_storageType()
          @test eltype(trafos[1]) == Complex{T}

          # SimpleSparsityParameters: Vector of strings with fill-up
          sp = SimpleSparsityParameters(["Wavelet"])
          trafos = sp(algoT, 3)
          @test length(trafos) == 3
          @test !isnothing(trafos[1])
          @test trafos[1] isa LinearOperator
          @test isnothing(trafos[2])
          @test isnothing(trafos[3])

          # SimpleSparsityParameters: Mixed vector (if supported)
          sp = SimpleSparsityParameters(["Wavelet", nothing, "Wavelet"])
          trafos = sp(algoT, 3)
          @test length(trafos) == 3
          @test !isnothing(trafos[1])
          @test isnothing(trafos[2])
          @test !isnothing(trafos[3])
          @test trafos[3] isa LinearOperator
        end

        @testset "CustomSparsityParameters" begin
          # CustomSparsityParameters: nothing -> nothing
          cp = CustomSparsityParameters(nothing)
          trafos = cp(algoT, 1)
          @test length(trafos) == 1
          @test isnothing(trafos[1])

          # CustomSparsityParameters: single operator
          op = opEye(ComplexF64, prod(shape))
          cp = CustomSparsityParameters(op)
          trafos = cp(algoT, 1)
          @test length(trafos) == 1
          @test trafos[1] === op

          # CustomSparsityParameters: vector of operators with fill-up
          op1 = opEye(ComplexF64, prod(shape))
          cp = CustomSparsityParameters([op1])
          trafos = cp(algoT, 3)
          @test length(trafos) == 3
          @test trafos[1] === op1
          @test isnothing(trafos[2])
          @test isnothing(trafos[3])

          # CustomSparsityParameters: string fallback to SparseOp
          cp = CustomSparsityParameters("Wavelet")
          trafos = cp(algoT, 1)
          @test length(trafos) == 1
          @test !isnothing(trafos[1])
        end
      end
    end
  end
end


function testRegularizationParameters(arrayType=Array, T=Float64)

  @testset "RegularizationParameters $arrayType, $T" begin
    shape = (8, 8, 8)
    ctx = MRIRecoContext(shape, T, arrayType)
    algoT = AbstractIterativeMRIRecoAlgorithm

    with(MRIRECO_CONTEXT => ctx) do
      S = ctx_storageType()

      @testset "Construction" begin
        # Default
        rp = RegularizationParameters()
        @test rp.reg === nothing
        @test rp.sparsity isa SimpleSparsityParameters

        # With single regularization
        reg = L1Regularization(1e-3)
        rp = RegularizationParameters(reg=reg)
        @test rp.reg === reg

        # With vector of regularizations
        regs = [L1Regularization(1e-3), L2Regularization(1e-4)]
        rp = RegularizationParameters(reg=regs)
        @test rp.reg == regs

        # With SimpleSparsityParameters
        rp = RegularizationParameters(reg=L1Regularization(1e-3), sparsity="Wavelet")
        @test rp.reg isa L1Regularization
        @test rp.sparsity.sparsity == "Wavelet"
        @test rp.sparsity isa SimpleSparsityParameters

        rp = RegularizationParameters(reg=[L1Regularization(1e-3), L1Regularization(1e-3)], sparsity=["Wavelet", "Wavelet"])
        @test rp.reg isa Vector
        @test rp.sparsity.sparsity == ["Wavelet", "Wavelet"]
        @test rp.sparsity isa SimpleSparsityParameters

        rp = RegularizationParameters(reg=[L1Regularization(1e-3), L1Regularization(1e-3)], sparsity=["Wavelet", nothing])
        @test rp.reg isa Vector
        @test rp.sparsity.sparsity == ["Wavelet", nothing]
        @test rp.sparsity isa SimpleSparsityParameters

        # With CustomSparsityParameters
        op = opEye(T, prod(shape); S=S)
        rp = RegularizationParameters(reg=L1Regularization(1e-3), sparsity=CustomSparsityParameters(op))
        @test rp.sparsity.trafo === op
      end

      @testset "Callable - simple reg" begin
        # No reg -> default L1(0)
        rp = RegularizationParameters()
        result = rp(algoT, FISTA)
        @test first(result) isa L1Regularization
        @test iszero(first(result).λ)

        # Single reg, nothing trafo -> returns reg directly
        rp = RegularizationParameters(reg=L1Regularization(1e-3), sparsity=nothing)
        result = rp(algoT, FISTA)
        @test first(result) isa L1Regularization
        @test length(result) == 1

        # Single reg with trafo -> returns TransformedRegularization
        rp = RegularizationParameters(reg=L1Regularization(1e-3), sparsity="Wavelet")
        result = rp(algoT, FISTA)
        @test first(result) isa TransformedRegularization
        @test RegularizedLeastSquares.innerreg(first(result)) isa L1Regularization
        @test length(result) == 1

        # Multiple regs -> returns vector of reg
        rp = RegularizationParameters(reg=[L1Regularization(1e-3), L2Regularization(1e-4)])
        result = rp(algoT, FISTA)
        @test result isa Vector
        @test result[1] isa L1Regularization
        @test result[2] isa L2Regularization
        @test length(result) == 2

        # Multiple regs with nothing trafo -> returns vector of reg
        rp = RegularizationParameters(reg=[L1Regularization(1e-3), L2Regularization(1e-4)], sparsity = [nothing, nothing])
        result = rp(algoT, FISTA)
        @test result isa Vector
        @test result[1] isa L1Regularization
        @test result[2] isa L2Regularization
        @test length(result) == 2

        # Multiple regs with trafo-> returns vector of reg
        rp = RegularizationParameters(reg=[L1Regularization(1e-3), L2Regularization(1e-4)], sparsity = ["Wavelet", "Wavelet"])
        result = rp(algoT, FISTA)
        @test result isa Vector
        @test result[1] isa TransformedRegularization
        @test RegularizedLeastSquares.innerreg(result[1]) isa L1Regularization
        @test result[2] isa TransformedRegularization
        @test RegularizedLeastSquares.innerreg(result[2]) isa L2Regularization
        @test length(result) == 2

        # Multiple regs with mixed trafo-> returns vector of reg
        rp = RegularizationParameters(reg=[L1Regularization(1e-3), L2Regularization(1e-4)], sparsity = ["Wavelet", nothing])
        result = rp(algoT, FISTA)
        @test result isa Vector
        @test result[1] isa TransformedRegularization
        @test RegularizedLeastSquares.innerreg(result[1]) isa L1Regularization
        @test result[2] isa L2Regularization
        @test length(result) == 2
      end

      @testset "Callable - reg, reg trafo" begin
        # No reg -> returns default reg and filled trafos (opEye)
        rp = RegularizationParameters()
        reg, trafos = rp(algoT, ADMM)
        @test length(reg) == 1
        @test length(trafos) == 1

        # Single reg, nothing trafo -> reg stays, trafo filled with opEye
        rp = RegularizationParameters(reg=L1Regularization(1e-3), sparsity=nothing)
        reg, trafos = rp(algoT, ADMM)
        @test length(reg) == 1
        @test reg[1] isa L1Regularization
        @test length(trafos) == 1

        # Single reg with trafo -> TransformedRegularization in reg
        rp = RegularizationParameters(reg=L1Regularization(1e-3), sparsity="Wavelet")
        reg, trafos = rp(algoT, ADMM)
        @test length(reg) == 1
        @test reg[1] isa L1Regularization
        @test length(trafos) == 1

        # Multiple regs with mixed trafos
        rp = RegularizationParameters(
          reg=[L1Regularization(1e-3), L2Regularization(1e-4)],
          sparsity=["Wavelet", nothing]
        )
        reg, trafos = rp(algoT, ADMM)
        @test length(reg) == 2
        @test reg[1] isa L1Regularization
        @test reg[2] isa L2Regularization
        @test length(trafos) == 2
      end
    end
  end
end

for arrayType in arrayTypes
  for T in [Float32, Float64]
    testSparsityParameters(arrayType, T)
    testRegularizationParameters(arrayType, T)
  end
end