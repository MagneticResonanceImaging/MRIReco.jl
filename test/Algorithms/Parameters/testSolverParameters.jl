function testSparsityOperators(arrayType=Array, T=ComplexF64)
  @testset "SolverParameters $arrayType, $T" begin

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

for arrayType in arrayTypes
  for T in [Float32, Float64]
    testSparsityOperators(arrayType, T)
  end
end