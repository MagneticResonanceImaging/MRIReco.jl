function testWeightingParameters(arrayType=Array, T=Float64)
  @testset "Weighting Parameters $arrayType, $T" begin
    @testset "Construction" begin
      # DensityWeightingParameters - default
      dwp = DensityWeightingParameters()
      @test dwp isa DensityWeightingParameters

      # UniformWeightingParameters - default  
      uwp = UniformWeightingParameters()
      @test uwp isa UniformWeightingParameters

      # CustomWeightingParameters - with weights
      weights = Vector{Complex{T}}[[T(1), T(2)], [T(3), T(4)]]
      cwp = CustomWeightingParameters(; weights=weights)
      @test cwp.weights === weights

      # CustomWeightingParameters - empty weights
      cwp_empty = CustomWeightingParameters(; weights=Vector{Complex{T}}[])
      @test cwp_empty.weights == Vector{Complex{T}}[]
    end



    # simulation
    N = 32
    x = shepp_logan(N)
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Cartesian"
    params[:numProfiles] = floor(Int64, N)
    params[:numSamplingPerProfile] = N
    acqData = simulation(x, params)
    
    ctx = MRIRecoContext(acqData, arrayType)
    algoT = AbstractIterativeMRIRecoAlgorithm

    with(MRIRECO_CONTEXT => ctx) do
      shape = ctx_reconSize()

      @testset "DensityWeightingParameters" begin
        dwp = DensityWeightingParameters()
        weights = dwp(algoT)
        
        @test weights isa Vector
        @test length(weights) > 0
        @test all(w -> eltype(w) == eltype(ctx_storageType()), weights)
      end

      @testset "UniformWeightingParameters" begin
        uwp = UniformWeightingParameters()
        weights = uwp(algoT)
        
        @test weights isa Vector
        @test length(weights) > 0
        
        # Check uniform weight value
        expected_weight = T(1) / sqrt(prod(shape))
        for weight_vec in weights
          @test all(w -> abs(w - Complex{T}(expected_weight)) < T(1e-10), weight_vec)
        end
      end

      @testset "CustomWeightingParameters" begin
        custom_weights = Vector{Complex{T}}[randn(T, N * N)]
        cwp = CustomWeightingParameters(; weights=custom_weights)
        result = cwp(algoT)
        
        @test result === custom_weights
      end

      @testset "CustomWeightingParameters validation" begin
        # Get acqData to determine correct size
        acqData = ctx_acqData()
        numContr = MRIBase.numContrasts(acqData)
        
        # Test: wrong number of contrasts
        wrong_contr_weights = Vector{Complex{T}}[fill(Complex{T}(1), 2) for _ in 1:(numContr + 1)]
        cwp_wrong = CustomWeightingParameters(; weights=wrong_contr_weights)
        @test_throws DimensionMismatch cwp_wrong(algoT)
        
        # Test: wrong number of samples per contrast
        wrong_samples_weights = Vector{Complex{T}}[fill(Complex{T}(1), 100) for _ in 1:numContr]
        cwp_wrong2 = CustomWeightingParameters(; weights=wrong_samples_weights)
        @test_throws DimensionMismatch cwp_wrong2(algoT)
      end
    end
  end
end

for arrayType in arrayTypes
  for T in [Float32, Float64]
    testWeightingParameters(arrayType, T)
  end
end