function testEncodingParameters(arrayType=Array, T=Float64)
  @testset "Encoding Parameters $arrayType, $T" begin
    
    @testset "Type definitions" begin
      @test AbstractMRIRecoEncodingParameters <: AbstractMRIRecoParameters
      @test EncodingParameters <: AbstractMRIRecoEncodingParameters
      @test CustomEncodingParameters <: AbstractMRIRecoEncodingParameters
    end

    @testset "Construction" begin
      # EncodingParameters - default
      ep = EncodingParameters()
      @test ep isa EncodingParameters
      @test isnotinhg(ep.correctionMap)
      @test isnotinhg(ep.method)
      @test isnotinhg(ep.toeplitz)
      @test isnotinhg(ep.oversamplingFactor)
      @test isnotinhg(ep.kernelSize)
      @test isnotinhg(ep.K)
      @test isnotinhg(ep.K_tol)

      # EncodingParameters - with correctionMap
      cmap = rand(Complex{T}, 8, 8, 8)
      ep_with_cmap = EncodingParameters(; correctionMap=cmap)
      @test ep_with_cmap.correctionMap === cmap

      # EncodingParameters - with some options
      ep_with_opts = EncodingParameters(; 
        method="fast",
        kernelSize=8,
        K=5
      )
      @test ep_with_opts.method == "fast"
      @test ep_with_opts.kernelSize == 8
      @test ep_with_opts.K == 5
      @test isnothing(ep_with_opts.correctionMap)

      # CustomEncodingParameters - with encoding operators
      ops = [opEye(Complex{T}, 64) for _ in 1:4]
      cep = CustomEncodingParameters(; encodingOps=ops)
      @test cep.encodingOps === ops
    end
  end
end

function testEncodingParametersWithContext(arrayType=Array, T=Float64)
  @testset "Encoding Parameters with Context $arrayType, $T" begin
    shape = (8, 8, 8)
    
    # Create sample acquisition data for testing
    acqData = acquisitionData([:cartesian], 1, shape[1:2]; numChannels=1, numSlices=1)
    
    ctx = MRIRecoContext(shape, acqData, arrayType)
    algoT = AbstractIterativeMRIRecoAlgorithm

    with(MRIRECO_CONTEXT => ctx) do
      
      @testset "EncodingParameters - default" begin
        ep = EncodingParameters()
        encOps = ep(algoT, 1)
        
        @test encOps isa Vector
        @test length(encOps) > 0
      end

      @testset "EncodingParameters - with method" begin
        ep = EncodingParameters(; method="fast")
        encOps = ep(algoT, 1)
        
        @test encOps isa Vector
      end

      @testset "EncodingParameters - with correctionMap" begin
        cmap = rand(Complex{T}, shape...)
        ep = EncodingParameters(; correctionMap=cmap)
        encOps = ep(algoT, 1)
        
        @test encOps isa Vector
      end

      @testset "CustomEncodingParameters" begin
        # Create mock encoding operators
        n = prod(shape)
        ops = [opEye(Complex{T}, n) for _ in 1:numContrasts(acqData)]
        
        cep = CustomEncodingParameters(; encodingOps=ops)
        result = cep(algoT, 1)
        
        @test result isa LinearOperator
      end
    end
  end
end

for arrayType in arrayTypes
  for T in [Float32, Float64]
    testEncodingParameters(arrayType, T)
    testEncodingParametersWithContext(arrayType, T)
  end
end