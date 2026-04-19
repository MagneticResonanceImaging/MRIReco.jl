function testEncodingParameters(arrayType=Array, T=Float64)
  @testset "Encoding Parameters $arrayType, $T" begin

    @testset "Construction" begin
      # EncodingParameters - default
      ep = EncodingParameters()
      @test ep isa EncodingParameters
      @test isnothing(ep.correctionMap)
      @test isnothing(ep.method)
      @test isnothing(ep.toeplitz)
      @test isnothing(ep.oversamplingFactor)
      @test isnothing(ep.kernelSize)
      @test isnothing(ep.K)
      @test isnothing(ep.K_tol)

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

    # simulation
    N = 32
    x = T.(shepp_logan(N))
    params = Dict{Symbol,Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Cartesian"
    params[:numProfiles] = floor(Int64, N)
    params[:numSamplingPerProfile] = N
    acqData = simulation(x, params)

    ctx = MRIRecoContext(acqData, arrayType)
    shape = encodingSize(acqData)
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
        ops = [rand(Complex{T}, n, n) for _ in 1:MRIBase.numContrasts(acqData)]
        cep = CustomEncodingParameters(; encodingOps=ops)
        result = cep(algoT, 1)

        @test first(result) isa Matrix
      end
    end
  end
end

for arrayType in arrayTypes
  for T in [Float64] # can't get Float32 simulation to work
    testEncodingParameters(arrayType, T)
  end
end