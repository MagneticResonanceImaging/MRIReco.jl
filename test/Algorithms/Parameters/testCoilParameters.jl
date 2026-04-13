function testCoilParameters(arrayType=Array, T=Float64)
  @testset "Coil Parameters $arrayType, $T" begin

    @testset "Construction" begin
      # Default constructor
      @test_throws ArgumentError CoilParameters()

      # Keyword arguments
      @test_throws ArgumentError CoilParameters(; senseMaps=nothing, noiseData=nothing)

      # With senseMaps only
      smaps = rand(Complex{T}, 8, 8, 1, 2)  # x, y, slice, coil
      cp2 = CoilParameters(; senseMaps=smaps)
      @test cp2.senseMaps === smaps
      @test isempty(cp2.noiseData)

      # With noiseData only
      ndata = rand(Complex{T}, 100, 100)
      cp3 = CoilParameters(; noiseData=ndata)
      @test isempty(cp3.senseMaps)
      @test cp3.noiseData === ndata

      # With both
      cp4 = CoilParameters(smaps, ndata)
      @test cp4.senseMaps === smaps
      @test cp4.noiseData === ndata
    end

    @testset "Coil Parameters with Context $arrayType, $T" begin
      N = 32
      Tc = Complex{T}
      x = shepp_logan(N)
      rmap = 20.0*ones(N,N)
      coilsens = Tc.(measured2DSensitivity(N, 8))
      params = Dict{Symbol, Any}()
      params[:simulation] = "fast"
      params[:trajName] = "Cartesian"
      params[:numProfiles] = floor(Int64, N)
      params[:numSamplingPerProfile] = N
      params[:r2map] = rmap
      params[:T_echo] = [2.e-2, 4.e-2]
      params[:seqName] = "ME"
      params[:refocusingAngles] = Float64[pi,pi]
      params[:senseMaps] = coilsens
      acqData = simulation( real(x), params )

      # Tests dont actually use real sense map data, but 
      # in case I've simulated an acqData from a sense test
      senseMaps = reshape(coilsens,N,N,1,8)
      shape = encodingSize(acqData)
      numChan = numChannels(acqData)

      ctx = MRIRecoContext(acqData, arrayType)

      with(MRIRECO_CONTEXT => ctx) do
        @testset "With senseMaps only" begin
          smaps = rand(Complex{T}, shape[1], shape[2], 1, numChan)
          cp = CoilParameters(; senseMaps = smaps)
          L_inv, senseMapsUnCorr = cp(AbstractIterativeMRIRecoAlgorithm)

          @test isnothing(L_inv)
          @test isapprox(senseMapsUnCorr, arrayType(smaps))
          @test senseMapsUnCorr isa arrayType
        end

        @testset "With noiseData only - creates decorrelation" begin
          ndata = rand(Complex{T}, 100, numChan)

          cp = CoilParameters(; noiseData=ndata)
          L_inv, senseMapsUnCorr = cp(AbstractIterativeMRIRecoAlgorithm)

          @test L_inv !== nothing
          @test isempty(senseMapsUnCorr)

          # L_inv should be lower triangular
          @test L_inv isa LowerTriangular
        end

        @testset "With both senseMaps and noiseData - decorrelates" begin
          smaps = rand(Complex{T}, shape[1], shape[2], 1, numChan)

          # Create noise samples
          ndata = rand(Complex{T}, 100, numChan)

          cp = CoilParameters(smaps, ndata)
          L_inv, senseMapsUnCorr = cp(AbstractIterativeMRIRecoAlgorithm)

          @test L_inv !== nothing
          @test !isempty(senseMapsUnCorr)
          @test size(senseMapsUnCorr) == size(smaps)
          @test senseMapsUnCorr isa arrayType
        end
      end
    end
  end
end

for arrayType in arrayTypes
  for T in [Float32, Float64]
    testCoilParameters(arrayType, T)
  end
end