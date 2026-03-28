@testset "MRIRecoContext" begin

  @testset "Basic construction and accessors" begin
    ctx = MRIRecoContext((128, 128), Float64, Array)
    with(MRIRECO_CONTEXT => ctx) do
      @test ctx_reconSize() == (128, 128)
      @test ctx_arrayType() == Array
      @test ctx_storageType() == Vector{ComplexF64}
      @test ctx_acqData() === nothing
    end
  end

  @testset "Constructor with acqData" begin
    # Create acquisition data via simulation
    N = 32
    x = shepp_logan(N)

    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Cartesian"
    params[:numProfiles] = floor(Int64, N)
    params[:numSamplingPerProfile] = N

    acqData = simulation(x, params)

    # Constructor with acqData (3-arg version)
    ctx = MRIRecoContext((N, N), acqData, Array)
    with(MRIRECO_CONTEXT => ctx) do
      @test ctx_reconSize() == (N, N)
      @test ctx_arrayType() == Array
      @test ctx_storageType() == Vector{ComplexF64}
      @test ctx_acqData() === acqData
    end
  end

  @testset "Default values" begin
    # 2-arg constructor: (reconSize, eltype) uses default arrayType = Array
    ctx = MRIRecoContext((64, 64), Float64)
    with(MRIRECO_CONTEXT => ctx) do
      @test ctx_reconSize() == (64, 64)
      @test ctx_arrayType() == Array
      @test ctx_storageType() == Vector{ComplexF64}
    end
  end

  @testset "4-argument constructor" begin
    # Create acqData for valid 4-arg constructor
    N = 32
    x = shepp_logan(N)
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Cartesian"
    params[:numProfiles] = floor(Int64, N)
    params[:numSamplingPerProfile] = N
    acqData = simulation(x, params)

    # Valid 4-arg constructor
    ctx = MRIRecoContext((N, N), Array, Vector{ComplexF64}, acqData)
    with(MRIRECO_CONTEXT => ctx) do
      @test ctx_reconSize() == (N, N)
      @test ctx_arrayType() == Array
      @test ctx_storageType() == Vector{ComplexF64}
      @test ctx_acqData() === acqData
    end

    # 4-arg constructor with nothing as acqData
    ctx2 = MRIRecoContext((N, N), Array, Vector{ComplexF64}, nothing)
    with(MRIRECO_CONTEXT => ctx2) do
      @test ctx_reconSize() == (N, N)
      @test ctx_acqData() === nothing
    end
  end
end