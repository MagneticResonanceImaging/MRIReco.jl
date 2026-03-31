"""
Minimal test algorithm for testing context parameters.
"""
@reconstruction mutable struct TestAlgorithm{P} <: AbstractIterativeMRIRecoAlgorithm
  @parameter parameter::P
end

"""
Minimal test parameters that track calls and capture context.
"""
@parameter struct TestIterativeParameters <: AbstractIterativeRecoParameters
  alloc_called::Ref{Bool} = Ref(false)
  loop_called::Ref{Int} = Ref(0)
  finalize_called::Ref{Bool} = Ref(false)
  captured_reconSize::Ref{Any} = Ref(nothing)
  captured_acqData::Ref{Any} = Ref(nothing)
  arrayType::Type = Array
  lock::ReentrantLock = ReentrantLock()
end

function (params::TestIterativeParameters)(
  algo::TestAlgorithm, 
  reconSize::NTuple{D, Int64}
) where D
  params.alloc_called[] = true
  params.captured_reconSize[] = reconSize
  
  acqData = ctx_acqData()
  params.captured_acqData[] = acqData
  
  numRep = numRepetitions(acqData)
  numSl = numSlices(acqData)
  
  Ireco = zeros(ComplexF64, prod(reconSize), numSl, numRep)
  indices = CartesianIndices((numRep, numSl))
  weights = ones(Float64, 1)
  
  return (Ireco, indices, weights)
end

function (params::TestIterativeParameters)(
  algoT::Type{<:TestAlgorithm},
  Ireco::Array{Complex{T}, 3},
  index::CartesianIndex{2}, 
  weights
) where T

  lock(params.lock) do
    params.loop_called[] += 1
  end
  # Can access context during iteration
  reconSize = ctx_reconSize()
  acqData = ctx_acqData()
  arrayType = ctx_arrayType()

  # can update Ireco
  Ireco[1, index[2], index[1]] = complex(params.loop_called[])
end

function (params::TestIterativeParameters)(
  algo::TestAlgorithm, 
  Ireco::Array{Complex{T}, 3}
) where T
  params.finalize_called[] = true
  
  acqData = ctx_acqData()
  numSl = numSlices(acqData)
  numRep = numRepetitions(acqData)
  
  # Return as AxisArray
  return reshape(Ireco, ctx_reconSize()..., numSl, numRep)
end


# Actual test functions

function testIterativeContextParametersSerial()
  @testset "Serial Iterative Context" begin
    # Create acquisition data
    x = shepp_logan(32)

    # simulation
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Cartesian"
    params[:numProfiles] = floor(Int64, N)
    params[:numSamplingPerProfile] = N

    acqData = simulation(x, params)
    numSl = numSlices(acqData)
    numRep = numRepetitions(acqData)
    
    # Create test parameters
    params = TestIterativeParameters()
    
    # Create context parameter
    contextParams = SerialIterativeMRIRecoContextParameter(
      reconSize = shape,
      parameter = params
    )
    
    algo = TestAlgorithm(contextParams)
    
    # Run reconstruction
    result = reconstruct(algo, acqData)
    
    # Verify allocation was called
    @test params.alloc_called[]
    
    # Verify loop was called correct number of times (numRep * numSl)
    @test params.loop_called[] == numRep * numSl
    
    # Verify finalization was called
    @test params.finalize_called[]
    
    # Verify reconSize was captured
    @test params.captured_reconSize[] == shape
    
    # Verify acqData was captured
    @test params.captured_acqData[] === acqData
    
    # Verify result is AxisArray
    @test size(result, 5) == numRep  # rep dimension
    @test !iszero(result)
  end
end

function testIterativeContextParametersThreaded()
  @testset "Threaded Iterative Context" begin
    # Create acquisition data
    x = shepp_logan(32)

    # simulation
    params = Dict{Symbol, Any}()
    params[:simulation] = "fast"
    params[:trajName] = "Cartesian"
    params[:numProfiles] = floor(Int64, N)
    params[:numSamplingPerProfile] = N

    acqData = simulation(x, params)
    numSl = numSlices(acqData)
    numRep = numRepetitions(acqData)

    # Test different schedulers
    schedulers = [
      :serial,
      :dynamic,
      MRIReco.OhMyThreads.SerialScheduler(),
      MRIReco.OhMyThreads.DynamicScheduler()
    ]
    
    for sched in schedulers
      @testset "Scheduler: $sched" begin
        # Create fresh test parameters
        params = TestIterativeParameters()
        
        # Create threaded context parameter
        contextParams = ThreadedIterativeMRIRecoContextParameter(
          reconSize = shape,
          parameter = params,
          scheduler = sched
        )
        
        algo = TestAlgorithm(contextParams)
        
        # Run reconstruction
        result = reconstruct(algo, acqData)
        
        # Verify allocation was called
        @test params.alloc_called[]
        
        # Verify loop was called correct number of times (numRep * numSl)
        @test params.loop_called[] == numRep * numSl
        
        # Verify finalization was called
        @test params.finalize_called[]
        
        # Verify result is AxisArray
        @test size(result, 5) == numRep  # rep dimension
        @test !iszero(result)
      end
    end
  end
end

@testset "Iterative Context Parameters" begin
  testIterativeContextParametersSerial()
  testIterativeContextParametersThreaded()
end