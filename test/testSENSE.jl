

#=function test_SENSE(N::Int64=64,reductionFactor::Int64=2)
    result = MRILib.test_cartesianSENSE(N,1,reductionFactor,4)
    @test result.NRMSE < 0.5
end=#
