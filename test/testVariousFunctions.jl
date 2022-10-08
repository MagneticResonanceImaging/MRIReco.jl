@testset "VariousFunctions" begin
    # test mergeChannels
    N = 32
    data = zeros(Float32,N,N,N,1,2,4);
    data[:,:,:,1,:,1] .= sqrt(2);
    dataMerge = mergeChannels(data);

    @test ndims(dataMerge) == 6
    @test abs(dataMerge[1,1,1,1,1,1] - 2) < 1e-6
    @test dataMerge[1,1,1,1,1,2] == 0


end
