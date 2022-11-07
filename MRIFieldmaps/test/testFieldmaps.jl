
function testFieldmap(N=32,slices=1)
    
    TEs = [4.5, 8.2] #[ms]
    off_res = quadraticFieldmap(N,N,300.0)
    p1 = off_res*TEs[1]/1000 #[rad]
    p2 = off_res*TEs[2]/1000 #[rad]

    sim_im_pic = complex(1 .+ shepp_logan(N))
    sim_im = repeat(sim_im_pic,outer=(1,1,slices,2,1))

    for slice = 1:slices

        sim_im[:,:,slice,1,1] = sim_im[:,:,slice,1,1] .* exp.(1im .* (p1 .+ 0.1*randn(N,N)))
        sim_im[:,:,slice,2,1] = sim_im[:,:,slice,2,1] .* exp.(1im .* (p2 .+ 0.1*randn(N,N)))

    end

    # perform cmap calculation
    cmaps = estimate_cmap(sim_im,1:slices,TEs[1],TEs[2],true; β = 0.5,reltol=1e-4)

    # check that cmap calculation is low error
    @test nrmsd(cmaps[:,:,1],off_res[:,:,1]) < 0.1 

    # TODO: More realistic fieldmap and image scenario than modified sl (sl + 1) and a quadratic field map

end

@testset "Fieldmaps" begin

testFieldmap()

end
