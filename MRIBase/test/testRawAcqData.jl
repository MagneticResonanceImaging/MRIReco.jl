using Random

using MRIBase


function get_default_cartesian_raw_acq_data(n_profiles, n_samples_per_profile, data=nothing)
    params = Dict{String, Any}()
    params["trajectory"] = "Cartesian"
    params["encodedSize"] = [n_profiles, n_samples_per_profile]
    
    traj_nodes = kspaceNodes(trajectory(Float32, "Cartesian", n_profiles, n_samples_per_profile))
    traj_nodes = reshape(traj_nodes, 2, n_samples_per_profile, n_profiles)
    if data == nothing
        data = rand(ComplexF32, n_samples_per_profile, n_profiles)
    end

    profiles = Vector{Profile}(undef, n_profiles)
    for i_profile in 1:n_profiles
        profile_head = AcquisitionHeader(idx=EncodingCounters(kspace_encode_step_1=i_profile-1))
        profiles[i_profile] = Profile(
            profile_head, 
            traj_nodes[:, :, i_profile], 
            reshape(data[:, i_profile], n_samples_per_profile, 1)
        )
    end

    f = RawAcquisitionData(
        params,
        profiles,
    )
    return f
end


@testset "MRIBase.rawdata with dense 2D cartesian data" begin
    n_profiles = 10
    n_samples_per_profile = 12
    data = rand(ComplexF32, n_samples_per_profile, n_profiles)
    f = get_default_cartesian_raw_acq_data(n_profiles, n_samples_per_profile, data)
    kdata = MRIBase.rawdata(f)
    kdata = reshape(kdata, n_samples_per_profile, n_profiles)
    @test isapprox(kdata, data)
end


@testset "MRIBase.rawdata with 2D cartesian data and shuffled profiles" begin
    n_profiles = 10
    n_samples_per_profile = 12
    data = rand(ComplexF32, n_samples_per_profile, n_profiles)
    f = get_default_cartesian_raw_acq_data(n_profiles, n_samples_per_profile, data)
    f.profiles = shuffle(f.profiles)
    kdata = MRIBase.rawdata(f)
    kdata = reshape(kdata, n_samples_per_profile, n_profiles)
    @test isapprox(kdata, data)
end