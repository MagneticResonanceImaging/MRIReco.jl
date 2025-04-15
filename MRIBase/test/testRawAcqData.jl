function get_default_cartesian_raw_acq_data(n_profiles, n_samples_per_profile, data=nothing)
    params = Dict{String, Any}()
    params["trajectory"] = "Cartesian"
    params["encodedSize"] = [n_profiles, n_samples_per_profile]
    params["encodedFOV"] = [n_profiles, n_samples_per_profile, 1]
    
    traj_nodes = kspaceNodes(trajectory(Float32, "Cartesian", n_profiles, n_samples_per_profile))
    traj_nodes = reshape(traj_nodes, 2, n_samples_per_profile, n_profiles)

    profiles = Vector{Profile}(undef, n_profiles)
    for i_profile in 1:n_profiles
        idx = EncodingCounters(kspace_encode_step_1=i_profile-1)
        profile_head = AcquisitionHeader(idx=idx)
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


@testset "AcquisitionData averages signal" begin
    n_profiles = 10
    n_samples_per_profile = 12
    data = ones(ComplexF32, n_samples_per_profile, n_profiles) 
    f1 = get_default_cartesian_raw_acq_data(n_profiles, n_samples_per_profile, data.*0.5f0)
    f2 = get_default_cartesian_raw_acq_data(n_profiles, n_samples_per_profile, data.*1.5f0)
    for (f_n, file) in enumerate([f1, f2])
        for profile in file.profiles
            profile.head.idx.average = f_n-1
        end
    end
    raw_file = RawAcquisitionData(
        f1.params,
        vcat(f1.profiles, f2.profiles),
    )

    acq = AcquisitionData(raw_file)
    kdata = reshape(acq.kdata[1, 1, 1], n_samples_per_profile, n_profiles)
    @test isapprox(kdata, data)
end
