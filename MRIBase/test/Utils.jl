module Utils
using MRIBase

export get_default_cartesian_raw_acq_data

function get_default_cartesian_raw_acq_data(n_profiles, n_samples_per_profile, data = nothing)
	params = Dict{String, Any}()
	params["trajectory"] = "Cartesian"
	params["encodedSize"] = [n_samples_per_profile, n_profiles]
	params["encodedFOV"] = [n_samples_per_profile, n_profiles, 1]

	traj_nodes = kspaceNodes(trajectory(Float32, "Cartesian", n_profiles, n_samples_per_profile))
	traj_nodes = reshape(traj_nodes, 2, n_samples_per_profile, n_profiles)
	if data == nothing
		data = rand(ComplexF32, n_samples_per_profile, n_profiles)
	end

	profiles = Vector{Profile}(undef, n_profiles)
	for i_profile in 1:n_profiles
		profile_head = AcquisitionHeader(idx = EncodingCounters(kspace_encode_step_1 = i_profile - 1), center_sample = div(n_samples_per_profile, 2))
		profiles[i_profile] = Profile(
			profile_head,
			traj_nodes[:, :, i_profile],
			reshape(data[:, i_profile], n_samples_per_profile, 1),
		)
	end

	f = RawAcquisitionData(
		params,
		profiles,
	)
	return f
end
end
