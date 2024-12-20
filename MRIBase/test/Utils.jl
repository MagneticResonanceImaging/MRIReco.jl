module Utils
using MRIBase

export get_default_cartesian_raw_acq_data

function get_default_cartesian_raw_acq_data_2d(n_profiles, n_samples_per_profile, data = nothing)
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


function get_default_cartesian_raw_acq_data_3d(n_slices, n_profiles, n_samples_per_profile, data = nothing, mask=nothing)
	params = Dict{String, Any}()
	params["trajectory"] = "Cartesian"
	params["encodedSize"] = [n_samples_per_profile, n_profiles, n_slices]
	params["encodedFOV"] = [n_samples_per_profile, n_profiles, n_slices]

	traj_nodes = kspaceNodes(trajectory(Float32, "Cartesian3D", n_profiles, n_samples_per_profile, numSlices=n_slices))
	traj_nodes = reshape(traj_nodes, 3, n_samples_per_profile, n_profiles, n_slices)
	if data == nothing
		data = rand(ComplexF32, n_samples_per_profile, n_profiles, n_slices)
	end

	n_elemens =  mask == nothing ? n_profiles*n_slices : sum(mask)
	profiles = Vector{Profile}(undef, n_elemens)
	i = 1
	for i_profile in 1:n_profiles
		for i_slice in 1:n_slices
			if mask != nothing && !mask[1, i_profile, i_slice]
				continue
			end
		
			profile_head = AcquisitionHeader(
				idx = EncodingCounters(
					kspace_encode_step_1 = i_profile - 1,
					kspace_encode_step_2 = i_slice - 1,
					), 
				center_sample = div(n_samples_per_profile, 2),
				trajectory_dimensions=3,
			)
			# profiles[(i_slice-1)*n_profiles+i_profile] = Profile(
			profiles[i] = Profile(
				profile_head,
				traj_nodes[:, :, i_profile, i_slice],
				reshape(data[:, i_profile, i_slice], n_samples_per_profile, 1),
			)
			i += 1
		end
	end

	f = RawAcquisitionData(
		params,
		profiles,
	)
	return f
end

end
