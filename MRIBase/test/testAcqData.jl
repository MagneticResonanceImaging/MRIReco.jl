using .Utils: get_default_cartesian_raw_acq_data


@testset "kDataCart restores the order when profiles are shuffled" begin
	n_profiles = 10
	n_samples_per_profile = 12
	data = rand(ComplexF32, n_samples_per_profile, n_profiles)
	f = get_default_cartesian_raw_acq_data(n_profiles, n_samples_per_profile, data)
	f.profiles = shuffle(f.profiles)
	acq = AcquisitionData(f)
	kdata = kDataCart(acq)[:, :, 1, 1, 1, 1]
	@test isapprox(kdata, data)
end
