using .Utils: get_default_cartesian_raw_acq_data_2d, get_default_cartesian_raw_acq_data_3d


@testset "kDataCart restores the order when profiles are shuffled" begin
	n_profiles = 10
	n_samples_per_profile = 12
	data = rand(ComplexF32, n_samples_per_profile, n_profiles)
	f = get_default_cartesian_raw_acq_data_2d(n_profiles, n_samples_per_profile, data)
	f.profiles = shuffle(f.profiles)
	acq = AcquisitionData(f)
	kdata = kDataCart(acq)[:, :, 1, 1, 1, 1]
	@test isapprox(kdata, data)
end


@testset "kDataCart restores the order when profiles are shuffled for 3D data" begin
	n_profiles = 10
	n_samples_per_profile = 12
	n_slices = 3
	data = rand(ComplexF32, n_samples_per_profile, n_profiles, n_slices)
	elipsoid_mask = zeros(Bool, 1, n_profiles, n_slices)
	for i in 1:n_profiles
		for j in 1:n_slices
			elipsoid_mask[1, i, j] = (i-5)^2/25 + (j-2)^2/4 < 0.5
		end
	end
	data .*= elipsoid_mask
	f = get_default_cartesian_raw_acq_data_3d(n_slices, n_profiles, n_samples_per_profile, data, elipsoid_mask)
	f.profiles = shuffle(f.profiles)
	acq = AcquisitionData(f)
	kdata = kDataCart(acq)[:, :, :, 1, 1, 1]

	@test isapprox(kdata, data)
end


@testset "profileData returns the sampled profiles of a slice or partition" begin
    kspace = rand(ComplexF32, 8, 6, 5, 2, 2, 3)  # x, y, z or slices, coils, contrasts, repetitions
    undersampled = copy(kspace)
    undersampled[:, 2, 1, :, 2, :] .= 0
    @test MRIBase.profileData(AcquisitionData(kspace; enc2D = true), 2, 4, 3, 5) == kspace[:, 5, 4, :, 2, 3]
    @test MRIBase.profileData(AcquisitionData(kspace), 2, 4, 3, 5) == kspace[:, 5, 4, :, 2, 3]
    # in contrast 2, partition 1 lacks profile 2, so its second sampled profile is profile 3
    @test MRIBase.profileData(AcquisitionData(undersampled), 2, 1, 3, 2) == kspace[:, 3, 1, :, 2, 3]
    @test MRIBase.profileData(AcquisitionData(undersampled), 2, 2, 3, 2) == kspace[:, 2, 2, :, 2, 3]
end
