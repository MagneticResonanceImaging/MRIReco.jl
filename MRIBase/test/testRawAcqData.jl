using Random
using MRIBase
using .Utils: get_default_cartesian_raw_acq_data_2d


@testset "MRIBase.rawdata with dense 2D cartesian data" begin
	n_profiles = 10
	n_samples_per_profile = 12
	data = rand(ComplexF32, n_samples_per_profile, n_profiles)
	f = get_default_cartesian_raw_acq_data_2d(n_profiles, n_samples_per_profile, data)
	kdata = MRIBase.rawdata(f)
	kdata = reshape(kdata, n_samples_per_profile, n_profiles)
	@test isapprox(kdata, data)
end


@testset "MRIBase.rawdata is not supposed to handle shuffles profiles" begin
	n_profiles = 10
	n_samples_per_profile = 12
	data = rand(ComplexF32, n_samples_per_profile, n_profiles)
	f = get_default_cartesian_raw_acq_data_2d(n_profiles, n_samples_per_profile, data)
    shuffle!(f.profiles)
	kdata = MRIBase.rawdata(f)
	kdata = reshape(kdata, n_samples_per_profile, n_profiles)
	@test !isapprox(kdata, data)
end
