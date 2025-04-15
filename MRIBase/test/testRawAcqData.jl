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


@testset "AcquisitionData averages signal" begin
    n_profiles = 10
    n_samples_per_profile = 12
    data = ones(ComplexF32, n_samples_per_profile, n_profiles) 
    f1 = get_default_cartesian_raw_acq_data_2d(n_profiles, n_samples_per_profile, data.*0.5f0)
    f2 = get_default_cartesian_raw_acq_data_2d(n_profiles, n_samples_per_profile, data.*1.5f0)
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