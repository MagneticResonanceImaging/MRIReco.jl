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


@testset "AcquisitionData reverses profiles without changing the raw data" begin
    n_profiles = 10
    n_samples_per_profile = 12
    data = rand(ComplexF32, n_samples_per_profile, n_profiles)
    f = get_default_cartesian_raw_acq_data_2d(n_profiles, n_samples_per_profile, data)
    flag_set!(f.profiles[2], "ACQ_IS_REVERSE")
    acq = AcquisitionData(f)
    kdata = reshape(acq.kdata[1, 1, 1], n_samples_per_profile, n_profiles)
    @test kdata[:, 2] == reverse(data[:, 2])
    @test flag_is_set(f.profiles[2], "ACQ_IS_REVERSE")
    @test AcquisitionData(f).kdata == acq.kdata
end


@testset "RawAcquisitionData(acqData) converts back to the same AcquisitionData" begin
    kspace = rand(ComplexF32, 8, 4, 3, 2, 2, 2)  # x, y, z or slices, coils, contrasts, repetitions
    undersampled = copy(kspace)
    undersampled[:, 2, 1, :, :, :] .= 0
    single_samples = kspace[1:1, :, :, :, :, :]  # one sample per profile
    for acq in (AcquisitionData(kspace; enc2D = true), AcquisitionData(kspace), AcquisitionData(undersampled),
                AcquisitionData(single_samples))
        raw = RawAcquisitionData(acq)
        @test raw.params["encodedSize"] == [get(acq.encodingSize, d, 1) for d in 1:3]
        @test raw.params["receiverChannels"] == 2
        @test raw.params["trajectory"] == "cartesian"
        # each profile stores the nodes of its samples, the same for every contrast, slice and repetition
        nodes = kspaceNodes(trajectory(acq, 1))[:, acq.subsampleIndices[1]]
        @test reduce(hcat, [p.traj for p in raw.profiles]) == repeat(nodes, 1, length(acq.kdata))
        acq2 = AcquisitionData(raw)
        @test isCartesian(trajectory(acq2, 1))
        @test acq2.kdata == acq.kdata
        @test acq2.subsampleIndices == acq.subsampleIndices
        @test acq2.encodingSize == acq.encodingSize
    end
end


@testset "AcquisitionHeader(acqData, ...) fills the encoding counters" begin
    acq = AcquisitionData(rand(ComplexF32, 8, 4, 1, 1, 2, 1))
    # every argument has its own value, so that each counter can be told apart
    idx = AcquisitionHeader(acq, 9, 3, 4, 2, 5, 20; phase = 6, set = 7, segment = 8, user = ntuple(Int16, 8)).idx
    @test (idx.kspace_encode_step_1, idx.kspace_encode_step_2, idx.average, idx.slice) == (4, 3, 0, 2)
    @test (idx.contrast, idx.phase, idx.repetition, idx.set, idx.segment) == (1, 5, 8, 6, 7)
    @test idx.user == ntuple(identity, 8)
    @test AcquisitionHeader(acq, 9, 3, 4, 2, 5, 20; repetition = 1).idx.repetition == 0
end
