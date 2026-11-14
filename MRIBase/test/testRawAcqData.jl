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
    reversed = RawAcquisitionData(AcquisitionData(kspace))
    reverse!(reversed.profiles)  # stored in another order than the profile numbers
    for acq in (AcquisitionData(kspace; enc2D = true), AcquisitionData(kspace), AcquisitionData(undersampled),
                AcquisitionData(single_samples), AcquisitionData(reversed))
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


@testset "RawAcquisitionData(acqData) stores the nodes of the sampled profiles of radial data" begin
    tr = trajectory(Float32, "Radial", 6, 8)
    idx = findall(in((1, 3, 4, 6)), repeat(1:6, inner = 8))  # the samples of four of the six profiles
    acq = AcquisitionData(tr, reshape([rand(ComplexF32, 32, 2)], 1, 1, 1); idx = [idx], encodingSize = (8, 8))
    raw = RawAcquisitionData(acq)
    @test reduce(hcat, [p.traj for p in raw.profiles]) == kspaceNodes(tr)[:, idx]
end


@testset "RawAcquisitionData(acqData) converts binned data back to the same AcquisitionData" begin
    tr = trajectory(Float32, "Radial", 6, 8)
    raw = RawAcquisitionData(AcquisitionData(tr, reshape([rand(ComplexF32, 48, 2)], 1, 1, 1); encodingSize = (8, 8)))
    # stored in reverse order, the profiles 2 and 1 form the first contrast and the profiles 6 to 3 the second;
    # a "custom" trajectory only holds the profiles of its contrast, in this order
    reverse!(raw.profiles)
    for (i, profile) in enumerate(raw.profiles)
        profile.head.idx.contrast = i > 4 ? 0 : 1
    end
    raw.params["trajectory"] = "custom"
    acq = AcquisitionData(raw)
    acq2 = AcquisitionData(RawAcquisitionData(acq))
    @test acq2.kdata == acq.kdata
    @test kspaceNodes.(acq2.traj) == kspaceNodes.(acq.traj)
    @test acq2.subsampleIndices == acq.subsampleIndices
end
