export RawAcquisitionData_3DUTE

function RawAcquisitionData_3DUTE(b::BrukerFile)
    T = Complex{acqDataType(b)}

    filename = joinpath(b.path, "fid")
    filenameTraj = joinpath(b.path, "traj")

    N = acqSize(b)
    # The data is padded in case it is not a multiple of 1024 
    # For multi-channel acquisition data at concatenate then padded to a multiple of 1024 bytes
    numChannel = pvmEncNReceivers(b)
    profileLength = Int((ceil(N[1]*numChannel*sizeof(T)/1024))*1024/sizeof(T))
    numAvailableChannel = pvmEncAvailReceivers(b)
    phaseFactor = acqPhaseFactor(b)
    numSlices = acqNumSlices(b)
    numEchos = acqNumEchos(b)
    numEncSteps2 = length(N) == 3 ? N[3] : 1
    numRep = acqNumRepetitions(b)
    sample_time = parse(Float64,b["PVM_DigDw"])/1000 #us

    I = open(filename,"r") do fd
      read!(fd,Array{T,7}(undef, profileLength,
                                     numEchos,
                                     phaseFactor,
                                     numSlices,
                                     div(N[2], phaseFactor),
                                     numEncSteps2,
                                     numRep))[1:N[1]*numChannel,:,:,:,:,:,:]
    end

    traj = open(filenameTraj,"r") do fd
        read!(fd,Array{Float64,3}(undef, 3,N[1],
                                       N[2]))
    end

    objOrd = acqObjOrder(b)
    objOrd = objOrd.-minimum(objOrd)

    gradMatrix = acqGradMatrix(b)

    offset1 = acqReadOffset(b)
    offset2 = acqPhase1Offset(b)
    offset3 = ndims(b) == 2 ? acqSliceOffset(b) : pvmEffPhase2Offset(b)

    profiles = Profile[]
    for nR = 1:numRep
      for nEnc2 = 1:numEncSteps2
        for nPhase2 = 1:div(N[2], phaseFactor)
          for nSl = 1:numSlices
            for nPhase1 = 1:phaseFactor
              for nEcho=1:numEchos
                  counter = EncodingCounters(kspace_encode_step_1=0,
                                             kspace_encode_step_2=0,
                                             average=0,
                                             slice=objOrd[nSl],
                                             contrast=nEcho-1,
                                             phase=0,
                                             repetition=nR-1,
                                             set=0,
                                             segment=0 )

                  G = gradMatrix[:,:,nSl]
                  read_dir = (G[1,1],G[2,1],G[3,1])
                  phase_dir = (G[1,2],G[2,2],G[3,2])
                  slice_dir = (G[1,3],G[2,3],G[3,3])

                  # Not sure if the following is correct...
                  pos = offset1[nSl]*G[:,1] +
                        offset2[nSl]*G[:,2] +
                        offset3[nSl]*G[:,3]

                  position = (pos[1], pos[2], pos[3])

                  head = AcquisitionHeader(number_of_samples=N[1], idx=counter,
                                          read_dir=read_dir, phase_dir=phase_dir,
                                          slice_dir=slice_dir, position=position,
                                          center_sample=0,
                                          available_channels = numAvailableChannel, #TODO
                                          active_channels = numChannel,
                                          trajectory_dimensions = 3,
                                          sample_time_us = sample_time)
                  
                  dat = map(T, reshape(I[:,nEcho,nPhase1,nSl,nPhase2,nEnc2,nR],:,numChannel))
                  push!(profiles, Profile(head,convert.(Float32,traj[:,:,nPhase2]),dat) )
              end
            end
          end
        end
      end
    end

    params = brukerParams(b)
    params["trajectory"] = "custom"
    params["encodedSize"] = pvmMatrix(b)
    return RawAcquisitionData(params, profiles)
end
