export RawAcquisitionDataFid_360

function RawAcquisitionDataFid_360(b::BrukerFile)
  T = Complex{MRIFiles.acqWordSize(b)}
  filename = joinpath(b.path, "rawdata.job0")

  N = pvmMatrix(b)
  factor_AntiAlias = pvmAntiAlias(b)
  
  N = round.(Int,N .* factor_AntiAlias)

  if length(N) < 3
    N_ = ones(Int,3)
    N_[1:length(N)] .= N
    N = N_
  end

  numChannel = parse.(Int,b["PVM_EncNReceivers"])
  numAvailableChannel = MRIFiles.pvmEncAvailReceivers(b)
  numSlices = MRIFiles.acqNumSlices(b)
  numEchos = MRIFiles.acqNumEchos(b)
  phaseFactor = MRIFiles.acqPhaseFactor(b) 
  numRep = MRIFiles.acqNumRepetitions(b)
  numEncSteps = MRIFiles.acqSpatialSize1(b)
  readoutLength =  parse.(Int,b["PVM_EncMatrix"])[1]

  centerSample = parse.(Int,b["PVM_EncPftOverscans"])[1]
  profileLength = readoutLength*numChannel#MRIFiles.acqSize(b)[1] #Int((ceil(N[1]*numChannel*sizeof(dtype)/1024))*1024/sizeof(dtype)) # number of points + zeros

  I = open(filename,"r") do fd
    read!(fd,Array{T,6}(undef, profileLength,
                                   numEchos,
                                   phaseFactor,
                                   numSlices,
                                   div(numEncSteps, phaseFactor),
                                   numRep))
  end

  encSteps1 = parse.(Int,b["PVM_EncGenSteps1"]).+floor(Int,N[2]/2)
  encSteps2 = parse.(Int,b["PVM_EncGenSteps2"]).+floor(Int,N[3]/2)

  objOrd = MRIFiles.acqObjOrder(b)
  objOrd = objOrd.-minimum(objOrd)

  gradMatrix = MRIFiles.acqGradMatrix(b)

  offset1 = MRIFiles.acqReadOffset(b)
  offset2 = MRIFiles.acqPhase1Offset(b)

  offset3 = ndims(b) == 2 ? MRIFiles.acqSliceOffset(b) : MRIFiles.pvmEffPhase2Offset(b)

 
  profiles = Profile[]
  for nR = 1:numRep
    for nEnc = 1:div(numEncSteps, phaseFactor)
      for nSl = 1:numSlices
        for nPhase = 1:phaseFactor
          for nEcho=1:numEchos
              counter = EncodingCounters(kspace_encode_step_1=encSteps1[nPhase+phaseFactor*(nEnc-1)],
                                          kspace_encode_step_2=encSteps2[nPhase+phaseFactor*(nEnc-1)],
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
                                        center_sample=centerSample,#div(N[1],2),
                                        available_channels = numChannel, #numAvailableChannel ?
                                        active_channels = numChannel)


              traj = Matrix{Float32}(undef,0,0)
              dat = map(T, reshape(I[:,nEcho,nPhase,nSl,nEnc,nR],:,numChannel))
              p = Profile(head,traj,dat)

              if (b["EchoAcqMode"] == "allEchoes") && mod(nEcho,2) == 0 # set reverse specific flags for MGE
                @info nEcho
              flag_set!(p,"ACQ_IS_REVERSE")
              end

              push!(profiles,p)

              
          end
        end
      end
    end
  end

  params = MRIFiles.brukerParams(b)
  params["trajectory"] = "cartesian"
  
  params["encodedSize"] = N

  return RawAcquisitionData(params, profiles)
end