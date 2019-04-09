export EPITrajectory, epiNodes, epiDensity

function EPITrajectory(numProfiles, numSamplingPerProfile
                  ; TE::Float64=0.0
                  , AQ::Float64=1.e-3
                  , EPI_factor::Int64=1
                  , profileOffset= :equispaced
                  , kargs...)
  nodes = epiNodes(numProfiles, numSamplingPerProfile; EPI_factor=EPI_factor, profileOffset=profileOffset)
  times = readoutTimes(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ)
  return Trajectory("EPI", nodes, times, TE, AQ, numProfiles, numSamplingPerProfile, 1, true, false)
end

function epiNodes(numProfiles::Int64
                    , numSamplingPerProfile::Int64
                    ; EPI_factor::Int64=1
                    , profileOffset= :equispaced
                    , kargs...)
  nodes = zeros(2,numSamplingPerProfile, numProfiles)

  if profileOffset == :equispaced
      samplesperline = numSamplingPerProfile/EPI_factor -1
      for l = 1:numProfiles
        for gamma1 = 1:EPI_factor
          for gamma2 = 1:numSamplingPerProfile/EPI_factor
            # Node index
            index = Int64( (gamma1-1) * (numSamplingPerProfile/EPI_factor) + gamma2 )
            # Nodes containing kx values, second index is the actual sample index, l profile index
            nodes[1,index,l] =  (-1)^(gamma1-1) * (gamma2-1) / samplesperline + ((gamma1-1) % 2 ) - 0.5
            # Nodes containing ky values, second index is the actual sample index, l profile index
            nodes[2,index,l] = ( (l-1) + (gamma1-1)*numProfiles)/(numProfiles*EPI_factor) - 0.5
          end
        end
      end
  elseif profileOffset == :random
      samplesperline = numSamplingPerProfile/EPI_factor -1
      for l = 1:numProfiles
        # Sorting the random samples
        angleOffsets = sort(rand(EPI_factor)-0.5)
        for gamma1 = 1:EPI_factor
          for gamma2 = 1:numSamplingPerProfile/EPI_factor
            # Node index
            index = Int64( (gamma1-1) * (numSamplingPerProfile/EPI_factor) + gamma2 )
            # Nodes containing kx values, second index is the actual sample index, l profile index
            nodes[1,index,l] =  (-1)^(gamma1-1) * (gamma2-1) / samplesperline + ((gamma1-1) % 2 ) - 0.5
            # Nodes containing ky values, second index is the actual sample index, l profile index
            nodes[2,index,l] = angleOffsets[gamma1]
          end
        end
      end
  end
  return reshape(nodes, 2, numSamplingPerProfile*numProfiles)
end

function epiDensity(numProfiles::Int64, numSamplingPerProfile::Int64)
  @error "not implement yet"
end
