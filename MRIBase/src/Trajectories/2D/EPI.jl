export EPITrajectory, epiNodes, epiDensity

"""
    EPITrajectory(T, numProfiles, numSamplingPerProfile
                  ; TE::T=0.0
                  , AQ::T=1.e-3
                  , EPI_factor::Int64=1
                  , profileOffset= :equispaced
                  , kargs...)

returns a 2d cartesian trajectory.

# Arguments
* `numProfiles::Int64`            - number of profiles
* `numSamplingPerProfile::Int64`  - number of sampling points per profile
* (`TE::T=0.0`)             - echo time in s
* (`AQ::T=1.e-3`)           - readout duration in s (per profile)
* (`EPI_factor::Int64=1`)         - EPI factor, e.g. how many profiles to acquire per shot
* (`profileOffset= :equispaced`)  - equispaced or random ordering of the shots
"""
function EPITrajectory(::Type{T}, numProfiles, numSamplingPerProfile
                  ; TE=0.0
                  , AQ=1.e-3
                  , EPI_factor::Int64=1
                  , profileOffset= :equispaced
                  , kargs...) where T
  nodes = epiNodes(T, numProfiles, numSamplingPerProfile; EPI_factor=EPI_factor, profileOffset=profileOffset)
  times = readoutTimes(T, numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ)
  return Trajectory("EPI", nodes, times, T(TE), T(AQ), numProfiles, numSamplingPerProfile, 1, true, false)
end

function epiNodes(::Type{T}, numProfiles::Int64
                    , numSamplingPerProfile::Int64
                    ; EPI_factor::Int64=1
                    , profileOffset= :equispaced
                    , kargs...) where T
  nodes = zeros(T, 2, numSamplingPerProfile, numProfiles)

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

function epiDensity(::Type{T}, numProfiles::Int64, numSamplingPerProfile::Int64) where T
  @error "not implement yet"
end
