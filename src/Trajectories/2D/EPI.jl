export EPITrajectory

mutable struct EPITrajectory <: Abstract2DTrajectory
  numProfiles::Int
  numSamplingPerProfile::Int
  TE::Real # echo time in ms
  AQ::Real # time for each spiral arm in ms
  EPI_factor # EPI factor, integer greater one
  profileOffset::Symbol
end


function EPITrajectory(numProfiles::Int64
                    , numSamplingPerProfile::Int64
                    ; TE::Real=0.0
                    , AQ::Real=1.0e-3
                    , EPI_factor::Int64=1
                    , profileOffset= :equispaced
                    , kargs...)

   EPITrajectory(numProfiles
                , numSamplingPerProfile
                , TE
                , AQ
                , EPI_factor
                , profileOffset
                )
end
string(tr::EPITrajectory) = "EPI"

function kspaceNodes(tr::EPITrajectory)
  nodes = zeros(2,tr.numSamplingPerProfile, tr.numProfiles)

  if tr.profileOffset == :equispaced
      samplesperline = tr.numSamplingPerProfile/tr.EPI_factor -1
      for l = 1:tr.numProfiles
        for gamma1 = 1:tr.EPI_factor
          for gamma2 = 1:tr.numSamplingPerProfile/tr.EPI_factor
            # Node index
            index = Int64( (gamma1-1) * (tr.numSamplingPerProfile/tr.EPI_factor) + gamma2 )
            # Nodes containing kx values, second index is the actual sample index, l profile index
            nodes[1,index,l] =  (-1)^(gamma1-1) * (gamma2-1) / samplesperline + ((gamma1-1) % 2 ) - 0.5
            # Nodes containing ky values, second index is the actual sample index, l profile index
            nodes[2,index,l] = ( (l-1) + (gamma1-1)*tr.numProfiles)/(tr.numProfiles*tr.EPI_factor) - 0.5
          end
        end
      end
  elseif tr.profileOffset == :random
      samplesperline = tr.numSamplingPerProfile/tr.EPI_factor -1
      for l = 1:tr.numProfiles
        # Sorting the random samples
        angleOffsets = sort(rand(tr.EPI_factor)-0.5)
        for gamma1 = 1:tr.EPI_factor
          for gamma2 = 1:tr.numSamplingPerProfile/tr.EPI_factor
            # Node index
            index = Int64( (gamma1-1) * (tr.numSamplingPerProfile/tr.EPI_factor) + gamma2 )
            # Nodes containing kx values, second index is the actual sample index, l profile index
            nodes[1,index,l] =  (-1)^(gamma1-1) * (gamma2-1) / samplesperline + ((gamma1-1) % 2 ) - 0.5
            # Nodes containing ky values, second index is the actual sample index, l profile index
            nodes[2,index,l] = angleOffsets[gamma1]
          end
        end
      end
  end
  return reshape(nodes, 2, tr.numSamplingPerProfile*tr.numProfiles)
end

function kspaceDensity(tr::EPITrajectory)

end
