export SpiralPerturbedTrajectory, spiralPerturbedTrajectoryNodes

function SpiralPerturbedTrajectory(numProfiles, numSamplingPerProfile
                  ; TE::Real=0.0
                  , AQ::Real=1.e-3
                  , windings::Real= 6.25
                  , alpha::Real=2.0
                  , kmax::Real=0.5
                  , angleOffset::String="equispaced"
                  , kargs...)
  nodes = spiralPerturbedTrajectoryNodes(numProfiles, numSamplingPerProfile;
            windings=windings, alpha=alpha, angleOffset=angleOffset, kmax=kmax)
  times = readoutTimes(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ)
  return  Trajectory("SpiralPerturbed", nodes, times, TE, AQ, numProfiles, numSamplingPerProfile, 1, false, true)
end

# Constructing the amplitude limited case of the variable density spiral Trajectory
function spiralPerturbedTrajectoryNodes(numProfiles::Int64
                        , numSamplingPerProfile::Int64
                        ; windings=6.25
                        , alpha=2.0
                        , kmax::Real=0.5
                        , angleOffset= "equispaced"
                        , kargs...)

  nodes = zeros(2,numSamplingPerProfile, numProfiles)
  angles = collect((0:numProfiles-1)/numProfiles) # Phase offset for multiple Profiles

  if angleOffset == "golden"
      angles = [i*(3-sqrt(5))/2  for i=0:numProfiles-1 ]
  elseif angleOffset == "random"
      angles = sort(rand(numProfiles))
  elseif angleOffset == "equispaced"
      angles = collect((0:numProfiles-1)/numProfiles)
  end


  T_ea = 1/((alpha+1)*1/(2pi*windings))
  times = range(0,stop=T_ea,length=numSamplingPerProfile)
  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      tau = (1/(2pi*windings) *(alpha+1)*times[k])^(1/(alpha+1))
      tmp = (1+0.025*cos(k*0.1))*(tau^alpha)*(kmax-0.025)*exp(1im*2*pi*( windings*tau + angles[l] ))
      # Thresholding if the petrubed spiral goes beyond |x_j| > [0.5,0.5]
      if real(tmp) >= 0.5
          nodes[1,k,l] = 0.5 - 1e-6
      elseif real(tmp) < -0.5
          nodes[1,k,l] = -0.5
      else
          nodes[1,k,l] = real(tmp)
      end

      if imag(tmp) >= 0.5
          nodes[2,k,l] = 0.5 - 1e-6
      elseif imag(tmp) < -0.5
          nodes[2,k,l] = -0.5
      else
          nodes[2,k,l] = imag(tmp)
      end
    end
  end
  return reshape(nodes, 2, numSamplingPerProfile*numProfiles)
end
