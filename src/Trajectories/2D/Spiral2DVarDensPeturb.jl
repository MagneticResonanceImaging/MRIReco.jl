export SpiralVarDensPerturbTrajectory, spiralVarDensPeturbNodes,spiralVarDensPeturbDensity

function SpiralVarDensPerturbTrajectory(numProfiles, numSamplingPerProfile
                  ; TE::Float64=0.0
                  , AQ::Float64=1.e-3
                  , windings::Float64= 6.25
                  , alpha=2.0
                  , angleOffset::String="equispaced"
                  , kargs...)
  nodes = spiralVarDensPerturbNodes(numProfiles, numSamplingPerProfile; windings=windings, alpha=alpha, angleOffset=angleOffset)
  times = readoutTimes(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ)
  return  Trajectory("SpiralVarDensPerturb", nodes, times, TE, AQ, numProfiles, numSamplingPerProfile, 1, false, true)
end

# Constructing the amplitude limited case of the variable density spiral Trajectory
function spiralVarDensPeturbNodes(numProfiles::Int64
                        , numSamplingPerProfile::Int64
                        ; windings=6.25
                        , alpha=2.0
                        , angleOffset= "equispaced"
                        , kargs...)

  nodes = zeros(2,numSamplingPerProfile, numProfiles)
  angles = collect((0:numProfiles-1)/numProfiles) # Phase offset for multiple Profiles
  A = 0.5 # Maximum radius is 0.5
  gamma = 42.58 *1e6 # 42,58 MHz/T : Gyromagnetic relation for protons

  if angleOffset == "golden"
      angles = [i*(3-sqrt(5))/2  for i=0:numProfiles-1 ]
  elseif angleOffset == "random"
      angles = sort(rand(numProfiles))
  elseif angleOffset == "equispaced"
      angles = collect((0:numProfiles-1)/numProfiles)
  end


  T_ea = 1/((alpha+1)*gamma/(2pi*windings))
  times = range(0,stop=T_ea,length=numSamplingPerProfile)
  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      tau = (gamma/(2pi*windings) *(alpha+1)*times[k])^(1/(alpha+1))
      tmp = (1+0.025*cos(k*0.1))*(tau^alpha)*(A-0.025)*exp(1im*2*pi*( windings*tau + angles[l] ))
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


function spiralVarDensPeturbDensity(numSamplingPerProfile::Int64, numProfiles::Int64)
  error("Not implemented!")
end
