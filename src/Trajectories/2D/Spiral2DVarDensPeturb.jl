export SpiralTrajectoryVarDensPeturb

mutable struct SpiralTrajectoryVarDensPeturb <: Abstract2DTrajectory
  numProfiles::Int
  numSamplingPerProfile::Int
  TE::Float64 # echo time in ms
  AQ::Float64 # time for each spiral arm in ms
  alpha::Float64 # densyfing factor
  windings::Float64 # Number of turns the spiral is taking
  angleOffset::String
end

# Brought by:  Simple Analytic Variable Density Spiral Design
# by Dong-hyun Kim, 1,2 * Elfar Adalsteinsson, 1 and Daniel M. Spielman
function SpiralTrajectoryVarDensPeturb(numProfiles::Int64
                        , numSamplingPerProfile::Int64
                        ; TE=0.0
                        , AQ=1.0
                        , windings=6.25
                        , alpha=2.0
                        , angleOffset= "equispaced"
                        , kargs...)

    SpiralTrajectoryVarDensPeturb(numProfiles
                            , numSamplingPerProfile
                            , TE
                            , AQ
                            , alpha
                            , windings
                            , angleOffset
                            )
end

string(tr::SpiralTrajectoryVarDensPeturb) = "SpiralVarDensPeturb"

# Constructing the amplitude limited case of the variable density spiral Trajectory
function kspaceNodes(tr::SpiralTrajectoryVarDensPeturb)
  nodes = zeros(2,tr.numSamplingPerProfile, tr.numProfiles)
  angles = collect((0:tr.numProfiles-1)/tr.numProfiles) # Phase offset for multiple Profiles
  A = 0.5 # Maximum radius is 0.5
  gamma = 42.58 *1e6 # 42,58 MHz/T : Gyromagnetic relation for protons

  if tr.angleOffset == "golden"
      angles = [i*(3-sqrt(5))/2  for i=0:tr.numProfiles-1 ]
  elseif tr.angleOffset == "random"
      angles = sort(rand(tr.numProfiles))
  elseif tr.angleOffset == "equispaced"
      angles = collect((0:tr.numProfiles-1)/tr.numProfiles)
  end


  T_ea = 1/((tr.alpha+1)*gamma/(2pi*tr.windings))
  times = range(0,stop=T_ea,length=tr.numSamplingPerProfile)
  for l = 1:tr.numProfiles
    for k = 1:tr.numSamplingPerProfile
      tau = (gamma/(2pi*tr.windings) *(tr.alpha+1)*times[k])^(1/(tr.alpha+1))
      tmp = (1+0.025*cos(k*0.1))*(tau^tr.alpha)*(A-0.025)*exp(1im*2*pi*( tr.windings*tau + angles[l] ))
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
  return reshape(nodes, 2, tr.numSamplingPerProfile*tr.numProfiles)
end


function kspaceDensity(tr::SpiralTrajectoryVarDensPeturb)
  error("Not implemented!")
end
