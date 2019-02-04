export SpiralTrajectoryVarDens

mutable struct SpiralTrajectoryVarDens <: Abstract2DTrajectory
  numProfiles::Int
  numSamplingPerProfile::Int
  TE::Float64 # echo time in ms
  AQ::Float64 # time for each spiral arm in ms
  alpha::Float64 # densyfing factor
  windings::Float64 # Number of turns the spiral is taking
  angleOffset::String #Symbol
end

# Brought by:  Simple Analytic Variable Density Spiral Design
# by Dong-hyun Kim, 1,2 * Elfar Adalsteinsson, 1 and Daniel M. Spielman
function SpiralTrajectoryVarDens(numProfiles::Int64
                        , numSamplingPerProfile::Int64
                        ; TE::Float64=0.0
                        , AQ::Float64=1.0
                        , windings::Float64=6.25
                        , alpha::Float64=2.0
                        , angleOffset::String= "equispaced" #:equispaced
                        , kargs...)

    SpiralTrajectoryVarDens(numProfiles
                            , numSamplingPerProfile
                            , TE
                            , AQ
                            , alpha
                            , windings
                            , angleOffset
                            )
end

string(tr::SpiralTrajectoryVarDens) = "SpiralVarDens"

# Constructing the amplitude limited case of the variable density spiral Trajectory
function kspaceNodes(tr::SpiralTrajectoryVarDens)
  nodes = zeros(2,tr.numSamplingPerProfile, tr.numProfiles)
  angles = collect((0:tr.numProfiles-1)/tr.numProfiles) # Phase offset for multiple Profiles
  A = 0.5 # Maximum radius is 0.5
  gamma = 42.58 *1e6 # 42,58 MHz/T : Gyromagnetic relation for protons

  # println("angleOffset=$(tr.angleOffset)")

  if tr.angleOffset == "golden" #:golden
      angles = [i*(3-sqrt(5))/2  for i=0:tr.numProfiles-1 ]
  elseif tr.angleOffset == "random" #:random
      # angles = sort(rand(tr.numProfiles))
      angles = collect((0:tr.numProfiles-1)/tr.numProfiles)
      Random.seed!(1234)
      angles = angles .+ 1.0/tr.numProfiles*rand(tr.numProfiles)
  elseif tr.angleOffset == "equispaced" #:equispaced
      angles = collect((0:tr.numProfiles-1)/tr.numProfiles)
  end


  T_ea = 1/((tr.alpha+1)*gamma/(2pi*tr.windings))
  # times = linspace(0,T_ea,tr.numSamplingPerProfile)
  times = range(0,stop=T_ea,length=tr.numSamplingPerProfile)
  for l = 1:tr.numProfiles
    for k = 1:tr.numSamplingPerProfile
      tau = (gamma/(2pi*tr.windings) *(tr.alpha+1)*times[k])^(1/(tr.alpha+1))
      nodes[1,k,l] = A*tau^tr.alpha*cos(2*pi*( tr.windings*tau + angles[l] ))
      nodes[2,k,l] = A*tau^tr.alpha*sin(2*pi*( tr.windings*tau + angles[l] ))
    end
  end
  return reshape(nodes, 2, tr.numSamplingPerProfile*tr.numProfiles)
end


function kspaceDensity(tr::SpiralTrajectoryVarDens)
  error("Not implemented!")
end

isCircular(tr::SpiralTrajectoryVarDens) = true
