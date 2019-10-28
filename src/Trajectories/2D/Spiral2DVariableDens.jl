export SpiralTrajectoryVarDens, spiralVarDensNodes, spiralVarDensDensity

"""
    SpiralTrajectoryVarDens(numProfiles, numSamplingPerProfile
                  ; TE::Float64=0.0
                  , AQ::Float64=1.e-3
                  , windings::Real= 6.25
                  , alpha=2.0
                  , angleOffset= :equispaced
                  , kargs...)

returns a 2d spiral trajectory with variable density

# Arguments
* `numProfiles::Int64`            - number of profiles
* `numSamplingPerProfile::Int64`  - number of sampling points per profile
* (`TE::Float64=0.0`)             - echo time in s
* (`AQ::Float64=1.e-3`)           - readout duration in s (per profile)
* (`windings::Real= 6.25`)        - number of windings of the spiral profiles
* (`alpha=2.0`)                   - exponent describing the evolution of the magnitude of the sampling points along the profiles
* (`angleOffset= :equispaced`)    - spacing of profile angles (`:equispaced` sampling, `:golden` angle sampling or `:random` sampling)
"""
function SpiralTrajectoryVarDens(numProfiles, numSamplingPerProfile
                  ; TE::Real=0.0
                  , AQ::Real=1.e-3
                  , windings::Real= 6.25
                  , alpha::Real=2.0
                  , kmax::Real=0.5
                  , angleOffset::String="equispaced"
                  , kargs...)
  nodes = spiralVarDensNodes(numProfiles, numSamplingPerProfile; windings=windings,
                             alpha=alpha, angleOffset=angleOffset, kmax=kmax)
  times = readoutTimes(numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ)
  return  Trajectory("SpiralVarDens", nodes, times, TE, AQ, numProfiles, numSamplingPerProfile, 1, false, true)
end

# Constructing the amplitude limited case of the variable density spiral Trajectory
function spiralVarDensNodes(numProfiles::Int64
                        , numSamplingPerProfile::Int64
                        ; windings::Real=6.25
                        , alpha::Real=2.0
                        , kmax::Real=0.5
                        , angleOffset::String= "equispaced" #:equispaced
                        , kargs...)

  nodes = zeros(2,numSamplingPerProfile, numProfiles)
  
  if angleOffset == "golden" #:golden
      angles = [i*(3-sqrt(5))/2  for i=0:numProfiles-1 ]
  elseif angleOffset == "random" #:random
      # angles = sort(rand(numProfiles))
      angles = collect((0:numProfiles-1)/numProfiles)
      Random.seed!(1234)
      angles = angles .+ 1.0/numProfiles*rand(numProfiles)
  elseif angleOffset == "equispaced" #:equispaced
      angles = collect((0:numProfiles-1)/numProfiles)
  end

  times = range(0,stop=1,length=numSamplingPerProfile)
  for l = 1:numProfiles
    for k = 1:numSamplingPerProfile
      tau = (times[k])^(1/(alpha+1))
      nodes[1,k,l] = kmax*tau^alpha*cos(2*pi*( windings*tau + angles[l] ))
      nodes[2,k,l] = kmax*tau^alpha*sin(2*pi*( windings*tau + angles[l] ))
    end
  end
  return reshape(nodes, 2, numSamplingPerProfile*numProfiles)
end


function spiralVarDensDensity(numProfiles::Int64, numSamplingPerProfile::Int64)
  error("Not implemented!")
end
