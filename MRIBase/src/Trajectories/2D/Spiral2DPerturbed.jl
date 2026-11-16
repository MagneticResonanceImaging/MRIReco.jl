export SpiralPerturbedTrajectory, spiralPerturbedTrajectoryNodes

"""
    SpiralPerturbedTrajectory(T, numProfiles, numSamplingPerProfile
                  ; TE::Float64=0.0
                  , AQ::Float64=1.e-3
                  , windings::Real= 6.25
                  , alpha::Real=2.0
                  , kmax::Real=0.5
                  , angleOffset::String="equispaced"
                  , kargs...)

returns a 2d spiral trajectory with variable density, whose radius at the `k`-th sample of each profile is that of
`SpiralTrajectoryVarDens` with the same arguments times `(1 + 0.025 cos(0.1 k)) (kmax - 0.025) / kmax`. Each
coordinate is limited to the interval [-0.5, 0.5), which changes the nodes only for `kmax` above about 0.51.

# Arguments
* `numProfiles::Int64`            - number of profiles
* `numSamplingPerProfile::Int64`  - number of sampling points per profile
* (`TE::Float64=0.0`)             - echo time in s
* (`AQ::Float64=1.e-3`)           - readout duration in s (per profile)
* (`windings::Real= 6.25`)        - number of windings of the spiral profiles
* (`alpha::Real=2.0`)             - exponent describing the evolution of the magnitude of the sampling points
                                    along the profiles
* (`kmax::Real=0.5`)              - the sampling points lie within this distance from the center of k-space
* (`angleOffset::String="equispaced"`) - spacing of profile angles (`"equispaced"` sampling, `"golden"` angle
                                    sampling or `"random"` sampling)
"""
function SpiralPerturbedTrajectory(::Type{T}, numProfiles, numSamplingPerProfile
                  ; TE=0.0
                  , AQ=1.e-3
                  , windings::Real= 6.25
                  , alpha::Real=2.0
                  , kmax::Real=0.5
                  , angleOffset::String="equispaced"
                  , kargs...) where T
  nodes = spiralPerturbedTrajectoryNodes(T, numProfiles, numSamplingPerProfile;
            windings=windings, alpha=alpha, angleOffset=angleOffset, kmax=kmax)
  times = readoutTimes(T, numProfiles, numSamplingPerProfile; TE=TE, AQ=AQ)
  return  Trajectory("SpiralPerturbed", nodes, times, T(TE), T(AQ), numProfiles, numSamplingPerProfile, 1, false, true)
end

# Constructing the amplitude limited case of the variable density spiral Trajectory
function spiralPerturbedTrajectoryNodes(::Type{T}, numProfiles::Int64
                        , numSamplingPerProfile::Int64
                        ; windings=6.25
                        , alpha=2.0
                        , kmax::Real=0.5
                        , angleOffset= "equispaced"
                        , kargs...) where T

  nodes = zeros(T, 2, numSamplingPerProfile, numProfiles)
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
