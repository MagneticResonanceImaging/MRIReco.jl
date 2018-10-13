export epgRotation, epgRelaxation, epgDephasing, rfRotation

#
# apply Bloch-rotation (<=> RF pulse) to a set of EPG states
#
function epgRotation(alpha::Float64, F::Vector{T}, Z::Vector{T}; statesConsidered=nothing) where T

  # apply rotation to all states per default
  numStates = length(Z)

  n = numStates
  statesConsidered!=nothing && (n=statesConsidered)

  R = rfRotation(alpha, 0.)

  rotatedF = fill!(copy(F), 0)
  rotatedZ = fill!(copy(Z), 0)

  for i = 0:n-1
    posIdx = numStates + i
    negIdx = numStates - i

    rotStates = R*[F[posIdx]; conj(F[negIdx]); Z[i+1]]
    rotatedF[posIdx] = rotStates[1]
    rotatedF[negIdx] = conj(rotStates[2])
    rotatedZ[i+1] = rotStates[3]
  end

  return rotatedF, rotatedZ
end

#
# calculate relaxation between two RF pulses for a set of EPG states
#
function epgRelaxation( R1::Float64, R2::Float64, t::Float64, F::Vector{T}, Z::Vector{T}) where T

  relaxedF = exp(-R2*t).*F
  relaxedZ = fill!(copy(Z), 0)
  relaxedZ[2:end] = exp(-R1*t)*Z[2:end]
  relaxedZ[1] = exp(-R1*t)*(Z[1]-1.0) + 1.0

  return relaxedF, relaxedZ
end

#
# calculate dephasing of EPG-states
#
epgDephasing(F::Vector{T}, n=1) where T = circshift(F[:],n)

#
# RF-rotation matrix
#
function rfRotation(alpha, phi=0.)
  R = [ cos(alpha/2.)^2   exp(2*im*phi)*sin(alpha/2.)^2   -im*exp(im*phi)*sin(alpha);
        exp(-2*im*phi)*sin(alpha/2.)^2   cos(alpha/2.)^2   im*exp(-im*phi)*sin(alpha);
        -im/2 .*exp(-im*phi)*sin(alpha)   im/2 .*exp(im*phi)*sin(alpha)   cos(alpha) ]
end
