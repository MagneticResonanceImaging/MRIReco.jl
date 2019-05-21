export epgRotation, epgRelaxation, epgDephasing, rfRotation

"""
    epgRotation(alpha::Float64, F::Vector{T}, Z::Vector{T}; statesConsidered=nothing, phi::Float64=0.0)

applies Bloch-rotation (<=> RF pulse) to a set of EPG states.

# Arguments
* `alpha::Float64`           - flip angle of the RF pulse
* `F::Vector{T}`             - transverse dephasing stats
* `Z::Vector{T}`             - longitudinal dephasing stats
* `statesConsidered=nothing` - number of dephasing states to consider (nothing means all states are taken into account)
* `phi::Float64=0.0`         - phase of the RF pulse
"""
function epgRotation(alpha::Float64, F::Vector{T}, Z::Vector{T}; statesConsidered=nothing, phi::Float64=0.0) where T
  # apply rotation to all states per default
  numStates = length(Z)

  n = numStates
  statesConsidered!=nothing && (n=statesConsidered)

  R = rfRotation(alpha, phi)

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

"""
    epgRelaxation( R1::Float64, R2::Float64, t::Float64, F::Vector{T}, Z::Vector{T}) where T

applies relaxation matrices to a set of EPG states.

# Arguments
* `R1::Float64`   - R1
* `R2::Float64`   - R2
* `t::Float64`    - length of time interval in s
* `F::Vector{T}`  - transverse dephasing stats
* `Z::Vector{T}`  - longitudinal dephasing stats
"""
function epgRelaxation( R1::Float64, R2::Float64, t::Float64, F::Vector{T}, Z::Vector{T}) where T

  relaxedF = exp(-R2*t).*F
  relaxedZ = fill!(copy(Z), 0)
  relaxedZ[2:end] = exp(-R1*t)*Z[2:end]
  relaxedZ[1] = exp(-R1*t)*(Z[1]-1.0) + 1.0

  return relaxedF, relaxedZ
end

"""
    epgDephasing(F::Vector{T}, n=1) where T = circshift(F[:],n)

shifts the transverse dephasing states `F` corresponding to n dephasing-cycles.
"""
epgDephasing(F::Vector{T}, n=1) where T = circshift(F[:],n)

"""
    rfRotation(alpha, phi=0.)

returns the rotation matrix for a pulse with flip angle `alpha` and phase `phi`.
"""
function rfRotation(alpha, phi=0.)
  R = [ cos(alpha/2.)^2   exp(2*im*phi)*sin(alpha/2.)^2   -im*exp(im*phi)*sin(alpha);
        exp(-2*im*phi)*sin(alpha/2.)^2   cos(alpha/2.)^2   im*exp(-im*phi)*sin(alpha);
        -im/2 .*exp(-im*phi)*sin(alpha)   im/2 .*exp(im*phi)*sin(alpha)   cos(alpha) ]
end
