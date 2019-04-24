export MESequence, echoAmplitudes, flipAngles, numEchoes, T_echo

#=
  General Multi-Echo sequence with variable flip angles and TR.
  -----------------------
  The phase of the preparation pulse has a phase of -90°
  in order to fulfill CPMG conditions

  echoes apperat at times T_echo, 2*T_echo,..., numEchoes*T_echo after the preparation pulse
=#
mutable struct MESequence <: AbstractSequence
  excitationAngle::Float64
  refocusingAngles :: Vector{Float64}   # flip angles
  T_rf::Vector{Float64}                 # times of refocusing pulses
  T_echo:: Vector{Float64}                 # echo times
end

function MESequence(; T_echo::Union{Float64,Vector{Float64}}=[0.0]
                    , T_rf::Union{Float64,Vector{Float64}}=Float64[]
                    , excitationAngle::Float64=pi/2.0
                    , refocusingAngles::Union{Float64,Vector{Float64}}=fill(Float64(pi),length(T_echo))
                    , kargs...)
  # turn arguments to vectors if necessary
  refAng_vec = ( typeof(refocusingAngles)==Float64 ? Float64[refocusingAngles] : refocusingAngles )
  T_rf_vec = ( typeof(T_rf)==Float64 ? Float64[T_rf] : T_rf )
  T_echo_vec = ( typeof(T_echo)==Float64 ? Float64[T_echo] : T_echo )

  # calculate rf-timings if none are given
  if isempty(T_rf_vec)
    T_rf_vec = zeros(length(T_echo_vec))
    T_rf_vec[1] = T_echo_vec[1]/2.0
    T_rf_vec[2:end] .= [0.5*(T_echo_vec[i]+T_echo_vec[i-1]) for i=2:length(T_echo_vec)]
  end

  return MESequence(excitationAngle,refAng_vec, T_rf_vec, T_echo_vec)
end

numEchoes(seq::MESequence)  = length(seq.T_echo)

echoTimes(seq::MESequence) = seq.T_echo

flipAngles(seq::MESequence) = seq.refocusingAngles

string(seq::MESequence) = "ME"

 #=
   calculate echo amplitudes using the extended phase graph method.
   echo times are given in units of TR
 =#
function echoAmplitudes(seq::MESequence, R1::Float64, R2::Float64, numStates=nothing)
  # calculate amplitudes of longitudinal and dephased states
  f, z = epgAmplitudes(seq, R1, R2)

  # evolution of the relevant states until T_echo
  T_rel = seq.T_echo .- seq.T_rf
  amp = f[2*numEchoes(seq)-1,2:end] .* exp.(-T_rel*R2)

  return amp
end

#=
 calculate amplitudes of dephased and longitudinal states of the magnetization
 using the extended phase graph method.
=#
function epgAmplitudes(seq::MESequence, R1::Real, R2::Real, numStates=nothing)
  # repetion time and number of pulses
  np=numEchoes(seq)

  # initial states for the dephased (f) and longitudinal (z) magnetization.
  # the magnetization after the RF pulse is described by f_deph and z_deph.
  f = zeros(ComplexF64,4*np-1,np+1)
  z = zeros(ComplexF64,2*np,np+1)

  # an initial 90 degree pulse fully populates F0
  z[1,1] = 1.

  fDeph, zDeph = epgRotation(seq.excitationAngle ,f[:,1],z[:,1],statesConsidered=1, phi=pi/2)

  fDeph = epgDephasing(fDeph,1)
  fDeph,zDeph = epgRelaxation(R1, R2, seq.T_rf[1], fDeph, zDeph)

  # loop over refocusing-pulses
  for i = 2:np
    f[:,i], z[:,i] = epgRotation(seq.refocusingAngles[i-1], fDeph, zDeph, statesConsidered=2*(i-1))

    fDeph = epgDephasing(f[:,i],2)
    fDeph,zDeph = epgRelaxation(R1, R2, seq.T_rf[i]-seq.T_rf[i-1], fDeph, z[:,i])
  end

  # last pulse without (omit dephasing after pulse)
  f[:,np+1], z[:,np+1] = epgRotation(seq.refocusingAngles[np], fDeph, zDeph, statesConsidered=2*np)

  return f, z
 end
