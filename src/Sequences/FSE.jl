export FSESequence, echoAmplitudes, flipAngles, numEchoes, TE

#=
  Fast Spin Echo Sequence with a variable flip angle scheme.
  -----------------------
  The phase of the preparation pulse has a phase of -90°
  in order to fulfill CPMG conditions

  echoes apperat at times TE, 2*TE,..., numEchoes*TE after the preparation pulse
=#
mutable struct FSESequence <: AbstractMultiEchoSequence
  flipAngles :: Vector{Float64}         # flip angles
  numEchoes :: Int64                    # number of RF pulses
  TE :: Float64                         # echo times relative to the previous RF pulse
  traj :: Vector{AbstractTrajectory}    # readout trajectories
end

function FSESequence(numEchoes::Int64, TE::Float64, N::Int64)
  flipAngles = fill(1.0*pi, numEchoes)
  traj = SimpleCartesianTrajectory(N,N,TE,1.e-3)
  FSESequence(flipAngles, numEchoes, TE, traj)
 end

 function FSESequence(tr::AbstractTrajectory; numEchoes=1, TE = 1.e-3, flipAngles=nothing, kargs...)
   alpha = zeros(numEchoes)
   if flipAngles==nothing
     alpha = fill(180, numEchoes)
   else
     alpha=flipAngles
   end

   traj = fill(tr, numEchoes)

   return FSESequence(alpha, numEchoes, TE, traj)
 end

 #=
   calculate echo amplitudes using the extended phase graph method.
   echo times are given in units of TR
 =#
 function echoAmplitudes(seq::FSESequence, R1::Float64, R2::Float64)

   # calculate amplitudes of longitudinal and dephased states
   f, z = epgAmplitudes(seq, R1, R2)

   # evolution of the relevant states until TE --> TODO: revise
   amp = f[2*seq.numEchoes-1,2:end].*exp(-0.5*seq.TE*R2)

   return amp
 end

 #=
   calculate amplitudes of dephased and longitudinal states of the magnetization
   using the extended phase graph method.
 =#
 function epgAmplitudes(pulse::FSESequence, R1::Real, R2::Real)
   # repetion time and number of pulses
   TR = pulse.TE
   np=pulse.numEchoes

   # initial states for the dephased (f) and longitudinal (z) magnetization.
   # the magnetization after the RF pulse is described by f_deph and z_deph.
   f = zeros(ComplexF64,4*np-1,np+1)
   z = zeros(ComplexF64,2*np,np+1)

   # an initial 90 degree pulse fully populates F0
   f[2*np,1] = 1

   fDeph = epgDephasing(f[:,1],1)
   fDeph,zDeph = epgRelaxation(R1, R2, 0.5*TR, fDeph, z[:,1])

   # loop over inversion-pulses
   for i = 2:np+1

     f[:,i], z[:,i] = epgRotation(pulse.flipAngles[i-1], fDeph, zDeph, statesConsidered=2*(i-1))

     fDeph = epgDephasing(f[:,i],2)
     fDeph,zDeph = epgRelaxation(R1, R2, TR, fDeph, z[:,i])
   end

   return f, z
 end


flipAngles(seq::FSESequence) = seq.flipAngles

numEchoes(seq::FSESequence)  = seq.numEchoes

TE(seq::FSESequence) = seq.TE

string(seq::FSESequence) = "FSE"
