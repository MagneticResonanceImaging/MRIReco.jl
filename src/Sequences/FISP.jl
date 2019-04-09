export FISPSequence, echoAmplitudes, flipAngles, numEchoes, TE

#=
  Fast Imaging with Steady State precession sequence with variable flip angles and TR.
  -----------------------
  The phase of the preparation pulse has a phase of -90°
  in order to fulfill CPMG conditions

  echoes apperat at times TE, 2*TE,..., numEchoes*TE after the preparation pulse
=#
mutable struct FISPSequence <: AbstractMultiEchoSequence
  flipAngles :: Vector{Float64}         # flip angles
  TR::Vector{Float64}
  numEchoes :: Int64                    # number of RF pulses
  TE :: Float64                         # echo times relative to the previous RF pulse
  TI :: Float64                         # inversion time
  traj :: Vector{Trajectory}    # readout trajectories
end

function FISPSequence(numEchoes::Int64, TE::Float64, TI::Float64,  N::Int64)
  flipAngles = fill(1.0*pi, numEchoes)
  traj = traj = trajectory("Cartesian",N,N;TE=TE) #SimpleCartesianTrajectory(N,N,TE,1.e-3)
  FISPSequence(flipAngles, numEchoes, TE, traj)
end

function FISPSequence(tr::Trajectory; numEchoes=1, TE = 1.e-2, TI = 7.e-3, flipAngles=nothing, TR=nothing, kargs...)
  alpha = zeros(numEchoes)
  if flipAngles==nothing
    alpha = fill(pi/2., numEchoes)
  else
    alpha=flipAngles
  end

  repTimes = zeros(numEchoes)
  if TR==nothing
    repTimes = fill(1.5e-2, numEchoes)
  else
    repTimes = TR
  end

  #  traj = fill(tr, numEchoes)
  trajectories = [deepcopy(tr) for i = 1:numEchoes]
  for i = 1:numEchoes
    trNew=deepcopy(tr)
    trNew.TE = i*TE+repTimes[i]
    push!(trajectories, trNew)
  end

  return FISPSequence(alpha, repTimes, numEchoes, TE, TI, trajectories)
end

trajectory(seq::FISPSequence, n::Int64=1) = seq.traj[n]

numEchoes(seq::FISPSequence)  = seq.numEchoes

# FIXME: Provide correct implementation
echoTimes(seq::FISPSequence) = [ i*seq.traj.TE for i=1:numEchoes(seq) ]

flipAngles(seq::FISPSequence) = seq.flipAngles

string(seq::FISPSequence) = "FISP"

 #=
   calculate echo amplitudes using the extended phase graph method.
   echo times are given in units of TR
 =#
function echoAmplitudes(seq::FISPSequence, R1::Float64, R2::Float64, numStates=nothing)

  # number of states to consider
  numStates!=nothing ? ns=numStates : ns=seq.numEchoes

  # calculate amplitudes of longitudinal and dephased states
  f, z = epgAmplitudes(seq, R1, R2, ns)

  # evolution of the relevant states until TE
  amp = f[ns,1:end-1].*exp(-seq.TE*R2)

  return amp
end

#=
 calculate amplitudes of dephased and longitudinal states of the magnetization
 using the extended phase graph method.
=#
function epgAmplitudes(seq::FISPSequence, R1::Real, R2::Real, numStates=nothing)
  # repetion time and number of pulses
  np=seq.numEchoes

  # initial states for the dephased (f) and longitudinal (z) magnetization.
  # the magnetization after the RF pulse is described by f_deph and z_deph.
  numStates!=nothing ? ns=numStates : ns=np
  f = zeros(ComplexF64,2*ns+1,np+1)
  z = zeros(ComplexF64,ns,np+1)
  # f = zeros(ComplexF64,2*np+1,np+1)
  # z = zeros(ComplexF64,np,np+1)

  # an initial 90 degree pulse fully populates F0
  z[1,1] = -1.

  fDeph = epgDephasing(f[:,1],1)
  # fDeph,zDeph = epgRelaxation(R1, R2, seq.TR[1], f[:,1], z[:,1])
  fDeph,zDeph = epgRelaxation(R1, R2, seq.TI, f[:,1], z[:,1])

  # loop over inversion-pulses
  for i = 1:np

    f[:,i], z[:,i] = epgRotation(seq.flipAngles[i], fDeph, zDeph, statesConsidered=min(i,ns))

    fDeph = epgDephasing(f[:,i],1)
    fDeph,zDeph = epgRelaxation(R1, R2, seq.TR[i], fDeph, z[:,i])
  end

  return f, z
 end

 function encoding(seq::FISPSequence)
   tr_type = typeof(trajectory(seq,1))
   if tr_type <: Abstract2DTrajectory
     return "2D"
   elseif tr_type <: Abstract3DTrajectory
     return "3D"
   end
   return "other"
 end
