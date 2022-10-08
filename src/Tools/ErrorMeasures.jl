export nrmsd, optimalScalingFactor


"""
    optimalScalingFactor(I,Ireco)

computes the optimal scaling factor α such that || I - α Ireco ||₂ is minimal
"""
function optimalScalingFactor(I, Ireco)
  α = norm(Ireco)>0 ? dot(vec(Ireco),vec(I)) / dot(vec(Ireco),vec(Ireco)) : one(eltype(I))
  return α
end

"""
    nrmsd(I,Ireco)

computes the normalized root mean squared error of the image `Ireco`
with respect to the image `I`.
"""
function nrmsd(I, Ireco; optimalScaling=true)
  N = length(I)

  if optimalScaling
    # This is a little trick. We usually are not interested in simple scalings
    # and therefore "calibrate" them away
    α = optimalScalingFactor(I, Ireco)
    Ireco = α*Ireco
  end

  RMS =  1.0/sqrt(N)*norm(vec(I)-vec(Ireco))
  NRMS = RMS/(maximum(abs.(I))-minimum(abs.(I)) )
  return NRMS
end
