export nrmsd

function nrmsd(I,Ireco)
  N = length(I)

  # This is a little trick. We usually are not interested in simple scalings
  # and therefore "calibrate" them away
  alpha = (dot(vec(I),vec(Ireco))+dot(vec(Ireco),vec(I))) /
          (2*dot(vec(Ireco),vec(Ireco)))
  Ireco[:] .*= alpha

  RMS =  1.0/sqrt(N)*norm(vec(I)-vec(Ireco))
  NRMS = RMS/(maximum(abs.(I))-minimum(abs.(I)) )
  return NRMS
end
