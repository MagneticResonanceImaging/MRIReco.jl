export birdcageSensitivity

"""
    birdcageSensitivity(N::Int64, ncoils::Int64, relative_radius::Float64)

Computes the sensitivity maps for each coils that are arranged
in a birdcage manner.
"""
function birdcageSensitivity(N::Int64,
                             ncoils::Int64,
                             relative_radius::Float64)

  out = zeros(ComplexF64, N, N, 1, ncoils)
  for c=0:ncoils-1
    coilx = relative_radius*cos(c*(2*pi/ncoils))
    coily = relative_radius*sin(c*(2*pi/ncoils))
    coil_phase = -c*(2*pi/ncoils)

    for y=0:N-1
      y_co = (y - (N/2))/(N/2) - coily
      for x=0:N-1
        x_co = (x - (N/2))/(N/2) - coilx
        rr = sqrt(x_co^2+y_co^2)
        phi = atan(x_co, -y_co) + coil_phase
        out[x+1,y+1, 1, c+1] = 1/(rr) * exp(1im*phi)
      end
    end
  end
 return out
end
