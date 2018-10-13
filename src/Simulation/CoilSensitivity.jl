export birdcageSensitivity
"""

Computes the sensitivity maps for each coils that arranged
in a birdcage manner.

...
# Arguments
* `matrix_size::Int64` : quadratic size of a sensitivity map
* `ncoils::Int64` : the number of coils used
* `relative_radius::Float64` : the relative radius of the birdcage
...
"""
function birdcageSensitivity(matrix_size::Int64
                            , ncoils::Int64
                            , relative_radius::Float64
                            )

  out = zeros(ComplexF64,matrix_size,matrix_size,ncoils)
  for c=0:ncoils-1
    coilx = relative_radius*cos(c*(2*pi/ncoils))
    coily = relative_radius*sin(c*(2*pi/ncoils))
    coil_phase = -c*(2*pi/ncoils)

    for y=0:matrix_size-1
      y_co = (y - (matrix_size/2))/(matrix_size/2) - coily
      for x=0:matrix_size-1
        x_co = (x - (matrix_size/2))/(matrix_size/2) - coilx
        rr = sqrt(x_co^2+y_co^2)
        phi = atan(x_co, -y_co) + coil_phase
        out[x+1,y+1,c+1] = 1/(rr) * exp(1im*phi)
      end
    end
  end
 return out
end
