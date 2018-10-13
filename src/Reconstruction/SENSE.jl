# import Base.LinAlg.SingularException
"""
  Recomputes the decomposition for
  every pixel with simplified logic (so far...)
"""
function sense(images
                , sensemaps
                , R::Real)

  Nx,Ny,L = size(sensemaps)
  coords = zeros(Int64,Int64(R),2)
  mask = falses(Nx,Ny)

  assert(L >= R)

  # reduced FOV
  nrx = Nx
  nry = floor(Int64,Ny / R)

  im = zeros(ComplexF64,Nx,Ny)
  s = zeros(ComplexF64,L)

  for ii=1:Nx
    for jj=1:Ny
      if mask[ii,jj] == false

        if abs(sensemaps[ii,jj,1]) < 1e-6
          im[ii,jj] = 0.
        else
          for lx=0:1-1
            for ly=0:R-1
              # Determine pixel positions that superimpose
              ndx = Int64( ((ii-1)+lx*nrx) % Nx +1 )
              ndy = Int64( ((jj-1)+ly*nry) % Ny +1 )

              coords[(ly+1)+R*(lx),1] = ndx
              coords[(ly+1)+R*(lx),2] = ndy

              # Vector of sensitivity at pixel pos. [ndx,ndy]
              # Construct the sensitivity matrix
              ct = vec(sensemaps[ndx,ndy,:])
              if lx == 0 && ly == 0
                s = ct
              elseif abs(ct[1]) >= 0
                s = hcat(s,ct)
              end
            end
          end
          scs = s' *  s
          try
            unfold_mat = inv(scs)*s'
            # Assemble in vector a all the complex img. values of the chosen pixel
            a = vec(images[ii,jj,:])
            mr = unfold_mat  * a

            for j=1:length(mr)
              im[coords[j,1],coords[j,2]] = mr[j]
              mask[coords[j,1],coords[j,2]] = true
            end

          catch ex
            if isa(ex,Base.LinAlg.SingularException)
              println("Position to flame: ", ii, " ", jj)
              println("Matrix to blame: ", s)
              println("Rank of matrix to be inverted: ", rank(scs))
              error("Unfolding failure, inversion is not possible")
            else
              println(ex)
            end
          end
        end
      end
    end
  end
  return im
end
