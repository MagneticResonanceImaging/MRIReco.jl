# import Base.LinAlg.SingularException

function gfactor(sensemaps,rx,ry)

  Nx,Ny,L = size(sensemaps)
  assert(L >= rx*ry )

  # reduced FOV
  nrx = floor(Int64,Nx / rx)
  nry = floor(Int64,Ny / ry)

  s = zeros(ComplexF64,L)
  nc = zeros(ComplexF64,Nx,Ny)
  g = zeros(ComplexF64,Nx,Ny)
  for ii=1:Nx
    for jj=1:Ny
      if abs(sensemaps[ii,jj,1]) < 1e-6
        g[ii,jj] = 0
      else
        for lx=0:rx-1
          for ly=0:ry-1
            # Determine pixel positions that superimpose
            ndx = Int64( ((ii-1)+lx*nrx) % Nx +1 )
            ndy = Int64( ((jj-1)+ly*nry) % Ny +1 )

            # Vector of sensitivity at pixel pos. [ndx,ndy]
            # Construct the sensitivity matrix
            ct = vec(sensemaps[ndx,ndy,:])
            if lx == 0 && ly == 0
              s = ct
            elseif abs(ct[1]) >= 0
              s = hcat(s,ct)
            end
            nc[ii,jj] = nc[ii,jj]+1;
          end
        end
        scs = s' * s
        # Unfolding matrix
        try
          s_pinv = inv(s'*s) * s'
          g[ii,jj] = sqrt(complex(scs[1,1]*s_pinv[1,1]))
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
  return g,nc

end
