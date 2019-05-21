export quadraticFieldmap, stairsFieldmap, polarFieldmap

"""
    quadraticFieldmap(Nx::Int64, Ny::Int64, maxOffresonance::Float64=125.0)

Computes a parabolic fieldmap.
"""
function quadraticFieldmap(Nx::Int64, Ny::Int64, maxOffresonance::Float64=125.0)
  xx = range(-1, stop=1, length=Nx)
  yy = range(-1, stop=1, length=Ny)
  fieldmap = zeros(Nx,Ny, 1)
  for nx=1:Nx
    for ny=1:Ny
      fieldmap[nx,ny,1] = maxOffresonance*(xx[nx]^2 + yy[ny]^2) - maxOffresonance
    end
  end
  return fieldmap
end

"""
Computes a polaroid fieldmap.
"""
function polarFieldmap(Nx::Int64, Ny::Int64,
                        maxOffresonance::Float64=125.0;
                        xCenter::Int64=0,
                        yCenter::Int64=0
                       )
  fieldmap = zeros(Nx,Ny,1)
  xx = range(-1, stop=1, length=Nx)
  yy = range(-1, stop=1, length=Ny)
  for nx=1:Nx
    for ny=1:Ny
      fieldmap[nx,ny,1] = maxOffresonance*sqrt(xx[nx]^2 + yy[ny]^2)- maxOffresonance
    end
  end
  return fieldmap
end


"""
Computes a fieldmap consisting of stairs, which has discrete linear
ascending/descending values
"""
function stairsFieldmap(Nx::Int64, Ny::Int64;
                        numStairs::Int64=8,
                        minOffresonance::Float64=-125.0,
                        maxOffresonance::Float64=125.0,
                        ascending=true)
  # Input check for the number of stairs
  if numStairs < 2 || numStairs > (Nx > Ny ? Ny : Nx)
    error("Number of stairs not used correctly")
  end
  # Allocating output
  fieldmap = zeros(Nx,Ny)
  # Determine the height and width of the stairs (e.g for NxN -> deltaX = deltaY)
  deltaX = Int64(floor((Nx / (numStairs*2) )))
  deltaY = Int64(floor((Ny / (numStairs*2) )))

  if ascending == true
    scaling = range(minOffresonance,stop=maxOffresonance,length=numStairs)
  else
    scaling = -range(minOffresonance,stop=maxOffresonance,length=numStairs)
  end

  # Filling (Quarter) Matrix for constructing the stairs
  partialFieldMap = zeros(Int64(floor((Nx/2))),Int64(floor((Ny/2))))

  for i=1:numStairs
    partialFieldMap[(i-1)*deltaX+1 : end     , (i-1)*deltaY+1 : i*deltaY  ] .= -1 * scaling[i]
    partialFieldMap[(i-1)*deltaX+1 : i*deltaX, (i-1)*deltaY+1 : end ]       .= -1 * scaling[i]
  end

  # Correction if numOfStairs is not divisible by 4 fill the lower right corner of
  # the upper left (quarter) partial matrix
  if numStairs % 4 != 0
    partialFieldMap[numStairs*deltaX:end , numStairs*deltaY:end] .= -1 * scaling[numStairs]
  end

  # If width or height is not divisible by two the middle row or column(or both)
  # must be filled manually
  partSizeX , partSizeY = size(partialFieldMap)
  xRemainder, yRemainder = Nx%2, Ny%2

  # Rotating and mirroring the filled quarter of the matrix to complete the stairs
  fieldmap[1:partSizeX                , 1: partSizeY]                 = partialFieldMap
  fieldmap[1:partSizeX                , partSizeY+1+yRemainder : end] = partialFieldMap[:,end:-1:1]
  fieldmap[partSizeX+1+xRemainder:end , 1: partSizeY]                 = partialFieldMap[end:-1:1,:]
  fieldmap[partSizeX+1+xRemainder:end , partSizeY+1+yRemainder : end] = rot180(partialFieldMap)

  # Filling the row or column,needed when number of rows/columns is odd
  if xRemainder != 0 && yRemainder != 0
    midXIndx, midYIndx  = Int64(ceil(Nx/2)), Int64(ceil(Ny/2))
    fieldmap[midXIndx,: ] = fieldmap[midXIndx-1,:]
    fieldmap[:,midYIndx ] = fieldmap[:,midYIndx-1]
  elseif xRemainder != 0
    midXIndx = Int64(ceil(Nx/2))
    fieldmap[midXIndx,: ] = fieldmap[midXIndx-1,:]
  elseif yRemainder != 0
    midYIndx = Int64(ceil(Ny/2))
    fieldmap[:,midYIndx ] = fieldmap[:,midYIndx-1]
  end

  return reshape(fieldmap,Nx,Ny,1)
end
