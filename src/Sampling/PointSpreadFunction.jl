
"""
  Alternative calculation of incoherence using
"""
function pointspreadfunctionAlt(A::AbstractLinearOperator,shape::Tuple)

  tres = zeros(ComplexF64,prod(shape))
  unitvec_i = zeros(ComplexF64,prod(shape))
  unitvec_i[floor(Int64,prod(shape)/2 - shape[1]/2) ] = 1.

  println(floor(Int64,prod(shape)/2))
  unitvec_j = zeros(ComplexF64,prod(shape))

    B = zeros(ComplexF64,prod(shape), prod(shape) )
  for j=1:prod(shape)
    unitvec_j[j] = 1.
    uno = A * unitvec_j
    B[:,j] = uno
  end

  C = B'*B
  D = 1 ./diag(C)
  C = diagm(D) * B
  return C

end

function pointspreadfunction(A::AbstractLinearOperator
                            , shape::Tuple
                            ; verbose::Bool=true)
    numOfPixel = prod(shape)
    unitvec_i = zeros(ComplexF64,numOfPixel)
    unitvec_j = zeros(ComplexF64,numOfPixel)
    psfunc = zeros(ComplexF64,numOfPixel,numOfPixel)

    if verbose == true
        progr = Progress(numOfPixel^2+1,dt=0.1,desc="Calc PSF...";barglyphs=BarGlyphs("[=> ]"),barlen=50);
    end

    for i=1:numOfPixel
        unitvec_i = zeros(ComplexF64,numOfPixel)
        unitvec_i[i] = 1
        tmp = (vec(A\ (A.density .* (A * unitvec_i))) )
        for j=1:numOfPixel
            unitvec_j = zeros(ComplexF64,numOfPixel)
            unitvec_j[j] = 1
            psfunc[i,j] = dot(unitvec_j,tmp)
            if verbose == true
                next!(progr)
            end
        end
    end
    peak = maximum(abs(psfunc))
    println(peak)
    for i=1:numOfPixel
        psfunc[i,i] = 0.0
    end
    sidelobe = maximum(abs(psfunc))
    println("sidelobe to peak ration ",sidelobe/peak)
    return psfunc
end


"""
  Direct calculation of the point spread function
"""
function pointspreadfunction2(A::AbstractLinearOperator,shape::Tuple)
  numOfPixel = prod(shape)
  tres = zeros(ComplexF64,prod(shape))
  unitvec_i = zeros(ComplexF64,prod(shape))
  unitvec_i[floor(Int64,prod(shape)/2 - shape[1]/2) ] = 1.

  println(floor(Int64,prod(shape)/2))
  unitvec_j = zeros(ComplexF64,prod(shape))

  for j=1:prod(shape)
    unitvec_j[j] = 1.
    uno = A * unitvec_j
    unitvec_j[j] = 0.
    unouno = A\(uno)
    tres[j] = dot(unitvec_i , vec(unouno))
    unitvec_j[j] = 0.
  end

  backup = tres[floor(Int64,prod(shape)/2 - shape[1]/2)]
  tres[floor(Int64,prod(shape)/2 - shape[1]/2)] = 0.
  println("sidelobe-to-peak ratio:",maximum(abs(tres/backup)))
  tres[floor(Int64,prod(shape)/2 - shape[1]/2) ] = backup
  return reshape(tres,shape)
end
