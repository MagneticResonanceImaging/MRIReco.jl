export scalePhantom

function scalePhantom(A::Matrix{T},factor::Float64,interpolation="inverse_nearestneighbour") where T<:Real
imageSize = size(A);
newImage = zeros(eltype(A),imageSize);
#newImageSize = size(newImage);
newImageSize=imageSize

newCornerI = round(Int64,imageSize[1]*(1-factor)/2)
newCornerJ = round(Int64,imageSize[2]*(1-factor)/2)

if interpolation == "direct"
    for i=1:imageSize[1]
        for j=1:imageSize[2]
            newImage[round(Int64,i*factor),round(Int64,j*factor)] = A[i,j];
        end
    end

elseif interpolation == "inverse_nearestneighbour"
    for i=1:newImageSize[1]
        for j=1:newImageSize[2]
            oldI = round(Int64,(i-newCornerI)/factor);
            oldJ = round(Int64,(j-newCornerJ)/factor);

            if oldI < 1
                newImage[i,j] = 0
            elseif oldJ < 1
                newImage[i,j] = 0
            elseif oldI > imageSize[1]
              newImage[i,j] = 0
            elseif oldJ > imageSize[2]
              newImage[i,j] = 0
            else
              newImage[i,j] = A[oldI,oldJ];
            end
        end
    end

elseif interpolation == "inverse_bilinear"
    for i=1:newImageSize[1]
        for j=1:newImageSize[2]
            bottomX = floor(Int64,(j-newCornerJ)/factor);
            bottomY = floor(Int64,(i-newCornerI)/factor);
            topX    = ceil(Int64,(j-newCornerJ)/factor);
            topY    = ceil(Int64,(i-newCornerI)/factor);

            if topX < 1
              newImage[i,j] = 0
              continue
            end

            if topY < 1
              newImage[i,j] = 0
              continue
            end

            if bottomX > imageSize[2]
              newImage[i,j] = 0
              continue
            end

            if bottomY > imageSize[2]
              newImage[i,j] = 0
              continue
            end

            if bottomX < 1
                bottomX = 1;
            end
            if bottomY < 1
                bottomY = 1;
            end
            if topX > imageSize[2]
                topX = imageSize[2];
            end
            if topY > imageSize[1]
                topY = imageSize[1];
            end

            alpha = (j-newCornerJ)/factor - bottomX;
            beta  = (i-newCornerI)/factor - bottomY;

            Vax = bottomX;
            Vay = bottomY;

            Vbx = topX;
            Vby = bottomY;

            Vcx = bottomX;
            Vcy = topY;

            Vdx = topX;
            Vdy = topY;

            newImage[i,j] = ((1-alpha)*(1-beta)*A[Vay,Vax]
                            +alpha*(1-beta)*A[Vby,Vbx]
                            +(1-alpha)*beta*A[Vcy,Vcx]
                            +alpha*beta*A[Vdy,Vdx])
        end
    end

else
    newImage = A;
end

return newImage

end

function scalePhantom(A::Array{T,3},factor::Float64,interpolation="inverse_nearestneighbour") where T
imageSize = size(A);
newImage = zeros(eltype(A),imageSize);
#newImageSize = size(newImage);
newImageSize=imageSize

newCornerI = round(Int64,imageSize[1]*(1-factor)/2)
newCornerJ = round(Int64,imageSize[2]*(1-factor)/2)
newCornerK = round(Int64,imageSize[3]*(1-factor)/2)

if interpolation == "direct"
  for k=1:imageSize[3]
    for j=1:imageSize[2]
        for i=1:imageSize[1]
            newImage[round(Int64,i*factor),round(Int64,j*factor),round(Int64,k*factor)] = A[i,j,k];
        end
    end
  end

elseif interpolation == "inverse_nearestneighbour"
    for k=1:imageSize[3]
      for j=1:newImageSize[2]
        for i=1:newImageSize[1]

            oldI = round(Int64,(i-newCornerI)/factor);
            oldJ = round(Int64,(j-newCornerJ)/factor);
            oldK = round(Int64,(k-newCornerK)/factor);

            if oldI < 1
                newImage[i,j,k] = 0
            elseif oldJ < 1
                newImage[i,j,k] = 0
            elseif oldK < 1
                newImage[i,j,k] = 0
            elseif oldI > imageSize[1]
              newImage[i,j,k] = 0
            elseif oldJ > imageSize[2]
              newImage[i,j,k] = 0
            elseif oldK >imageSize[3]
              newImage[i,j,k] = 0
            else
              newImage[i,j,k] = A[oldI,oldJ,oldK];
            end
        end
      end
    end

else
    newImage = A;
end

return newImage

end
