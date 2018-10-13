export rotatePhantom

function rotatePhantom(A::Matrix{T},angle::Float64,interpolate="inverse_nearestneighbour") where T<:Real

    imageSize = size(A);
    newImageSize = zeros(Int64,2)
    #newImageSize[1] = max(round(Int64,abs(imageSize[1]*sin(angle)+imageSize[2]*cos(angle))),size(A,1));
    #newImageSize[2] = max(round(Int64,abs(imageSize[2]*sin(angle)+imageSize[1]*cos(angle))),size(A,2));
    #newImageSize[1] = round(Int64,abs(imageSize[1]*sin(angle)+imageSize[2]*cos(angle)));
    #newImageSize[2] = round(Int64,abs(imageSize[2]*sin(angle)+imageSize[1]*cos(angle)));
    newImageSize = imageSize

    newImage = zeros(eltype(A),newImageSize[1],newImageSize[2]);

if interpolate == "inverse_nearestneighbour"
    for i=1:newImageSize[1]
        for j=1:newImageSize[2]
            oldI = (i-newImageSize[1]/2)*cos(angle)+(j-newImageSize[2]/2)*sin(angle);
            oldJ = (j-newImageSize[1]/2)*cos(angle)-(i-newImageSize[2]/2)*sin(angle);

            oldI = round(Int64,oldI+imageSize[1]/2);
            oldJ = round(Int64,oldJ+imageSize[2]/2);

            if oldI < 1
                continue;
            end
            if oldJ < 1
                continue;
            end
            if oldI > imageSize[1]
                continue;
            end
            if oldJ > imageSize[2]
                continue;
            end

            newImage[i,j] = A[oldI,oldJ];

        end
    end
elseif interpolate == "inverse_bilinear"

    for i=1:newImageSize[1]
        for j=1:newImageSize[2]
            oldX = (j-newImageSize[1]/2)*cos(angle)-(i-newImageSize[2]/2)*sin(angle)+imageSize[1]/2;
            oldY = (i-newImageSize[1]/2)*cos(angle)+(j-newImageSize[2]/2)*sin(angle)+imageSize[2]/2;

            bottomX = floor(Int64,oldX);
            bottomY = floor(Int64,oldY);
            topX    = ceil(Int64,oldX);
            topY    = ceil(Int64,oldY);

            if bottomX < 1
                continue;
            end
            if bottomY < 1
                continue;
            end

            if topX < 1
                continue;
            end
            if topY < 1
                continue;
            end

            if bottomX > imageSize[1]
                continue;
            end
            if bottomY > imageSize[2]
                continue;
            end

            if topX > imageSize[2]
                continue;
            end
            if topY > imageSize[1]
                continue;
            end

            alpha = oldX - bottomX;
            beta  = oldY - bottomY;

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

end

  return newImage;

end

function rotatePhantom(A::Array{T,3},angleX::Float64,angleY::Float64,angleZ::Float64,interpolate="inverse_nearestneighbour") where T<:Real
  imageSize = size(A);
  newImageSize = imageSize
  newImage = zeros(eltype(A),newImageSize[1],newImageSize[2],newImageSize[3]);
  newImage[:] = A[:]

  for i=1:size(A,1)
    newImage[i,:,:]= rotatePhantom(newImage[i,:,:],angleX,interpolate)
  end

  for i=1:size(A,2)
    newImage[:,i,:]= rotatePhantom(newImage[:,i,:],angleY,interpolate)
  end

  for i=1:size(A,3)
    newImage[:,:,i]= rotatePhantom(newImage[:,:,i],angleZ,interpolate)
  end

  return newImage

end
