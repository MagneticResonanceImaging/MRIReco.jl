export sumOfSquareSensitivities, matchedFilterSensitivities, matchedFilter

function sumOfSquareSensitivities(img::Array{ComplexF64})
  nx,ny,nz,nt,nc = size(img)

  img_sum = sqrt.( sum(abs.(img).^2,5) )

  coilsens = zeros(ComplexF64,size(img))
  for i=1:nc
    coilsens[:,:,:,:,i] = img[:,:,:,:,i]./img_sum
  end

  return coilsens
end

function matchedFilterSensitivities(img::Array{ComplexF64,5}, noise::Matrix{ComplexF64})
  nx,ny,nz,nt,nc = size(img)

  # sum over different TE-images
  img_sum = squeeze(sum(img,4),4)

  mf = matchedFilter(img_sum, noise)
  img_ref = zeros(ComplexF64, nx,ny,nz)
  for i = 1:nz
    for j = 1:ny
      for k = 1:nx
        img_ref[k,j,i] = (mf[k,j,i,:]')*img_sum[k,j,i,:]
      end
    end
  end

  coilsens = similar(img_sum)
  for i = 1:nc
    coilsens[:,:,:,i] = img_sum[:,:,i]./img_ref
  end
  return coilsens
end

function matchedFilter(img::Array{ComplexF64,4}, noise::Matrix{ComplexF64})
  # sum over different TE-images
  # img_sum = squeeze(sum(img,4),4)

  # noise correlation matrix
  Rn = (noise*noise')./size(noise,2)

  # coil combination coefficients
  mf = similar(img)
  for i = 1:size(img,3)
    for j = 1:size(img,2)
      for k = 1:size(img,1)
        # signal correlation matrix
        signal = vec(img[k,j,i,:])
        Rs = signal*signal'
        # matched filter
        P = Rn\Rs
        eigFac = eigfact(P)
        maxIdx = findmax( abs.(eigFac[:values]) )[2]
        m = P[:,maxIdx]
        # idx = sub2ind(size(coilSens), (k,j,i,1))
        # coilSens[idx:prod(size(coilSens)[1:3]):end] = conj.( m./sqrt(m'*(Rn\m)) )
        mf[k,j,i,:] = m./sqrt(m'*(Rn\m))
      end
    end
  end
  return mf
end
