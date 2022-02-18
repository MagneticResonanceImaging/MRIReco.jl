
"""
  Returns the a_j,k and c_k,p Coefficients considering the approximation of the
  exponential term exp(-t_j * z_p) = sum( k : a_j,k * c_k,p )
"""
function get_AC_coefficients(K::Int64
                            , times::Vector
                            , z_p::Vector{Complex{T}}
                            , alpha::Float64
                            , m::Float64
                            ; method = "nfft") where T

  if method == "leastsquare"
    A,C = get_AC_coefficients_lsqr(K,times,z_p)
  elseif method == "nfft"
    A,C = get_AC_coefficients_nfft(K,times,z_p,alpha,m)
  elseif method == "hist"
    A,C = get_AC_coefficients_hist_lsqr(K,times,z_p)
  else
    A = zeros(Complex{T},length(times),K)
    C = zeros(Complex{T},K,length(z_p))
    error("Not implemented yet")
  end
  return A,C
end


function _precompute_coeffs_nfft(K::Int64
                                , times::Vector     # readout times
                                , z_p::Vector{Complex{T}}      # Correction map
                                , alpha::Float64   # oversampling factor of NFFT
                                , m::Float64       # kernel size of NFFT
                                ; adaptK::Bool=true
                                , centering::Bool=true) where T

  # Centering time and correctionterm to lower computational effort
  t_hat = (times[1] + times[end])/2
  z_hat = minimum(real(z_p)) + maximum(real(z_p)) + 1im*(minimum(imag(z_p)) + maximum(imag(z_p)))
  z_hat *= 0.5
  if !centering
    t_hat=0.0
    z_hat=0.0
  end
  # Centering readouttimes and correctionmap
  times = times .- t_hat
  z_p = z_p .- z_hat

  # Calculating timescaling factor so it satisfies: t_j / T in [0.5,0.5)
  t_scale = maximum(abs.(times)) / (0.5-1e-7)

  # determining N and K,
  # K can be choosen freely with tradeoff comp. complexity vs precision
  N = Int64(4*ceil( maximum(abs.(z_p)) * maximum(abs.(times)) / (2*pi) ))
  if adaptK
    K = ceil(Int64, alpha*N + 2*m)
    # println("Calculated K for NFFT: ",K)
  end

  win, win_hat = NFFT.getWindow(:kaiser_bessel)

  return win, win_hat, K, N, t_scale, times, z_p
end

function get_A_coefficients_nfft(K::Int64
                                , times::Vector     # readout times
                                , z_p::Vector{Complex{T}}       # Correction map
                                , alpha::Float64    # oversampling factor(NFFT)
                                , m::Float64        # kernel size(NFFT)
                                ; adaptK::Bool=true
                                , centering::Bool=true) where T

  win, win_hat, K, N, t_scale, times, z_p = _precompute_coeffs_nfft(K, times, z_p, alpha, m; adaptK=adaptK,centering=centering)

  A = zeros(Complex{T}, length(times), K)

  for kappa=1:K
    for j=1:length(times)
        argument = times[j]/t_scale - (kappa - K/2)/(K - 2*m)
        A[j,kappa] = win(argument,N*alpha,m,alpha)
    end
  end
  return A
end

function get_C_coefficients_nfft(K::Int64
                                , times::Vector     # readout times
                                , z_p::Vector{Complex{T}}       # Correction map
                                , alpha::Float64    # oversampling factor(NFFT)
                                , m::Float64       # kernel size(NFFT)
                                ; adaptK::Bool=true
                                , centering::Bool=true) where T

  win, win_hat, K, N, t_scale, times, z_p = _precompute_coeffs_nfft(K, times, z_p, alpha, m; adaptK=adaptK,centering=centering)

  C = zeros(Complex{T}, K, length(z_p))

  for p=1:length(z_p)
    for kappa=1:K
        fourierarg = 1im * z_p[p]*t_scale/(2pi)
        exparg = -z_p[p]*t_scale * ( (kappa - K/2)/(K-2*m) )
        C[kappa,p] = (1/ ( (K-2*m)* 1/(N*alpha)*win_hat(fourierarg,N*alpha,m,alpha))) * exp(exparg)
    end
  end
  return C
end


"""
  Returns the a_j,k and c_k,p Coefficients considering the approximation of the
  exponential term exp(-t_j * z_p) = sum( k : a_j,k * c_k,p )

  Using the nfft method
"""
function get_AC_coefficients_nfft(K::Int64
                                , times::Vector     # readout times
                                , z_p::Vector{Complex{T}}       # Correction map
                                , alpha::Float64   # oversampling factor(NFFT)
                                , m::Float64        # kernel size(NFFT)
                                ; adaptK::Bool=true
                                , centering::Bool=true) where T

  A = get_A_coefficients_nfft(K, times, z_p, alpha, m; adaptK=adaptK, centering=centering)
  C = get_C_coefficients_nfft(K, times, z_p, alpha, m; adaptK=adaptK, centering=centering)
  return A, C
end


"""
  Returns the a_j,k and c_k,p Coefficients considering the approximation of the
  exponential term exp(-t_j * z_p) = sum( k : a_j,k * c_k,p )

  Using the least squares method and time-segmentation
"""
function get_AC_coefficients_lsqr(K::Int64, times::Vector{T}, z_p::Vector{Complex{T}}) where T
  A = zeros(Complex{T},length(times),K)
  C = zeros(Complex{T},K,length(z_p))
  G = zeros(Complex{T},length(z_p),K)

  t_hat = zeros(K)
  t_hat = collect(range(times[1], stop=times[length(times)], length=K))

  z_p = vec(z_p)
  for kappa=1:K
    for p=1:length(z_p)
      G[p,kappa] = exp(-t_hat[kappa]*z_p[p])
      C[kappa,p] = exp(-t_hat[kappa]*z_p[p])
    end
  end

  systemMat = G'*G
  b = [exp(-times[j]*z_p[p]) for p=1:length(z_p), j=1:length(times)]
  At = systemMat \ (G'*b)
  A .= transpose(At)
  
  return A,C
end

"""
  Returns the a_j,k Coefficients considering the approximation of the
  exponential term exp(-t_j * z_p) = sum( k : a_j,k * c_k,p )

  Using the least squares method with histogramms and time-segmentation
"""
function get_A_coefficients_hist_lsqr(K::Int64, times::Vector, z_p::Vector{Complex{T}}) where T

  numBins = 10*K
  if real.(z_p)== zeros(length(z_p))
    fmax_im = maximum(imag.(z_p))
    fmin_im = minimum(imag.(z_p))
    δf = (fmax_im-fmin_im)/(numBins-1)
    edges_im = [fmin_im-δf/2 + i*δf for i = 0:numBins]
    h = fit(Histogram, imag.(z_p),edges_im,closed=:left)
    z_center = 0.5*1im*vec( [edges_im[i]+edges_im[i+1] for i=1:numBins] )
    z_weights = h.weights
  elseif imag.(z_p)== zeros(length(z_p))
    fmax_re = maximum(real.(z_p))
    fmin_re = minimum(real.(z_p))
    δf = (fmax_re-fmin_re)/(numBins-1)
    edges_re = collect(fmin_re-δf/2.:δf:fmax_re+δf/2)
    h = fit(Histogram, real.(z_p),edges_re,closed=:left)
    z_center = 0.5*vec( [edges_re[i]+edges_re[i+1] for i=1:numBins] )
    z_weights = h.weights
  else
    nbins = ceil(Int, sqrt(numBins))
    fmax_re = maximum(real.(z_p))
    fmin_re = minimum(real.(z_p))
    δf_re = (fmax_re-fmin_re)/(nbins-1)
    edges_re = collect(fmin_re-δf/2.:δf:fmax_re+δf/2)
    fmax_im = maximum(imag.(z_p))
    fmin_im = minimum(imag.(z_p))
    δf = (fmax_im-fmin_im)/(nbins-1)
    edges_im = collect(fmin_im-δf/2.:δf:fmax_im+δf/2)
    h = fit(Histogram, (real.(z_p),imag.(z_p)),(edges_re,edges_im),closed=:left)
    z_center = 0.5*vec( [edges_re[i]+edges_re[i+1]+1im*(edges_im[j]+edges_im[j+1]) for i=1:numBins[1], j=1:numBins[2] ] )
    z_weights = vec(h.weights)
    numBins = nbins^2
  end

  A = zeros(Complex{T},length(times),K)
  G = zeros(Complex{T},numBins,K)

  t_hat = zeros(K)
  t_hat = collect(range(times[1], stop=times[end], length=K))

  z_p = vec(z_p)
  for p=1:numBins
    for kappa=1:K
      G[p,kappa] = exp(-t_hat[kappa]*z_center[p])
    end
  end

  systemMat = G'*Diagonal(z_weights)*G
  b = [z_weights[p]*exp(-times[j]*z_center[p]) for p=1:numBins, j=1:length(times)]
  At = systemMat \ (G'*b)
  A .= transpose(At)
    
  return A
end


"""
  Returns the c_k,p Coefficients considering the approximation of the
  exponential term exp(-t_j * z_p) = sum( k : a_j,k * c_k,p )

  Using the least squares method with histogramms
"""
function get_C_coefficients_hist_lsqr(K::Int64, times::Vector, z_p::Vector{Complex{T}}) where T

  C = zeros(Complex{T},K,length(z_p))
  t_hat = collect(range(times[1],stop=times[end],length=K))
  # z_p = vec(z_p)
  for p=1:length(z_p)
    for kappa=1:K
      C[kappa,p] = exp(-t_hat[kappa]*z_p[p])
    end
  end

  return C
end


"""
  Returns the a_j,k and c_k,p Coefficients considering the approximation of the
  exponential term exp(-t_j * z_p) = sum( k : a_j,k * c_k,p )

  Using the least squares method and frequency-segmentation
"""
function get_AC_coefficients_lsqr_fs(K::Int64, times::Vector, z_p::Vector{Complex{T}}) where T

  ω = imag.(z_p) # field map
  # ω_hat = collect(range(minimum(ω), stop=maximum(ω),length=K)) # translate frequencies
  ω_hat = loydMax(vec(ω),K, iterations=1000)

  A = zeros(Complex{T},length(times),K)
  G = zeros(Complex{T},length(times),K)
  
  for kappa=1:K
    for j=1:length(times)
      G[j,kappa] = exp(-1im*times[j]*ω_hat[kappa])
      A[j,kappa] = exp(-1im*times[j]*ω_hat[kappa])
    end
  end

  systemMat = G'*G
  b = [exp(-times[j]*z_p[p]) for j=1:length(times), p=1:length(z_p)]
  C = systemMat \ (G'*b)
    
  return A,C
end

"""
  Returns the a_j,k Coefficients considering the approximation of the
  exponential term exp(-t_j * z_p) = sum( k : a_j,k * c_k,p )

  Using the least squares method with histogramms and frequency-segmentation
"""
function get_A_coefficients_hist_lsqr_fs(K::Int64, times::Vector, z_p::Vector{Complex{T}}) where T
  ω = imag.(z_p) # field map
  # ω_hat = collect(range(minimum(ω), stop=maximum(ω),length=K)) # translate frequencies
  ω_hat = loydMax(vec(ω),K, iterations=1000)

  A = zeros(Complex{T},length(times),K)
  for kappa=1:K
    for j=1:length(times)
      A[j,kappa] = exp(-1im*times[j]*ω_hat[kappa])
    end
  end

  return A
end

function get_A_coefficients_svd(K::Int64, times::Vector, z_p::Vector{Complex{T}}; numBins=10*K, K_tol::T=1.e-4) where T
  ω = imag.(z_p)

  # create histrogram of field map
  fmax = maximum(ω)
  fmin = minimum(ω)
  δf = (fmax-fmin)/(numBins-1)
  edges = [fmin-δf/2 + i*δf for i = 0:numBins]
  h = fit( Histogram, ω, edges, closed=:left )
  z_center = 0.5*1im * ( edges[1:numBins].+edges[2:numBins+1] )
  z_weights = h.weights

  # from reduced exponential matrix
  A = zeros(Complex{T}, length(times), numBins)
  for j=1:numBins, i=1:length(times)
    A[i,j] = sqrt(z_weights[j])*exp(-z_center[j]*times[i])
  end

  # compute SVD
  U,S,V = psvd(A, rank=2*K)
  # estimate accuracy of the approximation and determine rank threshold
  acc = sqrt.(cumsum(S.^2)./sum(S.^2))
  K0 = findfirst(x->x>=(1-K_tol-1e-14), acc)
  @debug "accuracy: $(acc)"
  @info " K_tol=$(K_tol) => chose K= $(K0)"
  
  return U[:,1:K0]
end

"""
  Returns the c_k,p Coefficients considering the approximation of the
  exponential term exp(-t_j * z_p) = sum( k : a_j,k * c_k,p )

  Uses the LSQR-method for a given set of a_j,k-coefficients
"""
function get_C_coefficients_lsqr(times::Vector, z_p::Vector{Complex{T}}, A::Matrix{Complex{T}}, numSamp::Int64, step=10) where T

  systemMat = adjoint(A[1:step:numSamp,:])* A[1:step:numSamp,:]
  b = [exp(-times[j]*z_p[i]) for j=1:step:numSamp, i=1:length(z_p)]
  C = systemMat \ ( adjoint(A[1:step:numSamp,:])*b )
  
  return C
end
