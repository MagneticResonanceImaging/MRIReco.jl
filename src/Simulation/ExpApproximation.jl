
"""
  Returns the a_j,k and c_k,p Coefficients considering the approximation of the
  exponential term exp(-t_j * z_p) = sum( k : a_j,k * c_k,p )
"""
function getA_Ccoefficients(K::Int64
                            , numOfNodes::Int64
                            , numOfPixel::Int64
                            , times::Vector,z_p
                            , alpha::Float64
                            , m::Float64
                            ; method = "nfft"
                            , double = false)
  A = zeros(ComplexF64,numOfNodes,K)
  C = zeros(ComplexF64,K,numOfPixel)
  if method == "leastsquare"
    A,C = getA_Coefficients_least_Squares(K,numOfNodes,numOfPixel,times,z_p)
  elseif method == "nfft"
    A,C = getA_Ccoefficients_nfft(K,numOfNodes,numOfPixel,times,z_p,alpha,m)
  elseif method == "hist"
    A,C = getA_Coefficients_hist_lsqr(K,numOfNodes,numOfPixel,times,z_p)
  else
    error("Not implemented yet")
  end
  return A,C
end

"""
  Returns the a_j,k and c_k,p Coefficients considering the approximation of the
  exponential term exp(-t_j * z_p) = sum( k : a_j,k * c_k,p )

  Using the nfft method
"""
function getA_Ccoefficients_nfft(K::Int64
                                , numOfNodes::Int64
                                , numOfPixel::Int64
                                , times::Vector     # readout times
                                , z_p::Vector       # Correction map
                                , alpha::Float64    # oversampling factor(NFFT)
                                , m::Float64
                                ; double::Bool=false)       # kernel size(NFFT)

  # Centering time and correctionterm to lower computational effort
  t_hat = (times[1] + times[end])/2
  z_hat = minimum(real(z_p)) + maximum(real(z_p)) + 1im*(minimum(imag(z_p)) + maximum(imag(z_p)))
  z_hat *= 0.5
  # Centering readouttimes and correctionmap
  times = times .- t_hat
  z_p = z_p .- z_hat

  # Calculating timescaling factor so it satisfies: t_j / T in [0.5,0.5)
  T =  maximum(abs.(times)) / (0.5-1e-7)

  # determining N and K,
  # K can be choosen freely with tradeoff comp. complexity vs precision
  N = Int64(4*ceil( maximum(abs.(z_p)) * maximum(abs.(times)) / (2*pi) ))
  K = ceil(Int64, alpha*N + 2*m)
  # println("Calculated K for NFFT: ",K)

  # Preallocating output
  A = zeros(ComplexF64,numOfNodes, K)
  C = zeros(ComplexF64,K, numOfPixel)

  win, win_hat = NFFT.getWindow(:kaiser_bessel)

  for kappa=1:K
    for j=1:numOfNodes
        argument = times[j]/T - (kappa - K/2)/(K -2*m)
        A[j,kappa] = win(argument,N*alpha,m,alpha)
    end
  end

  for p=1:numOfPixel
    for kappa=1:K
        fourierarg = 1im * z_p[p]*T/(2pi)
        exparg = -z_p[p]*T * ( (kappa - K/2)/(K-2*m) )
        C[kappa,p] = (1/ ( (K-2*m)* 1/(N*alpha)*win_hat(fourierarg,N*alpha,m,alpha))) * exp(exparg)
    end
  end
  return A,C
end

function getA_coefficients_nfft(K::Int64
                                , numOfNodes::Int64
                                , numOfPixel::Int64
                                , times::Vector     # readout times
                                , z_p::Vector       # Correction map
                                , alpha::Float64    # oversampling factor(NFFT)
                                , m::Float64
                                ; double::Bool=false
                                , adaptK::Bool=false)       # kernel size(NFFT)

  # Centering time and correctionterm to lower computational effort
  t_hat = (times[1] + times[end])/2
  z_hat = minimum(real(z_p)) + maximum(real(z_p)) + 1im*(minimum(imag(z_p)) + maximum(imag(z_p)))
  z_hat *= 0.5
  # Centering readouttimes and correctionmap
  times = times .- t_hat
  z_p = z_p .- z_hat

  # Calculating timescaling factor so it satisfies: t_j / T in [0.5,0.5)
  T =  maximum(abs.(times)) / (0.5-1e-7)

  # determining N and K,
  # K can be choosen freely with tradeoff comp. complexity vs precision
  N = Int64(4*ceil( maximum(abs.(z_p)) * maximum(abs.(times)) / (2*pi) ))
  if adaptK
    K = ceil(Int64, alpha*N + 2*m)
  end
  # println("Calculated K for NFFT: ",K)

  # Preallocating output
  A = zeros(ComplexF64,numOfNodes, K)

  win, win_hat = NFFT.getWindow(:kaiser_bessel)

  for kappa=1:K
    for j=1:numOfNodes
        argument = times[j]/T - (kappa - K/2)/(K -2*m)
        A[j,kappa] = win(argument,N*alpha,m,alpha)
    end
  end
  return A
end

function getC_coefficients_nfft(K::Int64
                                , numOfNodes::Int64
                                , numOfPixel::Int64
                                , times::Vector     # readout times
                                , z_p::Vector       # Correction map
                                , alpha::Float64    # oversampling factor(NFFT)
                                , m::Float64
                                ; double::Bool=false
                                , adaptK::Bool=false)       # kernel size(NFFT)

  # Centering time and correctionterm to lower computational effort
  t_hat = (times[1] + times[end])/2
  z_hat = minimum(real(z_p)) + maximum(real(z_p)) + 1im*(minimum(imag(z_p)) + maximum(imag(z_p)))
  z_hat *= 0.5
  # Centering readouttimes and correctionmap
  times = times .- t_hat
  z_p = z_p .- z_hat

  # Calculating timescaling factor so it satisfies: t_j / T in [0.5,0.5)
  T =  maximum(abs.(times)) / (0.5-1e-7)

  # determining N and K,
  # K can be choosen freely with tradeoff comp. complexity vs precision
  N = Int64(4*ceil( maximum(abs.(z_p)) * maximum(abs.(times)) / (2*pi) ))
  if adaptK
    K = ceil(Int64, alpha*N + 2*m)
  end
  # println("Calculated K for NFFT: ",K)

  # Preallocating output
  C = zeros(ComplexF64,K, numOfPixel)

  win, win_hat = NFFT.getWindow(:kaiser_bessel)

  for p=1:numOfPixel
    for kappa=1:K
        fourierarg = 1im * z_p[p]*T/(2pi)
        exparg = -z_p[p]*T * ( (kappa - K/2)/(K-2*m) )
        C[kappa,p] = (1/ ( (K-2*m)* 1/(N*alpha)*win_hat(fourierarg,N*alpha,m,alpha))) * exp(exparg)
    end
  end
  return C
end

"""
  Returns the a_j,k and c_k,p Coefficients considering the approximation of the
  exponential term exp(-t_j * z_p) = sum( k : a_j,k * c_k,p )

  Using the least squares method
"""
function getA_Coefficients_least_Squares(K::Int64
                                        , numOfNodes::Int64
                                        , numOfPixel::Int64
                                        , times::Vector
                                        , z_p::Vector)
  A = zeros(ComplexF64,numOfNodes,K)
  C = zeros(ComplexF64,K,numOfPixel)
  G = zeros(ComplexF64,numOfPixel,K)

  progr = Progress(numOfNodes, 1, "Using leastsquares to get Coefficients...")

  t_hat = zeros(K)
  t_hat = collect(range(times[1], stop=times[numOfNodes], length=K))

  z_p = vec(z_p)
  for kappa=1:K
    for p=1:numOfPixel
      G[p,kappa] = exp(-t_hat[kappa]*z_p[p])
      C[kappa,p] = exp(-t_hat[kappa]*z_p[p])
    end
  end

  systemMat = G'*G
  err = zeros(numOfNodes)
  for j=1:numOfNodes
    b = [exp(-times[j]*z_p[p]) for p=1:numOfPixel ]
    A[j,:] = systemMat \ (G'*b)
    next!(progr)
  end
  return A,C
end

"""
  Returns the a_j,k Coefficients considering the approximation of the
  exponential term exp(-t_j * z_p) = sum( k : a_j,k * c_k,p )

  Using the least squares method with histogramms
"""
function getA_Coefficients_hist_lsqr(K::Int64
                                        , numOfNodes::Int64
                                        , numOfPixel::Int64
                                        , times::Vector
                                        , z_p::Vector)

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

  A = zeros(ComplexF64,numOfNodes,K)
  G = zeros(ComplexF64,numBins,K)

  # progr = Progress(numOfNodes, 1, "Using leastsquares to get Hist-coefficients...")

  t_hat = zeros(K)
  # t_hat = collect(linspace(times[1],times[numOfNodes],K))
    t_hat = collect(range(times[1], stop=times[end], length=K))

  z_p = vec(z_p)
  for p=1:numBins
    for kappa=1:K
      G[p,kappa] = exp(-t_hat[kappa]*z_center[p])
    end
  end

  # systemMat = G'*diagm(z_weights)*G
  systemMat = G'*Diagonal(z_weights)*G
  err = zeros(numOfNodes)
  @showprogress 1 "Using leastsquares to get Hist-coefficients..." for j=1:length(times) #numOfNodes
    b = [z_weights[p]*exp(-times[j]*z_center[p]) for p=1:numBins ]
    A[j,:] = systemMat \ (G'*b)
    # next!(progr)
  end
  return A
end


"""
  Returns the c_k,p Coefficients considering the approximation of the
  exponential term exp(-t_j * z_p) = sum( k : a_j,k * c_k,p )

  Using the least squares method with histogramms
"""
function getC_Coefficients_hist_lsqr(K::Int64
                                        , numOfNodes::Int64
                                        , numOfPixel::Int64
                                        , times::Vector
                                        , z_p::Vector)

  C = zeros(ComplexF64,K,numOfPixel)
  t_hat = collect(range(times[1],stop=times[end],length=K))
  # z_p = vec(z_p)
  for p=1:numOfPixel
    for kappa=1:K
      C[kappa,p] = exp(-t_hat[kappa]*z_p[p])
    end
  end

  return C
end

function getA_Coefficients_one_term(K::Int64
                                        , numOfNodes::Int64
                                        , numOfPixel::Int64
                                        , times::Vector
                                        , z_p::Vector)
  A = zeros(ComplexF64,numOfNodes,K)
  t_hat = collect(range(times[1],stop=times[end],length=K))
  for p=1:numOfNodes
    diff, idx = findmin(abs.(t_hat .- times[p]))
    A[p,idx]=1.0
  end
  return A
end

"""
  Returns the a_j,k Coefficients considering the approximation of the
  exponential term exp(-t_j * z_p) = sum( k : a_j,k * c_k,p )

  Using the least squares method with histogramms
"""
function getA_Coefficients_two_terms(K::Int64
                                        , numOfNodes::Int64
                                        , numOfPixel::Int64
                                        , times::Vector
                                        , z_p::Vector)

  # build histogram
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

  A = zeros(ComplexF64,numOfNodes,K)
  G = zeros(ComplexF64,numBins,K)

  progr = Progress(numOfNodes, 1, "Using leastsquares to get Hist-coefficients...")

  # t_hat = zeros(K)
  t_hat = collect(range(times[1],stop=times[end],length=K))

  z_p = vec(z_p)
  for p=1:numBins
    for kappa=1:K
      G[p,kappa] = exp(-t_hat[kappa]*z_center[p])
    end
  end

  for j=1:length(times) #numOfNodes
    # find closest translate with smaller sampling time
    diff, idx = findmin(abs.(t_hat .- times[j]))
    if (t_hat[idx] >= times[j]) && idx > 1
      idx -= 1
    end

    systemMat = G[:,idx:idx+1]'*diagm(0 => z_weights)*G[:,idx:idx+1]
    b = [z_weights[p]*exp(-times[j]*z_center[p]) for p=1:numBins ]
    A[j,idx:idx+1] = systemMat \ (G[:,idx:idx+1]'*b)
    next!(progr)
  end
  return A
end
