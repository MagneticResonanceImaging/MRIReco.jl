
"""
  With multiple Images as input, determines Amplitude and Phase of magnetisation
  , fieldmap and relaxationMap.

"""
function reconstructionWithEchoImages(echoImages
                                      , echoTimes
                                      ; ampFunc="minsquare"
                                      , phaseFunc="nonlinear"
                                      , verbose::Bool=true)
  # Get height width and number of echoimages
  height, width, numOfEchos = size(echoImages)
  if length(echoTimes) != numOfEchos
    error("Number of echo Images and length of echo times are not consistent!")
  end

  # Preallocate output
  amplResult = zeros(height,width)
  relaxResult = zeros(height,width)
  omegaResult = zeros(height,width)
  phiResult = zeros(height,width)

  phi = zeros(numOfEchos)
  k = zeros(numOfEchos)
  d = zeros(numOfEchos)

  if verbose==true
    progr = Progress(size(echoImages,1)*size(echoImages,2), 1, "Solving Nonlinear equations...")
  end

  # Solver nonlinear equations independently for each pixel
  for i=1:size(echoImages,1)
    for j=1:size(echoImages,2)
      # ----- Solving a nonlinear function(unconstrained) to determine amplitude of magentisation and relaxationmap
      y_amp = [abs(echoImages[i,j,l]) for l=1:numOfEchos]
      # Determine phase detection of echoimages
      for l=1:numOfEchos
        phi[l] = atan2( imag(echoImages[i,j,l]), real(echoImages[i,j,l]) )
      end
      # Amplitude of magentisation and relaxationrate
      amplResult[i,j], relaxResult[i,j] = determineAmplitudeAndRelaxation_optim(y_amp
                                                                    , echoTimes
                                                                    ; ampFunc=ampFunc)
      # Phase of magnetisation and offresonance
      omegaResult[i,j], phiResult[i,j] = determinePhaseAndOffset(phi
                                                                , echoTimes
                                                                , numOfEchos
                                                                ; method=phaseFunc)
      if verbose==true
          next!(progr)
      end
    end
  end

  # Returns |y|*exp(i*phi) , z
  return amplResult .* exp(1im*phiResult) , relaxResult + 1im*omegaResult

end


#
# Solves leastsquare problem of amplitude and Phase
#
function determineAmplitudeAndRelaxation_optim(y_amp::Vector
                                          , tau::Vector
                                          ; ampFunc="minsquare")

    potStartvectors = getStartPositions(y_amp, tau)
    # Guessing three possible start vectors for solving the nonlinear
    # equation
    if ampFunc == "minsquare"
        y = [minSquareAmplitude_optim(potStartvectors[1,:], y_amp, tau)
            ,minSquareAmplitude_optim(potStartvectors[2,:], y_amp, tau)
            ,minSquareAmplitude_optim(potStartvectors[3,:], y_amp, tau)]
    elseif ampFunc == "chisquare"
        sigma = 0.1 # choosing sigma at random...
        y = [chiSquareAmplitude_optim(potStartvectors[1,:], y_amp, tau, sigma)
            ,chiSquareAmplitude_optim(potStartvectors[2,:], y_amp, tau, sigma)
            ,chiSquareAmplitude_optim(potStartvectors[3,:], y_amp, tau, sigma)]
    else
        error("AmpFunc not found: ", ampFunc)
    end
    # Choose the startvector with the minimal cost
    minIndx = indmin(y)
    # Optimize function
    result = Optim.optimize(x->minSquareAmplitude_optim(x,y_amp,tau)
                          , potStartvectors[minIndx,:]
                          , LBFGS()
                          );
    return Optim.minimizer(result)
end


"""
  Determines the phase of the magnetisation and imaginary part (offResonance)
  of the fieldmap given multiple echo-images
"""
function determinePhaseAndOffset(phi::Vector
                                , tau::Vector
                                , numOfEchos
                                ; method="nonlinear")

  # Using start vector with the help of the first two given echo Images
  w0 = atan2( imag(exp(1im *(phi[1]-phi[2]))) ,real(exp(1im *(phi[1]-phi[2]))) ) / (tau[2]-tau[1])
  phip0 = atan2( imag( exp(1im*(phi[1]+tau[1]*w0)) ) , real( exp(1im*(phi[1]+tau[1]*w0)) ))

  k = zeros(numOfEchos)
  d = zeros(numOfEchos)
  for l=1:numOfEchos
    k[l] = (-phi[l]-tau[l]*w0+phip0)/(2*pi)
    d[l] =  1.# y_amp[l]^2 <- not good weight if image contains zeros, may lead to NaNs
  end

  if method == "linear"
    # Precalculation for solution
    nom = sum(d.*tau .* (phi + 2*pi*k)) * sum(d) - sum(d.* (phi + 2*pi*k)) * sum(d.*tau)
    den = sum(d.*tau)^2 - sum(d .* tau .* tau) * sum(d)
    # Direct linear solution
    omegaResult = nom / den
    phiResult = (sum(d .* (phi + 2*pi*k)) + omegaResult * sum(d.*tau)) / sum(d)
  elseif method == "nonlinear"
      result  = Optim.optimize(x->minSquarePhase_optim(x,d,tau,phi)
                            , [w0,phip0]
                            , LBFGS()
                            );
      omegaResult,phiResult = Optim.minimizer(result)
  end
  return omegaResult,phiResult
end

"""
  Joint estimation of magnetisation and fieldmap given multiple echo kdata
"""
function reconstructionWMultiEchoData(kdata
                                      , echoTimes::Vector
                                      , shape::Tuple
                                      ; aqTime::Float64=10e-3
                                      , alpha::Float64=1.75
                                      , m::Float64=4.
                                      , K::Int64=20
                                      , traj="cartesian"
                                      , coeffMethod="nfft")
  numOfPixel, numOfNodes = prod(shape), prod(shape)
  L = length(echoTimes)
  p = zeros(ComplexF64,shape)

  # Preallocating output
  z, y = zeros(ComplexF64,shape), zeros(ComplexF64,shape)
  z_old, y_old = zeros(ComplexF64,shape), zeros(ComplexF64,shape)
  echoImg = zeros(ComplexF64,shape[1],shape[2],L)

  loopCounter= 1
  r = Inf
  while true
    if loopCounter > 10
        println("Maximum of iterations reached, exiting calculation...")
        break
    end
    # Solving echo magentisation using CGNR and applying known correctionterm
    for i=1:L
      if loopCounter == 1
        # If very first iteration use NFFT Operator without any correctionterm
        # since it is unknown
        tr = getTrajectory(traj, shape, 0.0, 0.0) # no times required here since at z=0 time does not matter
        nftOp = NFFTOperator((shape[1],shape[2]),tr )
        f0 = Ac_mul_W_mul_b(nftOp,sdc(nftOp.plan), kdata[i,:])
        fL = cgnr(nftOp, reshape(f0,numOfPixel), kdata[i,:], sdc(nftOp.plan), 30,verboseRes=false)
        echoImg[:,:,i] = reshape(fL,shape)
      else
        # Important: no echoTimes given here since we want to reconstruct echoImages
        # and not the actual image
        tr = getTrajectory(traj,shape,0.0, aqTime)
        nftWCOp = NFFTWCorrectionOperator(ComplexF64,(shape[1],shape[2]),tr,z,method=coeffMethod)
        f0 = Ac_mul_W_mul_b(nftWCOp,sdc(nftWCOp.plan),kdata[i,:])
        fl = cgnr(nftWCOp, reshape(f0,numOfPixel), kdata[i,:], sdc(nftWCOp.plan), 30,verboseRes=false)
        echoImg[:,:,i] = reshape(fl,shape)
      end
    end
    ############################################################################
    println("Getting echoImages is finished ")
    # Second block of joint algorithm consists of (re-)estimating correctionterm
    # apmlutide and phase of magnetisation
    y_old = y
    z_old = z
    # Calculating Magentisation and Correctionterm with non linear equations
    if L == 2
        y,z = determineMagnetisaitonAndCorrection(shape,echoTimes,echoImg)
    else
        y,z = reconstructionWithEchoImages(echoImg,echoTimes;ampFunc="minsquare",phaseFunc="linear")
    end
    println("Non linear problems solved")
    # Filtering the Relaxationmap and unwrapping the offResonancemap
    relax_new  = median_filter(real(z))
    offres_new = unwrap_phase(imag(z))
    z = relax_new + 1im*offres_new

    # Calculating the residuum
    newRes = 0.0
    for i=1:L
        # P_z is an addtional vector to "echo" the calculated image
        P_z = exp(-echoTimes[i] .* z)   #'echofying' the original img
        echimg = y .* P_z
        tr = getTrajectory(traj,shape,0.0, aqTime)
        nftWCOp = NFFTWCorrectionOperator(ComplexF64,(shape[1],shape[2]),tr,z,method=coeffMethod)
        # Calculating kdata of calculated correctionmap and magnetisation
        calc_kdata = nftWCOp * vec(echimg)
        newRes +=  norm( (kdata[i,:] - calc_kdata) .* sqrt(nftWCOp.density) )^2
    end
    println("New Residuum: ", newRes)
    loopCounter +=1
    # Stopping criteria
    if newRes > r || abs(newRes-r) / (2*abs(newRes+r)) < 1e-4
      y = y_old
      z = z_old
      println("New calculated residuum is higher than the old one, or fixpoint found, exiting calculation.")
      break
    else
      r = newRes
    end
  end
  return y,z
end

"""
    Determine magentisaiton and correction with 2 echotimes directly

"""
function determineMagnetisaitonAndCorrection(imgSize::Tuple
                                            , echoTimes::Vector
                                            , echoImg)
    if length(echoTimes) != 2
        error("This function can only solve magentisaiton and correction if 2 echoImages are provided")
    end
    correction = zeros(ComplexF64,imgSize)
    magnetisation = zeros(ComplexF64,imgSize)
    for i=1:imgSize[1]
        for j=1:imgSize[2]
            correction[i,j] = log(echoImg[i,j,1]/echoImg[i,j,2])  / (echoTimes[2]-echoTimes[1])
            magnetisation[i,j] = echoImg[i,j,1] / exp(-echoTimes[1] * correction[i,j] )
        end
    end
    return magnetisation, correction
end


########################## Helper functions ####################################

# Median filter to eliminate ausreisser
function median_filter(im::Vector, filterSize=3)
 N = size(im)
 out = similar(im)
 K = floor(Int,filterSize/2)
 for x=1:N[1]
   x_min = max(1, x-K)
   x_max = min(N[1], x+K)

   s = im[x_min:x_max]
   out[x] = median(s[:])
 end
 out
end

function median_filter(im::Matrix, filterSize=3)
 N = size(im)
 out = similar(im)
 K = floor(Int,filterSize/2)
 for x=1:N[1], y=1:N[2]
   x_min = max(1, x-K)
       x_max = min(N[1], x+K)
   y_min = max(1, y-K)
       y_max = min(N[2], y+K)

   s = im[x_min:x_max, y_min:y_max]
   out[x,y] = median(s[:])
 end
 out
end

# Start positions for Amplutide and Phase
function getStartPositions(y_amp::Vector, tau::Vector)
    # Guess three potential startvectors
    relaxationGuess = [ log(y_amp[end]/ y_amp[1]) / (tau[1]-tau[end])
                      , log(y_amp[2] / y_amp[1]) / (tau[1]-tau[2])
                      , log(y_amp[end-2]+y_amp[end-1]+y_amp[end] / (y_amp[1]+y_amp[2]+y_amp[3])  ) / (tau[1] - tau[end-1])
                      ]

    potStartvectors = [ [y_amp[1]*exp(tau[1]*relaxationGuess[1]) relaxationGuess[1]]
                      ; [y_amp[1]*exp(tau[1]*relaxationGuess[2]) relaxationGuess[2]]
                      ; [1/3 * (y_amp[1]+y_amp[2]+y_amp[3]) * exp(tau[2] * relaxationGuess[3]) relaxationGuess[3]]
                      ]
    return potStartvectors
end


function getTrajectory(traj,shape::Tuple,echoTime::T,aqTime::T) where T<:Real
    if traj=="cartesian"
      tr = SimpleCartesianTrajectory(shape[1],shape[2], echoTime, aqTime)
    elseif traj=="radial"
      tr = RadialTrajectory(shape[1] * 2, shape[2], echoTime, aqTime)
    elseif traj=="spiral"
      tr = SpiralTrajectory(32, (shape[1]*shape[2]) >> 4, echoTime, aqTime)
    else
      error("Unknown traj")
    end
    return tr
end
