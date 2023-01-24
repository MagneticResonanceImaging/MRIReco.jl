"""
        noiseDeCorr!(acqData, noiseData, senseMaps, numContr, numSl, numRep)

    uncorrelates noise in the sense maps and acquisition data (noise pre-whitening)
"""

function noiseDeCorr!(acqData::AcquisitionData{T,D},
                    noiseData::Array{Complex{T},2},
                    senseMaps::Array{Complex{T},4},
                    numCont::Int64,
                    numSl::Int64,
                    numRep::Int64)

    L_inv = noiseDeCorrMat(noiseData)
    ApplytoMaps!(L_inv, senseMaps, numSl)
    #ApplytoData!(L_inv, acqData, numCont, numSl, numRep)

end


"""
    UnCorrSenseMaps(L_inv, senseMaps, numSl)

    multiplies the senseMaps by the noise uncorrelation matrix (L_inv).
"""

function UnCorrSenseMaps(L_inv::LowerTriangular{ComplexF64, Matrix{ComplexF64}},
                    senseMaps::Array{Complex{T},4},
                    numChan::Int64)
    
    sizeMaps = size(senseMaps)
    senseMapsUnCorr = reshape(senseMaps, :, numChan) * L_inv'
    senseMapsUnCorr = reshape(senseMapsUnCorr, sizeMaps)

    return senseMapsUnCorr
end


"""
    ApplytoData!(L_inv, acqData)

    multiplies the acqData by the noise uncorrelation matrix (L_inv).
"""

function ApplytoData!(L_inv::LowerTriangular{ComplexF64, Matrix{ComplexF64}},
                    acqData::AcquisitionData{T,D},
                    numCont::Int64,
                    numSl::Int64,
                    numRep::Int64)
    
    for l = 1:numRep, k = 1:numSl, j = 1:numCont
        acqdData.kdata[j, k, l] = acqdData.kdata[j, k, l] * L_inv'
    end

end


"""
        noiseDeCorrMat(noiseData)

    computes the noise uncorrelation matrix (L_inv) from the noise acquisition.
"""

function noiseDeCorrMat(noiseData::Array{Complex{T},2})

    psi = cov(noiseData)
    L = cholesky(psi, check = true)
    L_inv = inv(L.L)

    return L_inv
end


"""
        covariance(noiseData)

  computes the covariance of the noise acquisition.
"""

function covariance(noiseData::Array{Complex{T},2})

    N = size(noiseData, 1)
    cov = (1/(N-1)) .* (noiseData' * noiseData)

    return cov
end