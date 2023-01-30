export decorrelateSenseMaps, covariance

"""
        decorrelateSenseMaps(L_inv, senseMaps, numSl)

    multiplies the senseMaps by the noise uncorrelation matrix (L_inv).
"""

function decorrelateSenseMaps(L_inv::Union{LowerTriangular{Complex{T}, Matrix{Complex{T}}}, Nothing},
                    senseMaps::Array{Complex{T},4},
                    numChan::Int64) where {T}
    if isnothing(L_inv)
        senseMapsUnCorr = senseMaps
    else
        sizeMaps = size(senseMaps)
        senseMapsUnCorr = reshape(senseMaps, :, numChan) * L_inv'
        senseMapsUnCorr = reshape(senseMapsUnCorr, sizeMaps)
    end

    return senseMapsUnCorr
end


"""
        covariance(noiseData)

  computes the covariance of the noise acquisition.
"""

function covariance(noiseData::Array{Complex{T},2}) where {T}

    N = size(noiseData, 1)
    cov = (1/(N-1)) .* (noiseData' * noiseData)

    return cov
end