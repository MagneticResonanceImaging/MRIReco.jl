export decorrelateSenseMaps, covariance

"""
    decorrelateSenseMaps(L_inv, senseMaps, numChan)

multiplies the coil sensitivities `senseMaps`, an array of size (x, y, z, coils) with `numChan` coils, by `L_inv'`
from the right along the coil dimension, `reshape(senseMaps, :, numChan) * L_inv'`, which decorrelates them in the
same way as the data. `L_inv` is the inverse of the lower Cholesky factor of the matrix that `covariance` returns.
If `L_inv` is `nothing`, `senseMaps` is returned unchanged.
"""
function decorrelateSenseMaps(L_inv::Union{LowerTriangular{Complex{T}, <:AbstractMatrix{Complex{T}}}, Nothing},
                    senseMaps::AbstractArray{Complex{T},4},
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

computes the covariance matrix of the coils from the noise samples `noiseData`, an array of size (samples, coils),
as `noiseData' * noiseData / (N - 1)` for `N` samples.
"""
function covariance(noiseData::AbstractArray{Complex{T},2}) where {T}

    N = size(noiseData, 1)
    cov = (1/(N-1)) .* (noiseData' * noiseData)

    return cov
end