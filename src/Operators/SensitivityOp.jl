function prod_smap!(y::Vector{ComplexF64}, smaps::Matrix{ComplexF64}, x::Vector{ComplexF64}, numVox::Int64, numChan::Int64, numContr::Int=1)
    x = reshape(x,numVox,numContr)
    y = reshape(y,numVox,numContr,numChan)
    @inbounds for j=1:numChan
        sj = (@view smaps[:,j])
        @inbounds for i=1:numContr
            y[:,i,j] .= (@view x[:,i]).* sj
        end
    end
end

function ctprod_smap!(y::Vector{ComplexF64}, smapsC::Matrix{ComplexF64}, x::Vector{ComplexF64}, numVox::Int64, numChan::Int64, numContr::Int=1)
    x = reshape(x,numVox,numContr,numChan)
    y = reshape(y,numVox,numContr)
    y .= 0
    @inbounds for j=1:numChan
        sj = @view smapsC[:,j]
        @inbounds for i=1:numContr
            y[:,i] += (@view x[:,i,j]).* sj
        end
    end
end

function prod_smap(smaps::Matrix{ComplexF64}, x::Vector{ComplexF64}, numVox::Int64, numChan::Int64, numContr::Int=1)
    y = zeros(ComplexF64,numVox*numContr*numChan)
    prod_smap!(y,smaps,x,numVox,numChan,numContr)
    return y
end

function ctprod_smap(smapsC::Matrix{ComplexF64}, x::Vector{ComplexF64}, numVox::Int64, numChan::Int64, numContr::Int=1)
    y = zeros(ComplexF64,numVox*numContr)
    ctprod_smap!(y,smapsC,x,numVox,numChan,numContr)
    return y
end

"""
    SensitivityOp(sensMaps::Matrix{ComplexF64}, numEchoes::Int=1)

builds a `LinearOperator` which performs multiplication of a given image with
the coil sensitivities specified in `sensMaps`

# Arguments
* `sensMaps::Matrix{ComplexF64}`  - sensitivity maps ( 1. dim -> voxels, 2. dim-> coils)
* `numEchoes`                     - number of contrasts to which the opetaor will be applied
"""
function SensitivityOp(sensMaps::Matrix{ComplexF64}, numContr::Int=1)
    numVox, numChan = size(sensMaps)
    sensMapsC = conj.(sensMaps)
    return LinearOperator(numVox*numContr*numChan,numVox*numContr,false,false,
                         x->prod_smap(sensMaps,x,numVox,numChan,numContr),
                         nothing,
                         x->ctprod_smap(sensMapsC,x,numVox,numChan,numContr))
end

"""
    SensitivityOp(sensMaps::Array{T,4}, numContr::Int=1) where T

builds a `LinearOperator` which performs multiplication of a given image with
the coil sensitivities specified in `sensMaps`

# Arguments
* `sensMaps::Array{T,4}`  - sensitivity maps ( 1.-3. dim -> voxels, 4. dim-> coils)
* `numContr`             - number of contrasts to which the opetaor will be applied
"""
function SensitivityOp(sensMaps::Array{T,4}, numContr::Int=1) where T #{T,D}
  sensMaps_mat = reshape(sensMaps, div(length(sensMaps),size(sensMaps,4)),size(sensMaps,4))
  return SensitivityOp(sensMaps_mat, numContr)
end
