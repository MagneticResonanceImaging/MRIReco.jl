export regrid

"""
    `regrid(acqData::AcquisitionData{T}, kspaceSize::NTuple{2,Int64}
              ; cgnr_iter::Int64=3, correctionMap::Array{Complex{T}}=Complex{T}[]) where (D,T)`

Regrid non-cartesian k-space data in `acqData` to a cartesian grid of size `kspaceSize`.
Uses the CGNR method to invert the non-cartesian Fourier encoding.

# Arguments:
* `acqData::AcquisitionData{T}`        - AcquisitionData
* `ksize::NTuple{2,Int64}`          - size of the cartesian k-space grid
* `cgnr_iter::Int64=3`              - number of CGNR iterations
* `correctionMao::Array{Complex{T}}`- relaxation/b0 map
"""
function regrid(acqData::AcquisitionData{T,2}, kspaceSize::NTuple{2,Int64}; 
                cgnr_iter::Int64=3, correctionMap::Array{Complex{T}}=Complex{T}[]) where {T}
  dcf = samplingDensity(acqData, kspaceSize) 
  numContr, numChan, numSl = numContrasts(acqData), numChannels(acqData), numSlices(acqData)

  kdata_cart = [zeros(Complex{T}, prod(kspaceSize), numChan) for j=1:numContr, k=1:numSl, rep=1:1]
  F = T(sqrt(prod(kspaceSize)))*FFTOp(Complex{T}, shape=kspaceSize)
  img = zeros(Complex{T}, prod(kspaceSize))
  for k = 1:numSl
    E = encodingOps_simple(acqData, kspaceSize, slice=k, correctionMap=correctionMap)
    for j = 1:numContr
      W = WeightingOp(Complex{T}; weights=dcf[j])
      solver = RegularizedLeastSquares.CGNR(W*E[j], iterations=cgnr_iter)
      for i = 1:numChan
        kdata = kData(acqData,j,i,k).* dcf[j]
        img .= solve!(solver, kdata)
        kdata_cart[j,k,1][:,i] .= F*img
      end
    end
  end

  return AcquisitionData(CartesianTrajectory(T,kspaceSize[2],kspaceSize[1]), kdata_cart, encodingSize=kspaceSize, fov=acqData.fov)
end
