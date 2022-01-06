export regrid2d

"""
    `regrid2d(acqData::AcquisitionData, kspaceSize::NTuple{2,Int64}
              ; cgnr_iter::Int64=3, correctionMap::Array{ComplexF64}=ComplexF64[]) where D`

Regrid non-cartesian k-space data in `acqData` to a cartesian grid of size `kspaceSize`.
Uses the CGNR method to invert the non-cartesian Fourier encoding.

# Arguments:
* `acqData::AcquisitionData`        - AcquisitionData
* `ksize::NTuple{2,Int64}`          - size of the cartesian k-space grid
* `cgnr_iter::Int64=3`              - number of CGNR iterations
* `correctionMao::Array{ComplexF64}`- relaxation/b0 map
"""
function regrid2d(acqData::AcquisitionData, kspaceSize::NTuple{2,Int64}; cgnr_iter::Int64=3, correctionMap::Array{ComplexF64}=ComplexF64[]) where D
  dcf = samplingDensity(acqData,kspaceSize)

  nx,ny = kspaceSize
  numContr, numChan, numSl = numContrasts(acqData), numChannels(acqData), numSlices(acqData)

  kdata_cart = [zeros(ComplexF64,nx*ny,numChan) for j=1:numContr, k=1:numSl, rep=1:1]
  F = sqrt(prod(kspaceSize))*FFTOp(ComplexF64, kspaceSize)
  img = zeros(ComplexF64,nx*ny)
  for k = 1:numSl
    E = encodingOps_simple(acqData, kspaceSize, slice=k, correctionMap=correctionMap)
    for j = 1:numContr
      W = WeightingOp(dcf[j])
      solver = RegularizedLeastSquares.CGNR(W*E[j], iterations=cgnr_iter)
      for i = 1:numChan
        kdata = kData(acqData,j,i,k).* dcf[j]
        img .= solve(solver, kdata)
        kdata_cart[j,k,1][:,i] .= F*img
      end
    end
  end

  return AcquisitionData(CartesianTrajectory(ny,nx), kdata_cart, encodingSize=[nx,ny,1], fov=acqData.fov)
end
