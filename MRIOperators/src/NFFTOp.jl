
function LinearOperatorCollection.NFFTOp( tr::Trajectory{T}; shape::Tuple, toeplitz=false, oversamplingFactor=1.25, kernelSize=3, kargs...) where T
  return LinearOperatorCollection.NFFTOp( Complex{T}; shape, nodes=kspaceNodes(tr), toeplitz=toeplitz, oversamplingFactor=oversamplingFactor, kernelSize=kernelSize, kargs...)
end