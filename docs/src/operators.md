# Imaging Operators

The mapping between the proton density and the recorded signal is linear in MRI
and can be described in the continuous case as an integral equation and in the
discrete case as a matrix vector multiplication.

Depending on the imaging scenario, the MRI system matrix can have various different
forms. It may encode a Cartesian, or a spiral trajectory. It may take offresonance
into account, and it may also encode the sensitivity of the receive coil.

MRIReco implements various MRI imaging operators. In all cases, the operators
have a dedicated Julia type that acts as a matrix. The operator `E` thus can be
applied to a vector `x` by calling `E*x`. Similarly, the adjoint can be applied by
`adjoint(E)*x`. We note at this point that the adjoint operation is lazy in Julia
and thus the matrix `adjoint(E)` is never explicitly arranged.

MRIReco currently implements the following operators:
* `FFTOp`: A multidimensional FFT operator
* `NFFTOp`: A multidimensional operator for non-equidistant FFTs
* `FieldmapNFFTOp`: An operator that takes complex correction terms into account
* `SensitivityMapOp`: An operator for building a SENSE reconstruction. Has to be
  combined with one of the former encoding operators
* `SamplingOp`: An operator that describes the (sub)sampling of full trajectories.
  The operator is used for Compressed Sensing reconstruction
* `WaveletOp`: A multidimensional operator for applying Wavelet transformations

Each of these operators can be build by calling the corresponding constructor.
Alternatively one can use the `EncodingOp` constructor that allows for high-level
construction of the imaging operator.
