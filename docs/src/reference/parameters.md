# [Reconstruction Parameters](@id reconstruction-parameters)

[`reconstruction`](@ref) takes the parameters of a reconstruction as a dictionary with `Symbol` keys, for example `params[:iterations] = 20` for `params = Dict{Symbol, Any}()`. This page lists the keys, their default values and the reconstruction types that use them. Keys that the chosen reconstruction type does not use are ignored, and so are misspelled keys and names that older versions used, such as `:λ`, `:ρ`, `:regularization`, `:regTrafo`, which `:sparseTrafo` replaces, `:sparseTrafoName` and `:kwargWarning`. This happens without a warning, except for keys whose value is a parameter object, an algorithm or a plan, and for the four keys of the next paragraph with `"direct"`, which warns that it cannot set the parameter they select; with `ENV["JULIA_DEBUG"] = "MRIReco"`, `reconstruction` lists the other keys that the plan does not have, `:densityWeighting` included. A value that a parameter can't take, for example a solver given by name as in older versions, `:solver => "fista"`, is ignored with a warning, and the parameter keeps its default; this does not hold for the four keys of the next paragraph.

Apart from `:reco`, which selects the reconstruction plan, the keys set the parameters of the same name in the plan, see [Configure a Reconstruction with Plans](@ref configure-plans). For `:densityWeighting = false`, the `:weights` of the iterative reconstruction types, `:encodingOps` and transform operators in `:sparseTrafo`, `reconstruction` first replaces the parameter type of the weighting, the encoding or the transform in the plan. These parameter types can also be given directly, for example `:weightingParams => UniformWeightingParameters()` or `:encodingParams => CustomEncodingParameters(encodingOps = ops)`. [`reconstruct`](@ref) takes the same names as keywords and has the same defaults, for example `reconstruct("standard", acqData; solver = CGNR, reg = L2Regularization(1e-3))`, but it needs the parameter types: it ignores the four keys with a warning, as it does every keyword that the plan doesn't have.

The plan of a reconstruction type shows the parameters that the keys set, for example for `"multiCoil"`:

```@example parameters
using MRIReco
MRIRecoPlan("multiCoil")
```

## Reconstruction Type and Image

- `:reco` – the reconstruction type: `"direct"`, the default, `"standard"`, `"multiEcho"`, `"multiCoil"`, `"multiCoilMultiEcho"` or `"multiCoilMultiEchoSubspace"`, see [Reconstruction Types](../explanations/inverse_problem.md#Reconstruction-Types). It also takes the path of a saved plan, or its name if its folder was added with [`addMRIRecoPlanPath`](@ref).
- `:reconSize` – the size of the image as a tuple, by default the encoding size of the data, `encodingSize(acqData)`.
- `:arrayType` – the type of the arrays that the reconstruction computes with, by default `Array`. A GPU array type, such as `CuArray`, runs the reconstruction on the GPU, see [Run a Reconstruction on the GPU](@ref gpu-reconstruction).
- `:scheduler` – iterative types: how the independent problems of a reconstruction, for example of different slices, are distributed over the Julia threads. It takes `:dynamic`, `:static`, `:greedy`, `:serial` or a scheduler of [OhMyThreads.jl](https://juliafolds2.github.io/OhMyThreads.jl/stable/). The default is `:dynamic` on the CPU and `:serial` on the GPU.

## Encoding Operator

The iterative reconstruction types build their encoding operator from the data and these keys, see [Signal Model and Encoding Operators](@ref signal-model):

- `:correctionMap` – the complex correction map for off-resonance and relaxation, in the precision of the data, by default none, see [Correct Off-Resonance](@ref correct-off-resonance).
- `:method` – how a correction map is approximated: `"nfft"`, the default, `"leastsquare"`, `"hist"` or `"svd"`.
- `:K` – the number of terms of `"leastsquare"` and `"hist"` and the maximum number of terms of `"svd"`, by default 20. `"nfft"` determines the number itself.
- `:K_tol` – the tolerance with which `"svd"` chooses the number of terms, by default `1e-3`.
- `:toeplitz` – whether the normal operator uses the Toeplitz approach, by default `false`, see [The Normal Operator and the Toeplitz Approach](@ref).
- `:oversamplingFactor` and `:kernelSize` – the oversampling factor and the kernel size of the NFFT, by default `1.25` and `3`.
- `:encodingOps` – `"standard"` and `"multiCoil"`: encoding operators to use instead of the ones built from the data, as a matrix with one row per contrast and one column per slice, for example `reshape(encodingOps_simple(acqData, reconSize), :, 1)` for one slice; `"multiCoil"` takes the operators of [`encodingOps_parallel`](@ref) instead and still needs `:senseMaps` or `:noiseData`. The keys above then have no effect, and with `:noiseData`, `"multiCoil"` decorrelates the data but not these operators. [Use Encoding Operators Directly](@ref encoding-operators) describes the operators.

`:toeplitz`, `:oversamplingFactor` and `:kernelSize` only take effect for non-Cartesian data without a correction map. Without a correction map, Cartesian data uses the FFT; with one, the operator always uses the defaults. The direct reconstruction takes its correction map from `:cmap` instead, always with the approximation `"nfft"` ([#296](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/296)).

## Weighting

- `:densityWeighting` – iterative types: whether the data term is weighted with the square roots of the density compensation factors, by default `true`, or with a constant, see [Density Compensation](@ref).
- `:weights` – iterative types: weights to use instead of the density compensation, a vector with one vector per contrast and one weight per sample, in the form that [`samplingDensity`](@ref) returns.

## Coil Sensitivities and Noise

The reconstruction types `"multiCoil"`, `"multiCoilMultiEcho"` and `"multiCoilMultiEchoSubspace"` take:

- `:senseMaps` – the coil sensitivities, a complex array with the dimensions x, y, z or slice, and coil, in the precision of the data, as [`espirit`](@ref) returns them when called on an `AcquisitionData` with `nmaps = 1`. It is required unless `:encodingOps` and `:noiseData` are given, see [Multi-Coil Reconstruction](@ref multicoil-reconstruction).
- `:noiseData` – noise samples as a complex matrix of samples × coils in the precision of the data, by default none, see [Noise Decorrelation](@ref).

## Regularization

- `:reg` – a regularization term of [RegularizedLeastSquares.jl](https://juliaimagerecon.github.io/RegularizedLeastSquares.jl/stable/), for example `L2Regularization(λ)`, `L1Regularization(λ)` or `TVRegularization(λ; shape)`, or a vector of terms. ADMM and SplitBregman take several terms; FISTA, POGM and OptISTA take one term and CGNR an `L2Regularization`, and these four also take the projections `PositiveRegularization` and `RealRegularization`. The default `L1Regularization(0.0)` doesn't regularize; CGNR therefore needs this key.
- `:sparseTrafo` – the transform in whose domain a term applies, by default none: `"Wavelet"`, `"nothing"` for the identity, or an operator, for example `LinearOperator(M)` for a matrix `M`, see [Use a Custom Sparsifying Transform](@ref custom-transform). A single transform belongs to the first term; for several terms, a vector gives one transform per term, with `nothing` for a term without one. CGNR ignores the transform of its term.
- `:normalizeReg` – how the weight λ of the terms is scaled: `NoNormalization()`, the default, `MeasurementBasedNormalization()` or `SystemMatrixBasedNormalization()`.

## Solver

- `:solver` – the solver of RegularizedLeastSquares.jl: `ADMM`, the default, `CGNR`, `FISTA`, `POGM`, `OptISTA` or `SplitBregman`. `Kaczmarz`, `DirectSolver` and `PseudoInverse` don't work with the encoding operators.
- `:iterations` – the maximum number of (outer) iterations, by default 30.
- `:rho` – by default `0.05`: the penalty parameter of ADMM and SplitBregman and the step size of FISTA, POGM and OptISTA. These three would otherwise determine the step size themselves, but the reconstruction always passes `:rho` ([#305](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/305)).
- `:iterationsCG` – ADMM and SplitBregman: the maximum number of inner conjugate gradient iterations.
- `:iterationsInner` – SplitBregman: the maximum number of inner iterations.
- `:relTol` – the relative tolerance of the stopping criterion of all solvers above; `:absTol` – the absolute tolerance of ADMM and SplitBregman; `:tolInner` – the tolerance of their inner conjugate gradient iterations.
- `:restart` – FISTA and POGM: `:none`, the default, or `:gradient`.
- `:vary_rho` – ADMM: `:none`, the default, `:balance` or `:PnP`.
- `:verbose` – all solvers above except CGNR: whether the solver prints its progress, by default `false`.

Without a value, the tolerances and inner iterations keep the defaults of RegularizedLeastSquares.jl. A solver ignores the keys that it doesn't take, and other keyword arguments of the solvers, such as `precon` and `theta`, can't be set.

## Subspace Reconstruction

- `:basis` – required for `"multiCoilMultiEchoSubspace"`: the temporal basis as a complex matrix with one row per contrast and one column per basis function, in the precision of the data, see [Multi-Coil Multi-Echo and Subspace Reconstruction](@ref subspace-reconstruction).

## Limits

- With ADMM or SplitBregman, `"multiEcho"` and `"multiCoilMultiEcho"` with more than one contrast and `"multiCoilMultiEchoSubspace"` with more than one basis function fail for a regularization term without a transform; the transform `"nothing"`, for example `:sparseTrafo => "nothing"`, avoids this. Both solvers also fail for projection terms such as `PositiveRegularization` ([#392](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/392), [#396](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/396)).
- For data with several coils, `:densityWeighting = false` and `:weights` make the iterative reconstructions fail with a `DimensionMismatch` ([#392](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/392), [#396](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/396)).
- `:encodingOps` fails for `"multiEcho"`, `"multiCoilMultiEcho"` and `"multiCoilMultiEchoSubspace"` ([#396](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/396)).
- If the correction map of a slice is spatially constant, for example zero while the other slices have a map, the default `:method => "nfft"` gives `NaN` for this slice and the other methods throw an error, except `"leastsquare"` for a constant other than zero; a map that is zero in every slice is not used. If all samples have the same time, as for `AcquisitionData` created from a k-space array, `"nfft"` gives `NaN`, `"leastsquare"` and `"hist"` throw an error, and `"svd"` does so in some runs ([#360](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/360)).
- `"hist"` fails for maps with both a relaxation and an off-resonance part and for many maps with only a relaxation part, and `"svd"` for maps whose off-resonance part is spatially constant, zero included ([#360](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/360)). `"svd"` also fails for `Float32` data and is wrong without an error if the profiles have different sampling times or are too short ([#363](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/363)).
- The correction map of 3D data makes the reconstruction fail with a `MethodError` ([#361](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/361)).
