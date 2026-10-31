# [Reconstruction as an Inverse Problem](@id inverse-problem)

The function [`reconstruction`](@ref) computes an image from the measured data with the encoding operator described in [Signal Model and Encoding Operators](@ref signal-model). This page explains the problem it solves, the building blocks it uses and how the reconstruction types differ.

## The Inverse Problem

With the encoding operator $\mathbf{E}$, the measured samples $\mathbf{b}$ and the unknown image $\mathbf{x}$ satisfy $\mathbf{b} \approx \mathbf{E}\mathbf{x}$, up to noise and model errors. Computing $\mathbf{x}$ from $\mathbf{b}$ is an inverse problem. It is usually ill-posed: for undersampled data, many images fit the data equally well, and noise in the data can be strongly amplified. The iterative reconstructions therefore solve a regularized least-squares problem, a weighted form of equation (1) in the [MRIReco.jl paper](https://doi.org/10.1002/mrm.28792):

```math
\underset{\mathbf{x}}{\operatorname{argmin}} \; \frac{1}{2}\lVert \mathbf{W}(\mathbf{E}\mathbf{x} - \mathbf{b}) \rVert_2^2 + \mathrm{R}(\mathbf{x}).
```

The regularization term $\mathrm{R}$ expresses prior knowledge about the image, for example that it consists of regions with constant intensity or has few large wavelet coefficients, and the diagonal matrix $\mathbf{W}$ weights the samples, as described below. [Iterative Reconstruction](@ref iterative-reconstruction) shows how the choice of the regularization and of its weight λ affects the image.

The direct reconstruction, `"direct"`, doesn't solve this problem. It applies the adjoint of the encoding operator to the weighted data, $\mathbf{x} = \mathbf{E}^{\mathsf{H}}\mathbf{W}^2\mathbf{b}$. For fully sampled Cartesian data, this is the inverse FFT, and for non-Cartesian data, it is the gridding reconstruction.

## Density Compensation

Non-Cartesian trajectories sample k-space unevenly, radial trajectories for example much more densely at the center than at the edges. The density compensation factors on the diagonal of $\mathbf{W}^2$ compensate for this in the direct reconstruction. For non-Cartesian data, MRIReco.jl computes them with the [iterative method of Pipe and Menon](https://doi.org/10.1002/%28SICI%291522-2594%28199901%2941:1%3C179::AID-MRM25%3E3.0.CO;2-V), and for Cartesian data, they are constant. [`samplingDensity`](@ref) returns their square roots, the diagonal of $\mathbf{W}$. The iterative reconstructions weight the data term with $\mathbf{W}$ as well, which [speeds up the convergence](https://doi.org/10.1002/mrm.1241) of the solvers. With `:densityWeighting = false`, they use a constant $\mathbf{W}$ instead.

## Solvers and Regularization

MRIReco.jl takes the solvers and the regularization terms from [RegularizedLeastSquares.jl](https://juliaimagerecon.github.io/RegularizedLeastSquares.jl/stable/), as described in section 2.5 of the paper:

- The parameter `:solver` selects the solver, for example `CGNR`, the conjugate gradient method applied to the normal equations, which supports an $\ell_2$ regularization, or `ADMM` and `FISTA`, which also support non-smooth terms such as the $\ell_1$ norm and total variation. The solvers only need the encoding operator, its adjoint and the normal operator.
- The parameter `:reg` sets the regularization term, for example `L2Regularization(λ)`, `L1Regularization(λ)` or `TVRegularization(λ; shape)`. Each term provides its proximal map, which proximal solvers such as ADMM and FISTA use.
- The parameter `:sparseTrafo` applies the regularization to a transformed image. With `"Wavelet"`, an $\ell_1$ regularization becomes $\mathrm{R}(\mathbf{x}) = \lambda \lVert \mathbf{\Phi}\mathbf{x} \rVert_1$ with the wavelet transform $\mathbf{\Phi}$.
- The parameter `:normalizeReg` scales λ. By default, λ is used as given, and with `MeasurementBasedNormalization`, it is scaled with the magnitude of the data, so that the same value works for data of different scale.

The examples in the paper select the solver and the regularization by name, with `params[:solver] = "admm"` and `params[:regularization] = "TV"` together with `params[:λ]`. The current version takes them as types and objects, `params[:solver] = ADMM` and `params[:reg] = TVRegularization(λ; shape)`, and the other parameters of the solvers have the names of their keyword arguments in RegularizedLeastSquares.jl, for example `:rho` instead of `:ρ` for the penalty parameter of ADMM and `:iterationsCG` instead of `:iterationsInner` for its inner CG iterations. The old names, such as `:regularization`, `:λ` and `:ρ`, are ignored.

## Reconstruction Types

The reconstruction type, `:reco`, decides which data is reconstructed together and with which encoding operator. Different slices and repetitions are always reconstructed separately.

| `:reco` | Reconstructed together | Encoding operator |
|:---|:---|:---|
| `"direct"` | the data of one coil and contrast | Fourier encoding, only its adjoint is applied |
| `"standard"` | the data of one coil and contrast | Fourier encoding |
| `"multiEcho"` | all contrasts of one coil | Fourier encoding of each contrast |
| `"multiCoil"` | all coils of one contrast | coil sensitivities and Fourier encoding |
| `"multiCoilMultiEcho"` | all coils and contrasts | coil sensitivities and Fourier encoding of each contrast |
| `"multiCoilMultiEchoSubspace"` | all coils and contrasts | as for `"multiCoilMultiEcho"`, applied to the images of a temporal basis |

The reconstruction types with coil sensitivities, which `:senseMaps` provides, combine the data of all coils into one image. The subspace reconstruction represents the images of all contrasts by a few coefficient images in a temporal basis, which `:basis` provides, and reconstructs these coefficient images, as shown in [Subspace Reconstruction](@ref subspace-reconstruction).

## Noise Decorrelation

The noise of different receive coils is usually correlated. If a noise measurement is given with `:noiseData`, the reconstruction types `"multiCoil"` and `"multiCoilMultiEcho"` estimate the noise covariance matrix $\mathbf{\Psi}$ from it, compute its Cholesky factorization $\mathbf{\Psi} = \mathbf{L}\mathbf{L}^{\mathsf{H}}$ and apply $\mathbf{L}^{-1}$ to the data and to the coil sensitivities. This decorrelates the noise, which the least-squares problem assumes to be uncorrelated and of equal variance in all coils.

The same problem can also be solved on a GPU, as described in [Run a Reconstruction on the GPU](@ref gpu-reconstruction).
