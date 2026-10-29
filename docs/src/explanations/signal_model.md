# [Signal Model and Encoding Operators](@id signal-model)

A reconstruction needs a model of how the measured k-space data arises from the image. This page describes the signal model of MRIReco.jl, its discretization into the encoding operator and the operators that evaluate it for different trajectories, for several coils and with relaxation and off-resonance. [Reconstruction as an Inverse Problem](@ref inverse-problem) explains how a reconstruction uses this operator.

## The Signal Model

MRIReco.jl uses the signal model of equation (2) in the [MRIReco.jl paper](https://doi.org/10.1002/mrm.28792), written here with the symbols of the tutorials. The demodulated signal that coil $c$ receives at time $t$ is

```math
b_c(t) = \int S_c(\mathbf{r})\, x(\mathbf{r})\, \mathrm{e}^{-z(\mathbf{r})\, t}\, \mathrm{e}^{-2\pi\mathrm{i}\, \mathbf{k}(t) \cdot \mathbf{r}}\, \mathrm{d}\mathbf{r},
```

where $x(\mathbf{r})$ is the transverse magnetization at the position $\mathbf{r}$, which is the image to reconstruct, $S_c(\mathbf{r})$ is the sensitivity of the coil and $\mathbf{k}(t)$ is the k-space trajectory. The complex map $z(\mathbf{r}) = R_2^*(\mathbf{r}) + \mathrm{i}\,\Delta\omega(\mathbf{r})$ contains the relaxation rate $R_2^*$ in its real part and the off-resonance frequency $\Delta\omega$, in radians per second, in its imaginary part. MRIReco.jl calls it the correction map. Without the correction map and with a single coil of uniform sensitivity, the signal is the Fourier transform of the image at the k-space position $\mathbf{k}(t)$, which is why the direct reconstruction of fully sampled Cartesian data is an inverse FFT.

## The Encoding Operator

The signal is sampled at the times $t_j$, at the k-space positions $\mathbf{k}_j = \mathbf{k}(t_j)$, and the image is represented by its values $x_n$ at the positions $\mathbf{r}_n$ of a grid. The integral then becomes a sum,

```math
b_{c,j} = \sum_n S_c(\mathbf{r}_n)\, x_n\, \mathrm{e}^{-z(\mathbf{r}_n)\, t_j}\, \mathrm{e}^{-2\pi\mathrm{i}\, \mathbf{k}_j \cdot \mathbf{r}_n},
```

so the samples of all coils, stacked into the vector $\mathbf{b}$, depend linearly on the image vector $\mathbf{x}$: $\mathbf{b} = \mathbf{E}\mathbf{x}$ with the encoding operator $\mathbf{E}$.

MRIReco.jl never stores $\mathbf{E}$ as a matrix, which for an image of 256 × 256 voxels and as many samples would have more than four billion entries. Like the operators of MRIOperators and [LinearOperatorCollection.jl](https://juliaimagerecon.github.io/LinearOperatorCollection.jl/stable/), it is matrix-free: it only implements how $\mathbf{E}$ and its adjoint $\mathbf{E}^{\mathsf{H}}$ act on a vector, using fast algorithms such as the FFT. Operators are applied with `*`, their adjoint is `adjoint(E)`, and `∘` combines them, which is how MRIOperators builds $\mathbf{E}$ from the operators described in the following sections.

## Fourier Encoding

Without the correction map and the coil sensitivities, $\mathbf{E}$ is a Fourier transform at the sampled k-space positions. MRIOperators chooses the operator by the trajectory, using operators of LinearOperatorCollection.jl:

- For Cartesian trajectories, the samples lie on a grid, and the operator is an FFT, `FFTOp`. If only some of the grid points were measured, a `SamplingOp` keeps the measured samples.
- For non-Cartesian trajectories, such as radial and spiral ones, the operator is a non-uniform FFT, `NFFTOp`, computed by [NFFT.jl](https://juliamath.github.io/NFFT.jl/stable/). The NFFT approximates the sum with an FFT on a grid that is finer than the image grid by the oversampling factor σ and a convolution with a small window. The reconstruction parameters `:oversamplingFactor` and `:kernelSize`, with the defaults 1.25 and 3, set σ and the size of the window, and with them the balance between accuracy and speed.

MRIOperators also provides `ExplicitOp`, which evaluates the sum directly, including the correction map. It is exact, but its cost grows with the number of samples times the number of voxels. The reconstruction therefore uses the fast operators, and `ExplicitOp` mainly serves as a reference for testing them.

## Coil Sensitivities

With several receive coils, the image is first multiplied by the sensitivities of all coils, which `SensitivityOp` does, and the image of each coil is then Fourier encoded:

```math
\mathbf{E} = \begin{pmatrix} \mathbf{F} & & \\ & \ddots & \\ & & \mathbf{F} \end{pmatrix} \begin{pmatrix} \mathbf{S}_1 \\ \vdots \\ \mathbf{S}_C \end{pmatrix},
```

where $\mathbf{F}$ is the Fourier operator from above and $\mathbf{S}_c$ multiplies the image by the sensitivity of coil $c$. The block-diagonal operator is a `DiagOp`, which applies $\mathbf{F}$ to the image of each coil. This is the encoding operator of [SENSE](https://doi.org/10.1002/%28SICI%291522-2594%28199911%2942:5%3C952::AID-MRM16%3E3.0.CO;2-S), which the reconstruction types with coil sensitivities use.

## Relaxation and Off-Resonance

The factor $\mathrm{e}^{-z(\mathbf{r})\, t}$ depends on both position and time, so it can't be moved out of the Fourier transform. `FieldmapNFFTOp` therefore approximates it by a sum of $K$ products of a function of time and a function of position,

```math
\mathrm{e}^{-z(\mathbf{r}_n)\, t_j} \approx \sum_{\kappa=1}^{K} a_{j\kappa}\, c_{\kappa n},
```

which turns the operator into a weighted sum of $K$ NFFTs. The parameter `:method` chooses how the coefficients are computed. The default `"nfft"` approximates the exponential analytically and determines $K$ from the ranges of the correction map and of the readout times, while data-driven methods such as `"leastsquare"` and `"hist"` fit the coefficients for $K$ time segments, which the parameter `:K` sets, to the correction map or to its histogram. Section 2.2 of the paper gives references for these approximations.

The iterative reconstructions take the correction map from the parameter `:correctionMap`. They count the times from the echo time, so they reconstruct the magnetization at the echo time.

## The Normal Operator and the Toeplitz Approach

Iterative solvers such as CGNR apply the normal operator $\mathbf{E}^{\mathsf{H}}\mathbf{E}$ in every iteration, in MRIReco.jl together with the weighting described in [Reconstruction as an Inverse Problem](@ref inverse-problem). MRIReco.jl builds it with `normalOperator`, which can use a faster implementation than applying $\mathbf{E}$ and $\mathbf{E}^{\mathsf{H}}$ one after the other. For the `NFFTOp`, the normal operator is a convolution of the image with a kernel that depends only on the trajectory and the weights. With `:toeplitz = true`, MRIReco.jl computes this convolution with two FFTs on a grid of twice the image size in each dimension, so that the window convolutions of the NFFT are no longer needed in the iterations. Section 3.2 of the paper compares both approaches for a radial dataset: they were similarly fast for fully sampled data, while for fourfold undersampled data, the NFFT with σ = 1.25 was clearly faster, because its FFTs are smaller and fewer samples have to be gridded.
