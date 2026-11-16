export SolverParameter, SolverParameters, LeastSquaresSolverParameter
export RegularizationParameters
export SimpleSparsityParameters, CustomSparsityParameters
"""
    AbstractSparsityParameters

Abstract base type for sparsity/transformation operators used in regularized reconstruction.
"""
abstract type AbstractSparsityParameters <: AbstractMRIRecoParameters end

"""
    SolverParameters <: AbstractMRIRecoParameters

Abstract base type for solver parameters in MRIReco. Subtypes define how linear solvers
are configured for iterative reconstruction.
"""
abstract type SolverParameters <: AbstractMRIRecoParameters end

"""
    SimpleSparsityParameters{S} <: AbstractSparsityParameters

Parameters specifying sparsifying transforms by name, see [`SparseOp`](@ref).

# Fields
* `sparseTrafo::S` - The transform of the first regularization term (`"Wavelet"`, or `"nothing"` for the
                     identity) or a vector with one entry per term. With `""` (the default) or `nothing`, no term
                     has a transform; in a vector, `nothing` or a missing entry leaves that term without one. A term
                     without a transform applies to the image itself.

# Examples
```julia
SimpleSparsityParameters("Wavelet")               # Wavelet for the first term
SimpleSparsityParameters(["Wavelet", "nothing"])  # one transform per term
SimpleSparsityParameters()                        # no transform
```
"""
@parameter struct SimpleSparsityParameters{S <: Union{Nothing, String, Vector{<:Union{String, Nothing}}}} <: AbstractSparsityParameters
  sparseTrafo::S = ""
end

"""
    (sparsity::SimpleSparsityParameters)(algo, numTerms::Int64) -> Vector

Build the transforms of `numTerms` regularization terms: a [`SparseOp`](@ref) for each name and `nothing` for each
term without a transform. Reconstructions that solve for all contrasts or basis functions at once repeat each
transform for every contrast or basis function.

# Note
- Gets the reconstruction size and the storage type from the [`MRIRecoContext`](@ref) of the reconstruction
"""
function (sparsity::SimpleSparsityParameters)(::Type{<:AbstractMRIRecoAlgorithm}, numTerms::Int64)
  reconSize = ctx_reconSize()
  S = ctx_storageType()
  T = eltype(S)
  sp = sparsity.sparseTrafo
  
  # Normalize to vector
  if isnothing(sp) || isempty(sp)
    sp_vec = fill(nothing, numTerms)
  elseif sp isa String
    sp_vec = [sp]
  else
    sp_vec = collect(sp)
  end
  
  # Fill up with nothing to match numTerms
  if length(sp_vec) < numTerms
    sp_vec = vcat(sp_vec, fill(nothing, numTerms - length(sp_vec)))
  end
  
  # Build operators: nothing stays nothing, non-nothing strings get converted
  return [isnothing(s) ? nothing : SparseOp(T, s, reconSize; S = S) for s in sp_vec[1:numTerms]]
end

"""
    CustomSparsityParameters{L} <: AbstractSparsityParameters

Parameters specifying custom sparsifying transforms, see [Use a Custom Sparsifying Transform](@ref custom-transform).

# Fields
* `sparseTrafo::L` - The transform of the first regularization term (an `AbstractLinearOperator` such as
                     `LinearOperator(M)` for a matrix `M`, or a name as for [`SimpleSparsityParameters`](@ref)) or a
                     vector with one entry per term (default: `nothing`). With `nothing` or a missing entry, a term
                     has no transform.

# Examples
```julia
CustomSparsityParameters(op)                # op for the first term
CustomSparsityParameters([op, nothing])     # op for the first term, none for the second
CustomSparsityParameters()                  # no transform
```
"""
@parameter struct CustomSparsityParameters{L} <: AbstractSparsityParameters
  sparseTrafo::L = nothing
end

"""
    (sparsity::CustomSparsityParameters)(algo, numTerms::Int64) -> Vector

Build the transforms of `numTerms` regularization terms: the given operators, a [`SparseOp`](@ref) for each name
and `nothing` for each term without a transform. Reconstructions that solve for all contrasts or basis functions at
once repeat each transform for every contrast or basis function, so an operator has to act on a single image.

# Note
- Gets the reconstruction size and the storage type from the [`MRIRecoContext`](@ref) of the reconstruction
"""
function (sparsity::CustomSparsityParameters{L})(::Type{<:AbstractIterativeMRIRecoAlgorithm}, numTerms::Int64) where L
  reconSize = ctx_reconSize()
  S = ctx_storageType()
  T = eltype(S)
  trafo = sparsity.sparseTrafo
  
  # Normalize to vector
  trafo_vec = if isnothing(trafo)
    fill(nothing, numTerms)
  elseif trafo isa AbstractLinearOperator
    [trafo]
  elseif trafo isa String
    [trafo]
  else
    collect(trafo)
  end
  
  # Fill up with nothing to match numTerms
  if length(trafo_vec) < numTerms
    trafo_vec = vcat(trafo_vec, fill(nothing, numTerms - length(trafo_vec)))
  end
  
  # Build operators: nothing stays nothing
  result = []
  for t in trafo_vec[1:numTerms]
    if isnothing(t)
      push!(result, nothing)
    elseif t isa String
      push!(result, SparseOp(T, t, reconSize; S = S))
    elseif t isa AbstractLinearOperator
      push!(result, t)
    else
      push!(result, t)
    end
  end
  return result
end

"""
    RegularizationParameters{R, S} <: AbstractImageReconstructionParameters

Parameters specifying regularization terms and their sparsity transforms.

# Type Parameters
* `R` - Regularization type (nothing, AbstractRegularization, or Vector)
* `S` - Sparsity parameters type

# Fields
* `reg::R` - Regularization term(s). Can be `nothing`, a single term, or a vector (default: `nothing`). `nothing`
             stands for `L1Regularization(0.0)`, which does not regularize and which CGNR does not accept.
* `sparsityParams::S` - Sparsifying transforms of the terms, see [`SimpleSparsityParameters`](@ref) and
                        [`CustomSparsityParameters`](@ref) (default: `SimpleSparsityParameters()`).

# Examples
```julia
# L1 with wavelet
RegularizationParameters(reg = L1Regularization(1e-3),
                         sparsityParams = SimpleSparsityParameters("Wavelet"))

# Multiple terms
RegularizationParameters(reg = [L1Regularization(1e-3), L2Regularization(1e-4)],
                         sparsityParams = SimpleSparsityParameters(["Wavelet", "nothing"]))
```
"""
@parameter struct RegularizationParameters{
    R <: Union{Nothing, AbstractRegularization, Vector{<:AbstractRegularization}},
    S <: AbstractSparsityParameters
} <: AbstractImageReconstructionParameters
  reg::R = nothing
  sparsityParams::S = SimpleSparsityParameters()
end

"""
    (params::RegularizationParameters)(algo, solver::Type{<:AbstractLinearSolver})

Transform the regularization parameters for a specific solver type. ADMM and SplitBregman take several terms; FISTA,
POGM and OptISTA take one term; CGNR takes one `L2Regularization` and ignores its transform. FISTA, POGM, OptISTA and
CGNR also take the projections `PositiveRegularization` and `RealRegularization`; ADMM and SplitBregman currently fail
with them ([#392](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/392),
[#396](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/396)).

# Arguments
- `algo::AbstractIterativeMRIRecoAlgorithm` - The reconstruction algorithm
- `solver::Type{<:AbstractLinearSolver}` - The solver type (for example ADMM or FISTA)

# Returns
- For ADMM and SplitBregman, which take the transforms separately: the tuple `(reg, regTrafo)` with one transform
  per term, the identity of one image for a term without a transform. `"multiEcho"` and `"multiCoilMultiEcho"`
  with more than one contrast and `"multiCoilMultiEchoSubspace"` with more than one basis function therefore
  currently fail for such a term with these solvers; the transform `"nothing"` avoids this
  ([#392](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/392),
  [#396](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/396))
- For the other solvers: the vector of terms, each term with a transform as a `TransformedRegularization`

# Note
- Gets the reconstruction size and the storage type from the [`MRIRecoContext`](@ref) of the reconstruction
"""
function (params::RegularizationParameters)(
    algo::AbstractIterativeMRIRecoAlgorithm,
    ::Type{SL}, 
) where SL <: AbstractLinearSolver
  reg = params.reg
  T = real(eltype(ctx_storageType()))

  reg_vec = if isnothing(reg)
    [L1Regularization(zero(T))]
  elseif reg isa AbstractRegularization
    [reg]
  else
    reg
  end

  trafos = params.sparsityParams(algo, length(reg_vec))
  return params(algo, reg_vec, trafos, SL)
end

function (params::RegularizationParameters)(
    ::Type{<:AbstractIterativeMRIRecoAlgorithm},
    reg::Vector{<:AbstractRegularization}, 
    trafos::Vector, 
    ::Type{SL}
) where SL
  return [isnothing(t) ? r : TransformedRegularization(r, t) for (r, t) in zip(reg, trafos)]
end


function (params::RegularizationParameters)(
    ::Type{<:AbstractIterativeMRIRecoAlgorithm},
    reg::Vector{<:AbstractRegularization}, 
    trafos::Vector, 
    ::Type{SL}
  ) where SL <: Union{ADMM, SplitBregman}
  # Fill up nothing trafos with opEyes
  reconSize = ctx_reconSize()
  S = ctx_storageType()
  T = eltype(S)
  
  filled_trafos = []
  for t in trafos
    if isnothing(t)
      push!(filled_trafos, opEye(T, prod(reconSize); S = S))
    else
      push!(filled_trafos, t)
    end
  end
  
  # return reg and reg_trafo
  return (reg, filled_trafos)
end


"""
    LeastSquaresSolverParameter{SL, R, T} <: SolverParameters

Parametric solver parameters for iterative MRI reconstruction using RegularizedLeastSquares.

# Type Parameters
* `SL` - Solver type: ADMM, CGNR, FISTA, POGM, OptISTA or SplitBregman
* `R` - RegularizationParameters type
* `T` - Floating-point type of `rho` and the tolerances

# Fields
Except for `solver` and `regParams`, a field is passed to the solver as the keyword argument of the same name, but
only if it is not `nothing` and the solver takes that keyword; otherwise the solver uses its own default.
* `solver::Type{SL}` - Solver type (default: `ADMM`)
* `regParams::R` - Regularization terms and their transforms (default: `RegularizationParameters()`)
* `normalizeReg::AbstractRegularizationNormalization` - Scaling of the weights λ of the terms: `NoNormalization()`,
  `MeasurementBasedNormalization()` or `SystemMatrixBasedNormalization()` (default: `NoNormalization()`)
* `iterations::Int` - Maximum iterations; for SplitBregman, the outer iterations (default: 30)
* `rho::T` - Penalty parameter of ADMM and SplitBregman and step size of FISTA, POGM and OptISTA (default: 5e-2).
             FISTA, POGM and OptISTA always receive it, so they do not determine the step size themselves
             ([#305](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/305)).
* `iterationsInner::Union{Int, Nothing}` - Inner iterations of SplitBregman (default: `nothing`)
* `iterationsCG::Union{Int, Nothing}` - CG iterations of ADMM and SplitBregman (default: `nothing`)
* `absTol::Union{T, Nothing}` - Absolute tolerance of ADMM and SplitBregman (default: `nothing`)
* `relTol::Union{T, Nothing}` - Relative tolerance (default: `nothing`)
* `tolInner::Union{T, Nothing}` - Tolerance of the CG iterations of ADMM and SplitBregman (default: `nothing`)
* `verbose::Bool` - Print progress in all solvers except CGNR (default: `false`)
* `restart::Symbol` - Restart of FISTA and POGM: `:none` or `:gradient` (default: `:none`)
* `vary_rho::Symbol` - Adaptation of `rho` in ADMM: `:none`, `:balance` or `:PnP` (default: `:none`)

# Examples
```julia
# With L1 + Wavelet regularization
params = LeastSquaresSolverParameter(solver = FISTA,
  regParams = RegularizationParameters(
    reg = L1Regularization(1e-3),
    sparsityParams = SimpleSparsityParameters("Wavelet")
  ),
  iterations = 100
)

# Multiple regularization terms (ADMM)
params = LeastSquaresSolverParameter(solver = ADMM,
  regParams = RegularizationParameters(
    reg = [L1Regularization(1e-3), L2Regularization(1e-4)],
    sparsityParams = SimpleSparsityParameters(["Wavelet", "nothing"])
  ),
  rho = 1e-1
)
```
"""
@parameter constructor = false struct LeastSquaresSolverParameter{SL <: AbstractLinearSolver, R <: RegularizationParameters, T <: AbstractFloat} <: SolverParameters
  solver::Type{SL} 
  regParams::R 
  normalizeReg::AbstractRegularizationNormalization 
  iterations::Int 
  rho::T
  iterationsInner::Union{Int, Nothing} 
  iterationsCG::Union{Int, Nothing} 
  absTol::Union{T, Nothing} 
  relTol::Union{T, Nothing} 
  tolInner::Union{T, Nothing} 
  verbose::Bool
  restart::Symbol 
  vary_rho::Symbol
  
  @validate begin
    @assert iterations > 0 "iterations must be positive"
    @assert rho >= 0 "rho must be positive"
    @assert isnothing(iterationsInner) || iterationsInner > 0 "iterationsInner must be positive"
    @assert isnothing(iterationsCG) || iterationsCG > 0 "iterationsCG must be positive"
    @assert isnothing(absTol) || absTol > 0 "absTol must be positive"
    @assert isnothing(relTol) || relTol > 0 "relTol must be positive"
    @assert isnothing(tolInner) || tolInner > 0 "tolInner must be positive"
    @assert restart in (:none, :gradient) "restart must be :none or :gradient"
    @assert vary_rho in (:none, :balance, :PnP) "vary_rho must be :none, :balance, or :PnP"
  end
end
function LeastSquaresSolverParameter(;
    solver::Type{SL} = ADMM,
    regParams::R = RegularizationParameters(),
    normalizeReg::AbstractRegularizationNormalization = NoNormalization(),
    iterations::Int = 30,
    rho::Real = 5e-2,
    iterationsInner::Union{Int, Nothing} = nothing,
    iterationsCG::Union{Int, Nothing} = nothing,
    absTol = nothing,
    relTol = nothing,
    tolInner = nothing,
    verbose::Bool = false,
    restart::Symbol = :none,
    vary_rho::Symbol = :none
) where {SL <: AbstractLinearSolver, R <: RegularizationParameters}
  # Convert T-typed fields
  T = typeof(rho)
  rho_T = rho
  absTol_T = isnothing(absTol) ? nothing : T(absTol)
  relTol_T = isnothing(relTol) ? nothing : T(relTol)
  tolInner_T = isnothing(tolInner) ? nothing : T(tolInner)
  
  params = LeastSquaresSolverParameter{SL, R, T}(
    solver, regParams, normalizeReg, iterations, rho_T,
    iterationsInner, iterationsCG, absTol_T, relTol_T, tolInner_T,
    verbose, restart, vary_rho
  )
  validate!(params)
  return params
end
"""
    (params::LeastSquaresSolverParameter)(algo, b::AbstractVector, A, AHA=normalOperator(A))

Solve an inverse problem using the configured solver.

# Arguments
- `algo::AbstractIterativeMRIRecoAlgorithm` - The reconstruction algorithm
- `b::AbstractVector` - k-space data
- `A` - Forward operator
- `AHA` - Normal operator (optional, computed from A if not provided)

# Returns
Reconstructed image vector

# Note
- Gets the reconstruction size and the storage type from the [`MRIRecoContext`](@ref) of the reconstruction
"""
function (params::LeastSquaresSolverParameter{SL, R, T})(
    algo::AbstractIterativeMRIRecoAlgorithm,
    b::AbstractVector, 
    A, 
    AHA = normalOperator(A)
) where {SL <: AbstractLinearSolver, R, T}
  solver = params(algo, A, AHA)
  return solve!(solver, b)
end

"""
    (params::LeastSquaresSolverParameter)(algo, A, AHA=normalOperator(A))

Create a configured solver for the given problem.

# Arguments
- `algo::AbstractIterativeMRIRecoAlgorithm` - The reconstruction algorithm
- `A` - Forward operator
- `AHA` - Normal operator (optional, computed from A if not provided)

# Returns
A configured `AbstractLinearSolver` ready for use with `solve!`

# Note
- Gets the reconstruction size and the storage type from the [`MRIRecoContext`](@ref) of the reconstruction
"""
function (params::LeastSquaresSolverParameter{SL, R, T})(
    algo::AbstractIterativeMRIRecoAlgorithm,
    A, 
    AHA = normalOperator(A)
) where {SL <: AbstractLinearSolver, R, T}
  
  # Get regularization - returns reg or (reg, regTrafo)
  solver = params.solver
  reg_result = params.regParams(algo, solver)
  if reg_result isa Tuple
    reg, regTrafo = reg_result
  else
    reg = reg_result
    regTrafo = nothing
  end
  
  # Get solver's accepted kwargs
  solver_kw = union(Base.kwarg_decl.(methods(SL))...)
  
  # Build kwargs dict
  kwargs = (;
    reg = reg,
    normalizeReg = params.normalizeReg,
    iterations = params.iterations,
  )
  
  # Conditionally add params
  for (sym, val) in [
    :regTrafo => regTrafo,
    :iterationsInner => params.iterationsInner,
    :iterationsCG => params.iterationsCG,
    :absTol => params.absTol,
    :relTol => params.relTol,
    :tolInner => params.tolInner,
    :verbose => params.verbose,
    :rho => params.rho,
    :restart => params.restart,
    :vary_rho => params.vary_rho
  ]
    if !isnothing(val) && sym in solver_kw
      kwargs = merge(kwargs, (; sym => val))
    end
  end
  
  # Create solver
  return createLinearSolver(SL, A; AHA=AHA, kwargs...)
end