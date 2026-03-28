export SolverParameter, SolverParameters, LeastSquaresSolverParameter
export RegularizationParameters
export SimpleSparsityParameters, CustomSparsityParameters
"""
    AbstractSparsityParameters

Abstract base type for sparsity/transformation parameters used in regularized reconstruction.
"""
abstract type AbstractSparsityParameters end

"""
    SolverParameters <: AbstractImageReconstructionParameters

Abstract base type for solver parameters in MRIReco. Subtypes define how linear solvers
are configured for iterative reconstruction.
"""
abstract type SolverParameters <: AbstractImageReconstructionParameters end

"""
    SimpleSparsityParameters{S <: Union{Nothing, String, Vector{String}}} <: AbstractSparsityParameters

Parameters specifying sparsifying transforms via strings.

# Fields
* `sparsity::S` - String or vector of strings specifying the transform (e.g., "Wavelet", "nothing").
                  Empty string or "" defaults to identity transformation.

# Examples
```julia
SimpleSparsityParameters("Wavelet")           # single wavelet transform
SimpleSparsityParameters(["Wavelet", "nothing"])  # two transforms
SimpleSparsityParameters("")                  # identity (no transform)
```
"""
@parameter struct SimpleSparsityParameters{S <: Union{Nothing, String, Vector{String}}} <: AbstractSparsityParameters
  sparsity::S = ""
end

"""
    (sparsity::SimpleSparsityParameters)(reconSize, T) -> Vector

Build sparse transformation operators from the sparsity specification.

# Arguments
- `T::Type{<:Complex}` - Complex number type

# Returns
Vector of transformation operators

# Behavior
- Empty string or "" returns identity operator (`opEye`)
- Single string creates single SparseOp
- Vector of strings creates multiple SparseOps

# Note
- Gets reconSize from MRIRECO_CONTEXT ScopedValue
- Gets storage type S from MRIRECO_CONTEXT ScopedValue
"""
function (sparsity::SimpleSparsityParameters)(::Type{<:AbstractMRIRecoAlgorithm}, numTerms::Int64)
  reconSize = ctx_reconSize()
  S = ctx_storageType()
  T = eltype(S)
  sp = sparsity.sparsity
  # TODO fill up misisng terms with nothing
  if isempty(sp) || sp == ""
    return [opEye(T, prod(reconSize); S = S) for i in 1:numTerms]
  elseif sp isa String
    return [SparseOp(T, sp, reconSize; S = S)]
  else
    return [SparseOp(T, s, reconSize; S = S) for s in sp]
  end
end

"""
    CustomSparsityParameters{L} <: AbstractSparsityParameters

Parameters specifying custom sparsifying LinearOperators or any applicable transform.

# Fields
* `trafo::L` - A transform operator, vector of operators, or nothing.

# Examples
```julia
CustomSparsityParameters(wavelet_op)              # single custom operator
CustomSparsityParameters([op1, op2])              # multiple operators
CustomSparsityParameters(nothing)                 # identity (default)
CustomSparsityParameters(some_matrix)             # matrix as transform
```
"""
@parameter struct CustomSparsityParameters{L} <: AbstractSparsityParameters
  trafo::L = nothing
end

"""
    (sparsity::CustomSparsityParameters)(T) -> Vector

Build transformation operators from custom sparsity specification.

# Arguments
- `T::Type{<:Complex}` - Complex number type

# Returns
Vector of transformation operators

# Behavior
- nothing returns identity operator (`opEye`)
- Single operator returns single-element vector
- Iterable of operators returns collected vector

# Note
- Gets reconSize from MRIRECO_CONTEXT ScopedValue
- Gets storage type S from MRIRECO_CONTEXT ScopedValue
"""
function (sparsity::CustomSparsityParameters{L})(::Type{<:AbstractIterativeMRIRecoAlgorithm}, numTerms::Int64) where L
  reconSize = ctx_reconSize()
  S = ctx_storageType()
  trafo = sparsity.trafo
  # TODO: same as SimpleSparsityParameters, but also if element is a string fallback to SparseOp
  if trafo === nothing
    return [opEye(T, prod(reconSize); S = S)]
  elseif trafo isa AbstractLinearOperator
    return [trafo]
  else
    return collect(trafo)
  end
end

"""
    RegularizationParameters{R, S} <: AbstractImageReconstructionParameters

Parameters specifying regularization terms and their sparsity transforms.

# Type Parameters
* `R` - Regularization type (nothing, AbstractRegularization, or Vector)
* `S` - Sparsity parameters type

# Fields
* `reg::R` - Regularization term(s). Can be nothing, a single term, or a vector.
* `sparsity::S` - Sparsity transform specification.

# Callable Interface
```julia
(params::RegularizationParameters)(::Type{<:AbstractLinearSolver}, reconSize, T) -> (reg, regTrafo)
```

Transforms the regularization for a specific solver. Returns either:
- A `TransformedRegularization` (for proximal gradient solvers like FISTA)
- A tuple `(reg, regTrafo)` (for primal-dual solvers like ADMM)

# Examples
```julia
# L1 with wavelet
RegularizationParameters(reg=L1Regularization(1e-3), sparsity="Wavelet")

# Multiple terms
RegularizationParameters(reg=[L1Regularization(1e-3), L2Regularization(1e-4)], 
                         sparsity=["Wavelet", "nothing"])
```
"""
@parameter constructor = false struct RegularizationParameters{
    R <: Union{Nothing, AbstractRegularization, Vector{<:AbstractRegularization}},
    S <: AbstractSparsityParameters
} <: AbstractImageReconstructionParameters
  reg::R
  sparsity::S
end
RegularizationParameters(; reg = nothing, sparsity = "") = RegularizationParameters(reg, sparsity)
RegularizationParameters(reg::AbstractRegularization, sparsity::String) = RegularizationParameters(reg, SimpleSparsityParameters(sparsity))
RegularizationParameters(reg::Vector{<:AbstractRegularization}, sparsity::String) = RegularizationParameters(reg, [sparsity])
RegularizationParameters(reg::Vector{<:AbstractRegularization}, sparsity::Vector{String}) = RegularizationParameters(reg, SimpleSparsityParameters(sparsity))

"""
    (params::RegularizationParameters)(::Type{<:AbstractIterativeMRIRecoAlgorithm}, ::Type{<:AbstractLinearSolver}, T) -> (reg, regTrafo)

Transform the regularization parameters for a specific solver type.

# Arguments
- `::Type{<:AbstractLinearSolver}` - The solver type (e.g., ADMM, FISTA)
- `T::Type{<:Complex}` - Complex number type

# Returns
- For proximal gradient solvers (FISTA, POGM, OptISTA): `TransformedRegularization` or raw regularization
- For primal-dual solvers (ADMM, SplitBregman): Tuple `(reg, regTrafo)`

# Throws
- Error if proximal gradient solver is used without regularization
- Error if multiple regularization terms are used with proximal gradient solvers

# Note
- Gets reconSize from MRIRECO_CONTEXT ScopedValue
"""
function (params::RegularizationParameters)(
    algoT::Type{<:AbstractIterativeMRIRecoAlgorithm},
    ::Type{SL}, 
    T::Type{<:Complex}
) where SL <: AbstractLinearSolver
  reg = params.reg

  reg_vec = if isnothing(reg)
    [L1Regularization(zero(real(T)))]
  elseif reg isa AbstractRegularization
    [reg]
  else
    reg
  end

  trafos = params.sparsity(algo.T, length(reg))
  return params(algoT, reg_vec, trafos, SL)
end

function (params::RegularizationParameters)(
    algoT::Type{<:AbstractIterativeMRIRecoAlgorithm},
    reg::Vector{<:AbstractRegularization}, 
    trafos::Vector, 
    ::Type{SL}
) where SL <: Union{CGNR, FISTA, POGM, OptISTA}
  if length(reg) != 1
    error("$(SL) supports only a single regularization term, got $(length(reg))")
  end
  r = reg[1]
  t = trafos[1]
  return isnothing(t) ? r : TransformedRegularization(r, t)
end


function (params::RegularizationParameters)(
    ::Type{<:AbstractIterativeMRIRecoAlgorithm},
    reg::Vector{<:AbstractRegularization}, 
    trafos::Vector, 
    ::Type{SL}
) where SL <: Union{ADMM, SplitBregman}
  # TODO: fill up nothing trafos with opEyes
  return (reg, trafos)
end


"""
    LeastSquaresSolverParameters{SL, R} <: SolverParameters

Parametric solver parameters for iterative MRI reconstruction using RegularizedLeastSquares.

# Type Parameters
* `SL` - Solver type (ADMM, FISTA, SplitBregman, POGM, OptISTA, CGNR, etc.)
* `R` - RegularizationParameters type

# Fields
* `solver::Type{SL}`       - Solver type (default: ADMM)
* `regularization::R`      - Regularization configuration
* `normalizeReg`           - Regularization normalization (default: NoNormalization())
* `iterations::Int`        - Maximum iterations (default: 30)
* `iterationsInner::Int`   - Inner iterations for ADMM/SplitBregman (default: 10)
* `iterationsCG::Int`      - CG iterations for ADMM/SplitBregman (default: 10)
* `rho::Real`              - Penalty parameter (default: 5e-2)
* `absTol::Real`           - Absolute tolerance (default: eps(Float64))
* `relTol::Real`           - Relative tolerance (default: eps(Float64))
* `tolInner::Real`         - Inner CG tolerance (default: 1e-5)
* `verbose::Bool`          - Print progress (default: false)
* `restart::Symbol`        - FISTA restart strategy: :none or :gradient (default: :none)
* `vary_rho::Symbol`       - ADMM rho adaptation: :none, :balance, :PnP (default: :none)

# Callable Interface

    (params::LeastSquaresSolverParameter)(solverType, b, A, AHA)
    (params::LeastSquaresSolverParameter)(solverType, A, AHA)

Solve an inverse problem or get a configured solver. Uses MRIRECO_CONTEXT ScopedValue
for reconSize and storage type.

    solver = params(ADMM, A, AHA)
    x = solve!(solver, b)
    # or directly:
    x = params(ADMM, b, A, AHA)

# Examples
```julia
# Inside MRIRECO_CONTEXT
# With L1 + Wavelet regularization
params = LeastSquaresSolverParameter(solver = FISTA
  regularization = RegularizationParameters(
    reg = L1Regularization(1e-3),
    sparsity = "Wavelet"
  ),
  iterations = 100
)

# Multiple regularization terms (ADMM)
params = LeastSquaresSolverParameter(solver = ADMM
  regularization = RegularizationParameters(
    reg = [L1Regularization(1e-3), L2Regularization(1e-4)],
    sparsity = ["Wavelet", "nothing"]
  ),
  rho = 1e-1
)
```
"""
@parameter struct LeastSquaresSolverParameter{SL <: AbstractLinearSolver, R <: RegularizationParameters} <: LeastSquaresSolverParameter
  solver::Type{SL} = ADMM
  regularization::R = RegularizationParameters()
  normalizeReg::AbstractRegularizationNormalization = NoNormalization()
  
  iterations::Int = 30
  iterationsInner::Int = 10
  iterationsCG::Int = 10
  rho::Real = 5e-2
  absTol::Real = eps(Float64)
  relTol::Real = eps(Float64)
  tolInner::Real = 1e-5
  verbose::Bool = false
  restart::Symbol = :none
  vary_rho::Symbol = :none
  
  @validate begin
    @assert iterations > 0 "iterations must be positive"
    @assert rho >= 0 "rho must be positive"
    @assert restart in (:none, :gradient) "restart must be :none or :gradient"
    @assert vary_rho in (:none, :balance, :PnP) "vary_rho must be :none, :balance, or :PnP"
  end
end

"""
    (params::LeastSquaresSolverParameter)(::Type{SL}, b::AbstractVector, A, AHA=normalOperator(A)) where SL

Solve an inverse problem using the configured solver.

# Arguments
- `::Type{SL}` - Solver type (must match the param's solver type)
- `b::AbstractVector` - k-space data
- `A` - Forward operator
- `AHA` - Normal operator (optional, computed from A if not provided)

# Returns
Reconstructed image vector

# Throws
Error if reconSize is not provided. Use the signature with reconSize instead.

# Example
```julia
params = LeastSquaresSolverParameter(ADMM; iterations=50)
x = params(ADMM, kdata, A)  # Gets reconSize from MRIRECO_CONTEXT
```
"""
function (params::LeastSquaresSolverParameter{SL, R})(
    algo::Type{<:AbstractIterativeMRIRecoAlgorithm},
    b::AbstractVector, 
    A, 
    AHA = normalOperator(A)
) where {SL <: AbstractLinearSolver, R}
  solver = params(A, AHA)
  return solve!(solver, b)
end

"""
    (params::LeastSquaresSolverParameter)(::Type{SL}, A, AHA=normalOperator(A)) where SL

Create a configured solver for the given problem.

# Arguments
- `::Type{SL}` - Solver type (must match the param's solver type)
- `A` - Forward operator
- `AHA` - Normal operator (optional, computed from A if not provided)

# Returns
A configured `AbstractLinearSolver` ready for use with `solve!`

# Note
- Gets reconSize from MRIRECO_CONTEXT ScopedValue

# Example
```julia
params = LeastSquaresSolverParameter(ADMM; iterations=50)
solver = params(ADMM, A, AHA)
x = solve!(solver, kdata)

# Reuse solver for different data
x2 = solve!(solver, kdata2)
```
"""
function (params::LeastSquaresSolverParameter{SL, R})(
    ::Type{SL}, 
    A, 
    AHA = normalOperator(A)
) where {SL <: AbstractLinearSolver, R}
  return solver_from_params(params, A, AHA)
end

"""
    solver_from_params(params, A, AHA) -> AbstractLinearSolver

Create a solver from parameters.

# Arguments
- `params::LeastSquaresSolverParameter` - Configuration parameters
- `A` - Forward operator
- `AHA` - Normal operator

# Note
- Gets reconSize from MRIRECO_CONTEXT ScopedValue

# Returns
Configured solver instance
"""
function solver_from_params(params, A, AHA)
  T = eltype(AHA)
  reg, regTrafo = params.regularization(params.solver, Complex{T})
  kwargs = build_solver_kwargs(params, reg, regTrafo)
  return createLinearSolver(params.solver, A; AHA=AHA, kwargs...)
end

"""
    build_solver_kwargs(params, reg, regTrafo) -> NamedTuple

Build common keyword arguments for solver configuration.

# Arguments
- `params::LeastSquaresSolverParameter` - Configuration parameters
- `reg` - Regularization terms
- `regTrafo` - Regularization transformations

# Returns
NamedTuple with common solver kwargs
"""
function build_solver_kwargs(params, reg, regTrafo)
  (;
    reg = reg,
    regTrafo = regTrafo,
    normalizeReg = params.normalizeReg,
    iterations = params.iterations,
    absTol = params.absTol,
    relTol = params.relTol,
    verbose = params.verbose
  )
end

function build_solver_kwargs(params::LeastSquaresSolverParameter{FISTA, R}, reg, regTrafo) where R
  merge(build_solver_kwargs(params, reg, regTrafo), (;
    restart = params.restart,
    rho = params.rho
  ))
end

function build_solver_kwargs(params::LeastSquaresSolverParameter{ADMM, R}, reg, regTrafo) where R
  merge(build_solver_kwargs(params, reg, regTrafo), (;
    iterationsInner = params.iterationsInner,
    iterationsCG = params.iterationsCG,
    tolInner = params.tolInner,
    rho = params.rho,
    vary_rho = params.vary_rho
  ))
end

function build_solver_kwargs(params::LeastSquaresSolverParameter{SplitBregman, R}, reg, regTrafo) where R
  merge(build_solver_kwargs(params, reg, regTrafo), (;
    iterationsInner = params.iterationsInner,
    iterationsCG = params.iterationsCG,
    tolInner = params.tolInner,
    rho = params.rho
  ))
end

function build_solver_kwargs(params::LeastSquaresSolverParameter{POGM, R}, reg, regTrafo) where R
  merge(build_solver_kwargs(params, reg, regTrafo), (;
    restart = params.restart,
    rho = params.rho
  ))
end

function build_solver_kwargs(params::LeastSquaresSolverParameter{OptISTA, R}, reg, regTrafo) where R
  merge(build_solver_kwargs(params, reg, regTrafo), (;
    restart = params.restart,
    rho = params.rho
  ))
end

function build_solver_kwargs(params::LeastSquaresSolverParameter{CGNR, R}, reg, regTrafo) where R
  build_solver_kwargs(params, reg, regTrafo)
end