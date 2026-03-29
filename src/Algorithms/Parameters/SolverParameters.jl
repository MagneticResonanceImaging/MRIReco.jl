export SolverParameter, SolverParameters, LeastSquaresSolverParameter
export RegularizationParameters
export SimpleSparsityParameters, CustomSparsityParameters
"""
    AbstractSparsityParameters

Abstract base type for sparsity/transformation operators used in regularized reconstruction.
"""
abstract type AbstractSparsityParameters <: AbstractMRIRecoParameters end

"""
    SolverParameters <: AbstractImageReconstructionParameters

Abstract base type for solver parameters in MRIReco. Subtypes define how linear solvers
are configured for iterative reconstruction.
"""
abstract type SolverParameters <: AbstractMRIRecoParameters end

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
@parameter struct SimpleSparsityParameters{S <: Union{Nothing, String, Vector{<:Union{String, Nothing}}}} <: AbstractSparsityParameters
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
- nothing returns nothing (to be filled later)
- Single operator returns single-element vector
- String gets converted to SparseOp
- Iterable of operators returns collected vector
- Supports mixed vectors like [op, nothing, op]

# Note
- Gets reconSize from MRIRECO_CONTEXT ScopedValue
- Gets storage type S from MRIRECO_CONTEXT ScopedValue
"""
function (sparsity::CustomSparsityParameters{L})(::Type{<:AbstractIterativeMRIRecoAlgorithm}, numTerms::Int64) where L
  reconSize = ctx_reconSize()
  S = ctx_storageType()
  T = eltype(S)
  trafo = sparsity.trafo
  
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
* `reg::R` - Regularization term(s). Can be nothing, a single term, or a vector.
* `sparsity::S` - Sparsity transform specification.

# Callable Interface
```julia
(params::RegularizationParameters)(::Type{<:AbstractLinearSolver}, reconSize, T) -> reg | (reg, regTrafo)
```

Transforms the regularization for a specific solver. Returns either:
- A `regularization term` or a `TransformedRegularization` (for proximal gradient solvers like FISTA)
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
RegularizationParameters(reg, sparsity::Union{String, Nothing}) = RegularizationParameters(reg, SimpleSparsityParameters(sparsity))
RegularizationParameters(reg::Vector{<:AbstractRegularization}, sparsity::Vector{<:Union{String, Nothing}}) = RegularizationParameters(reg, SimpleSparsityParameters(sparsity))

"""
    (params::RegularizationParameters)(::Type{<:AbstractIterativeMRIRecoAlgorithm}, ::Type{<:AbstractLinearSolver}) -> reg | (reg, regTrafo)

Transform the regularization parameters for a specific solver type.

# Arguments
- `::Type{<:AbstractLinearSolver}` - The solver type (e.g., ADMM, FISTA)

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

  trafos = params.sparsity(algoT, length(reg_vec))
  return params(algoT, reg_vec, trafos, SL)
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
@parameter constructor = false struct LeastSquaresSolverParameter{SL <: AbstractLinearSolver, R <: RegularizationParameters, T <: AbstractFloat} <: SolverParameters
  solver::Type{SL} 
  regularization::R 
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
    regularization::R = RegularizationParameters(),
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
    solver, regularization, normalizeReg, iterations, rho_T,
    iterationsInner, iterationsCG, absTol_T, relTol_T, tolInner_T,
    verbose, restart, vary_rho
  )
  validate!(params)
  return params
end
"""
    (params::LeastSquaresSolverParameter)(algoT, ::Type{SL}, b::AbstractVector, A, AHA=normalOperator(A)) where SL

Solve an inverse problem using the configured solver.

# Arguments
- `algoT::Type{<:AbstractIterativeMRIRecoAlgorithm}` - Algorithm type
- `::Type{SL}` - Solver type (must match the param's solver type)
- `b::AbstractVector` - k-space data
- `A` - Forward operator
- `AHA` - Normal operator (optional, computed from A if not provided)

# Returns
Reconstructed image vector

# Example
```julia
params = LeastSquaresSolverParameter(ADMM; iterations=50)
x = params(IterativeMRIReco, ADMM, kdata, A)  # Gets reconSize from MRIRECO_CONTEXT
```
"""
function (params::LeastSquaresSolverParameter{SL, R, T})(
    algoT::Type{<:AbstractIterativeMRIRecoAlgorithm},
    ::Type{SL}, 
    b::AbstractVector, 
    A, 
    AHA = normalOperator(A)
) where {SL <: AbstractLinearSolver, R, T}
  solver = params(algoT, SL, A, AHA)
  return solve!(solver, b)
end

"""
    (params::LeastSquaresSolverParameter)(algoT, ::Type{SL}, A, AHA=normalOperator(A)) where SL

Create a configured solver for the given problem.

# Arguments
- `algoT::Type{<:AbstractIterativeMRIRecoAlgorithm}` - Algorithm type
- `::Type{SL}` - Solver type
- `A` - Forward operator
- `AHA` - Normal operator (optional, computed from A if not provided)

# Returns
A configured `AbstractLinearSolver` ready for use with `solve!`

# Note
- Gets reconSize from MRIRECO_CONTEXT ScopedValue

# Example
```julia
params = LeastSquaresSolverParameter(ADMM; iterations=50)
solver = params(IterativeMRIReco, ADMM, A, AHA)
x = solve!(solver, kdata)

# Reuse solver for different data
x2 = solve!(solver, kdata2)
```
"""
function (params::LeastSquaresSolverParameter{SL, R, T})(
    algoT::Type{<:AbstractIterativeMRIRecoAlgorithm},
    ::Type{SL}, 
    A, 
    AHA = normalOperator(A)
) where {SL <: AbstractLinearSolver, R, T}
  
  # Get regularization - returns reg or (reg, regTrafo)
  reg_result = params.regularization(algoT, SL)
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