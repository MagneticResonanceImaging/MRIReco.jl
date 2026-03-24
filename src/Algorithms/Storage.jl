export addRecoPlanPath, getRecoPlanList
const DEFAULT_PLANS_PATH = joinpath(@__DIR__, "..", "config")
const recoPlanPaths = AbstractString[DEFAULT_PLANS_PATH]
"""
    addRecoPlanPath(path)

Add all `RecoPlans` within the given directory as potential 
"""
addRecoPlanPath(path) = !(path in recoPlanPaths) ? pushfirst!(recoPlanPaths, path) : nothing
"""
    getRecoPlanList(; full = false)

Retrieve a list of currently known `RecoPlan`s. If `full` is `true` then full paths are returned. 
"""
function getRecoPlanList(; full = false)
  result = String[]
  for path in recoPlanPaths
    if isdir(path)
      push!(result, [plan for plan in filter(a->contains(a,".toml"),readdir(path, join =  full))]...)
    end
  end
  return result
end

export addRecoPlanModule, getRecoPlanModules
const recoPlanModules::Vector{Module} = [AbstractImageReconstruction, MRIOperators, MRIBase, MRIReco, RegularizedLeastSquares]
"""
    addRecoPlanPath(mod::Module)

Add the `mod` module to the list of modules which are used for plan loading 
"""
addRecoPlanModule(mod::Module) = !(mod in recoPlanModules) ? push!(recoPlanModules, mod) : nothing
"""
    getRecoPlanModules()

Retrieve a list of currently used modules for plan loading 
"""
getRecoPlanModules() = recoPlanModules

function planpath(name::AbstractString)
  if isfile(name)
    return name
  end

  for dir in recoPlanPaths
    filename = joinpath(dir, string(name, ".toml"))
    if isfile(filename)
      return filename
    end
  end
  throw(ArgumentError("Could not find a suitable MRI reconstruction plan with name $name. Custom plans can be stored in the following directories $(join(recoPlanPaths, ", "))."))
end

export reconstruct
"""
    reconstruct(name::AbstractString, data::MRIFile, cache::Bool = false, modules = getRecoPlanModules(); kwargs...)

Perform a reconstruction with the `RecoPlan` specified by `name` and given `data`.
Additional keyword arguments can be passed to the reconstruction plan.

`RecoPlans` can be stored in the in the MRIReco package config folder or in a custom folder. New folder can be added with `addRecoPlanPath`. The first plan found is used.
Alternatively, name can be a path to specific plan file.

If `cache` is `true` the reconstruction plan is cached and reused if the plan file has not changed. If a keyword argument changes the structure of the plan the cache is also bypassed.
The cache considers the last modification time of the plan file and can be manually be emptied with `emptyRecoCache!()`.
The cache size (in number of plans) can be changed with by setting the `MRIRECO_CACHE_SIZE` environment variable.

# Examples
```julia
# TODO:
```
"""
function AbstractImageReconstruction.reconstruct(name::AbstractString, acqData::AcquisitionData; cache::Bool = false, modules = getRecoPlanModules(), kwargs...)
  plan = loadRecoPlan(name, cache, modules; kwargs...)
  return AbstractImageReconstruction.reconstruct(build(plan), data)
end
"""
    reconstruct(f::Base.Callable, args...; kwargs...)

Attach a callback `f(solver, frame, iteration)` to the reconstruction. Callbacks are passed as a `callbacks` keyword argument to the underlying algorithm
For more information on solver callbacks, see the RegularizedLeastSquares documentation.
"""
function AbstractImageReconstruction.reconstruct(f, name::AbstractString, args...; kwargs...)
  frame = 0
  function frame_counter_callback(solver, iteration)
    if iteration == 0
      frame += 1
    end
    f(solver, frame, iteration)
  end
  reconstruct(name, args...; kwargs..., callbacks = frame_counter_callback)
end
# Load plan with RecoCache consideration
function loadRecoPlan(name::AbstractString, cache::Bool, modules; kwargs...)
  planfile = planpath(name)

  plan = loadRecoPlan(planfile, modules; kwargs...)
  return plan
end
function setKwargs!(plan; kwargs...)
  setFirst = filter(kw->kw[2] isa Union{<:AbstractImageReconstructionAlgorithm, <:AbstractImageReconstructionParameters, <:AbstractRecoPlan}, kwargs)
  setAll!(plan; setFirst...)
  setAll!(plan; [kw for kw in kwargs if kw ∉ setFirst]...)
end
# Load plan from a .toml file
function loadRecoPlan(planfile::AbstractString, modules; kwargs...)
  return open(planfile, "r") do io
    return loadRecoPlan(io, modules; kwargs...)
  end
end
# Load plan from an io (could be file or iobuffer backed string)
function loadRecoPlan(io, modules; kwargs...)
  plan = loadPlan(io, modules; field_style = MRIRecoStyle())
  setKwargs!(plan; kwargs...)
  return plan
end
export MRIRecoPlan
"""
    MRIRecoPlan(value, modules = getRecoPlanModules(); kwargs...)

Load a `RecoPlan` and set the keyword arguments if applicable. Value can be the name of a plan or path to plan file or an MDF or a path to one.
"""
function MRIRecoPlan(value::AbstractString, modules = getRecoPlanModules(); kwargs...)
  if isfile(value) && endswith(value, ".toml")
    return loadRecoPlan(value, modules; kwargs...)
  elseif isfile(value)
    return MRIRecoPlan(MDFFile(value), modules; kwargs...)
  else
    return loadRecoPlan(planpath(value), modules; kwargs...)
  end
end

"""
    emptyRecoCache!(plan::RecoPlan)

Empty all caches within the `RecoPlan`
"""
function emptyRecoCache!(plan::RecoPlan)
  for p in PostOrderDFS(plan)
    if p isa ProcessResultCache || p isa  RecoPlan{ProcessResultCache}
      empty!(p)
    end
  end
end
