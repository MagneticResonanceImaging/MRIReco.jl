export addMRIRecoPlanPath, getMRIRecoPlanList
const DEFAULT_PLANS_PATH = @path joinpath(@__DIR__, "..", "..", "reco_plans")
const recoPlanPaths = AbstractString[DEFAULT_PLANS_PATH]

"""
    addMRIRecoPlanPath(path)

Add all `RecoPlans` within the given directory as potential 
"""
addMRIRecoPlanPath(path) = !(path in recoPlanPaths) ? pushfirst!(recoPlanPaths, path) : nothing
"""
    getMRIRecoPlanList(; full = false)

Retrieve a list of currently known `RecoPlan`s. If `full` is `true` then full paths are returned. 
"""
function getMRIRecoPlanList(; full = false)
  result = String[]
  for path in recoPlanPaths
    if isdir(path)
      push!(result, [plan for plan in filter(a->contains(a,".toml"),readdir(path, join =  full))]...)
    end
  end
  return result
end

export addMRIRecoPlanModule, getMRIRecoPlanModules
const recoPlanModules::Vector{Module} = [AbstractImageReconstruction, MRIOperators, MRIBase, MRIReco, RegularizedLeastSquares]
"""
    addMRIRecoPlanPath(mod::Module)

Add the `mod` module to the list of modules which are used for plan loading 
"""
addMRIRecoPlanModule(mod::Module) = !(mod in recoPlanModules) ? push!(recoPlanModules, mod) : nothing
"""
    getMRIRecoPlanModules()

Retrieve a list of currently used modules for plan loading 
"""
getMRIRecoPlanModules() = recoPlanModules

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
    reconstruct(name::AbstractString, acqData::AcquisitionData; modules = getMRIRecoPlanModules(), kwargs...)

Perform a reconstruction with the `RecoPlan` specified by `name` and given `data`.
Additional keyword arguments can be passed to the reconstruction plan.

`RecoPlans` can be stored in the in the MRIReco package config folder or in a custom folder. New folder can be added with `addMRIRecoPlanPath`. The first plan found is used.
Alternatively, name can be a path to specific plan file.

# Examples
```julia
Ireco = reconstruct("standard", acqData; kwargs...)
```
"""
function AbstractImageReconstruction.reconstruct(name::AbstractString, acqData::AcquisitionData; modules = getMRIRecoPlanModules(), kwargs...)
  plan = loadRecoPlan(name, modules; kwargs...)
  return AbstractImageReconstruction.reconstruct(build(plan), data)
end
# Load plan from a .toml file
function loadRecoPlan(name::AbstractString, modules; kwargs...)
  planfile = planpath(name)
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
# First update plan structure, then set values
function setKwargs!(plan; kwargs...)
  setFirst = filter(kw->kw[2] isa Union{<:AbstractImageReconstructionAlgorithm, <:AbstractImageReconstructionParameters, <:AbstractRecoPlan}, kwargs)
  setAll!(plan; setFirst...)
  setAll!(plan; [kw for kw in kwargs if kw ∉ setFirst]...)
end


export MRIRecoPlan
"""
    MRIRecoPlan(value, modules = getMRIRecoPlanModules(); kwargs...)

Load a `RecoPlan` and set the keyword arguments if applicable. Value can be the name of a plan or path to plan file.
"""
function MRIRecoPlan(value::AbstractString, modules = getMRIRecoPlanModules(); kwargs...)
  if isfile(value) && endswith(value, ".toml")
    return loadRecoPlan(value, modules; kwargs...)
  elseif isfile(value)
    return MRIRecoPlan(MDFFile(value), modules; kwargs...)
  else
    return loadRecoPlan(planpath(value), modules; kwargs...)
  end
end