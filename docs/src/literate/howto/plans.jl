# # [Configure a Reconstruction with Plans](@id configure-plans)
#
# [`reconstruction`](@ref) builds every reconstruction from a reconstruction plan of
# [AbstractImageReconstruction.jl](https://juliaimagerecon.github.io/AbstractImageReconstruction.jl/dev/), as
# described in [Reconstruction Parameters](@ref reconstruction-parameters). The plan itself shows all parameters of a
# reconstruction type, and with it, a configured reconstruction can be run on several datasets, and its structure can
# be saved as a new reconstruction type. This guide shows how. MRIReco re-exports the functions for plans, such as
# `setAll!` and `build`; the documentation of AbstractImageReconstruction.jl explains plans in general.
#
# ## Configuring a Plan
#
# [`MRIRecoPlan`](@ref) loads the plan of a reconstruction type by its name. The loaded plan has the structure of the
# reconstruction, but none of its parameters has a value yet:

using MRIReco, MRISimulation
using ImageUtils: shepp_logan

plan = MRIRecoPlan("standard")

# `setAll!` sets a parameter wherever it occurs in the plan, by the name that the plan shows, and only warns about a
# name that the plan doesn't have; [Reconstruction Parameters](@ref reconstruction-parameters) describes the
# parameters. A parameter can also be set at its place in the plan:

setAll!(plan; solver = CGNR, reg = L2Regularization(1e-3))
plan.parameter.parameter.solverParams.iterations = 10
plan.parameter.parameter.solverParams

# `MRIRecoPlan` takes values as keyword arguments as well, for example `MRIRecoPlan("standard"; solver = CGNR)`.
# Parameters without a value take their defaults when the plan is built. A few, such as the coil sensitivities of
# `"multiCoil"`, have no default and have to be set.
#
# ## Running the Reconstruction
#
# `build` turns the plan into a reconstruction algorithm, and `reconstruct` applies it to data, here to simulated
# radial data of a phantom:

N = 128
phantom = shepp_logan(N)
simParams = Dict{Symbol, Any}(:simulation => "fast", :trajName => "Radial",
                              :numProfiles => 64, :numSamplingPerProfile => 2N)
acqData = simulation(phantom, simParams)

setAll!(plan; reconSize = (N, N))
algo = build(plan)
image = reconstruct(algo, acqData)
size(image)

# `reconstruction` with the same keys computes the same image:

params = Dict{Symbol, Any}(:reco => "standard", :reconSize => (N, N), :solver => CGNR,
                           :reg => L2Regularization(1e-3), :iterations => 10)
image == reconstruction(acqData, params)

# The algorithm reconstructs further data with the same settings, for example `reconstruct(algo, otherData)`.
# [`reconstruct`](@ref) with the name of a reconstruction type, `reconstruct("standard", acqData; kwargs...)`, loads,
# configures, builds and runs a plan in one call.
#
# ## Replacing a Part of the Reconstruction
#
# Some parameters of a plan are nested plans, for example the weighting of the data. Replacing such a plan changes how
# the reconstruction works. With a plan of `UniformWeightingParameters`, all samples have the same weight, as with
# `:densityWeighting = false` in `reconstruction`; like that key, it currently works only for data with one coil, see
# the limits in [Reconstruction Parameters](@ref reconstruction-parameters). The part is set as a plan, not as
# `UniformWeightingParameters()`, so that it stays when the plan is cleared below:

setAll!(plan; weightingParams = RecoPlan(UniformWeightingParameters))
uniform = reconstruct(build(plan), acqData)
params[:densityWeighting] = false
uniform == reconstruction(acqData, params)

# ## Saving a Plan as a Reconstruction Type
#
# A plan without parameter values, like the plans in the folder `reco_plans` of MRIReco.jl, serves as the template of
# a reconstruction type. `clear!` removes the values and keeps the structure, that is, the nested plans, and
# `savePlan` writes the plan to a TOML file. After [`addMRIRecoPlanPath`](@ref) with the folder of the file, the new
# type can be used by the name of the file without `.toml`:

clear!(plan)
dir = mktempdir()
savePlan(joinpath(dir, "uniformStandard.toml"), plan)
addMRIRecoPlanPath(dir)
reconstruct("uniformStandard", acqData; reconSize = (N, N), solver = CGNR,
            reg = L2Regularization(1e-3), iterations = 10) == uniform

# `reconstruction` takes the name as `:reco` as well, and the path of the file works without `addMRIRecoPlanPath`.
# Added folders are searched first, the last one added first, and stay added until Julia is restarted, so a plan file
# there can also replace a built-in type. [`getMRIRecoPlanList`](@ref) lists the plan files of all folders. A plan
# with parameter types from another module loads after [`addMRIRecoPlanModule`](@ref) with that module, which has to
# export the types. Save plans without values, as above, and pass the values when the plan is used: not all values
# can be saved and loaded yet, for example the solver, the regularization terms and `reconSize`
# ([#397](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/397)), and `reconstruction` replaces
# stored values such as `:iterations` with its defaults.
