# # [Configure a Reconstruction with Plans](@id configure-plans)
#
# [`reconstruction`](@ref) builds every reconstruction from a reconstruction plan of
# [AbstractImageReconstruction.jl](https://juliaimagerecon.github.io/AbstractImageReconstruction.jl/dev/), as
# described in [Reconstruction Parameters](@ref reconstruction-parameters). The plan itself shows all parameters of a
# reconstruction type, and with it, a configured reconstruction can be run on several datasets, its structure can be
# saved as a new reconstruction type, and a part of it can be replaced by another built-in part or by one of your own.
# This guide shows how. MRIReco re-exports the functions for plans, such as `setAll!` and `build`; the documentation
# of AbstractImageReconstruction.jl explains plans in general.
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
#
# ## Implementing a Part of Your Own
#
# The parts of a plan are what AbstractImageReconstruction.jl calls parameters: callable objects that each carry out
# one step of the reconstruction, as its documentation describes under
# [Interface](https://juliaimagerecon.github.io/AbstractImageReconstruction.jl/dev/generated/example/1_interface/).
#
# A weighting of the iterative reconstruction types is a subtype of [`AbstractMRIRecoWeightingParameters`](@ref),
# defined with `@parameter`, with a method `(weighting::MyWeighting)(::Type{<:AbstractMRIRecoAlgorithm})`. The
# reconstruction calls the method and expects one vector per contrast with one weight per sample, the same for all
# coils, in the form that [`samplingDensity`](@ref) returns. The vectors must have the type that
# [`ctx_storageType`](@ref) returns, complex and in the precision of the data. The reconstruction weights the encoding
# operator and the data with them, so the weights take the place of the square roots of the density compensation
# factors, not of the factors themselves, see [Density Compensation](@ref). Their scale matters as well: with the
# default `:normalizeReg`, weights `c` times larger scale the data term by `c^2` relative to the regularization.
#
# The following weighting has one field: with `exponent = 1`, the default, it returns the weights of
# `DensityWeightingParameters`, with `exponent = 0` the constant weight of `UniformWeightingParameters`, and with a
# value in between it applies the density compensation only in part, with weights between those of the two.
# [`ctx_acqData`](@ref) and [`ctx_reconSize`](@ref) give its method the data and the image size of the
# running reconstruction:

@parameter struct PartialDensityWeightingParameters <: AbstractMRIRecoWeightingParameters
    exponent::Float64 = 1.0
end

function (weighting::PartialDensityWeightingParameters)(::Type{<:AbstractMRIRecoAlgorithm})
    exponent = weighting.exponent
    uniformWeight = 1 / sqrt(prod(ctx_reconSize()))
    return map(samplingDensity(ctx_acqData(), ctx_reconSize())) do densityWeights
        weights = uniformWeight^(1 - exponent) .* densityWeights .^ exponent
        ctx_storageType()(weights)
    end
end

# With the default, the reconstruction computes the same image as with the density weighting above:

ownPlan = MRIRecoPlan("standard"; reconSize = (N, N), solver = CGNR,
                      reg = L2Regularization(1e-3), iterations = 10,
                      weightingParams = RecoPlan(PartialDensityWeightingParameters))
reconstruct(build(ownPlan), acqData) ≈ image

# The field of the type is a parameter of the plan, which `setAll!` sets like the built-in ones, so it needs a name
# that no other parameter of the plan has:

setAll!(ownPlan; exponent = 0.5)
ownPlan.parameter.parameter.weightingParams

# `reconstruction` takes an object of the type as `:weightingParams`. `:densityWeighting = false` from above would
# replace it and is removed first. Both ways compute the same image, which now differs from the one with the density
# weighting:

delete!(params, :densityWeighting)
params[:weightingParams] = PartialDensityWeightingParameters(exponent = 0.5)
partial = reconstruct(build(ownPlan), acqData)
(partial ≈ reconstruction(acqData, params), partial ≈ image)

# For weights that are fixed numbers, the key `:weights` is enough, currently only for data with one coil; a type
# computes them for the data that the plan is applied to. A plan with the type is saved like the one above and, for a
# type that a script defines, loads after `addMRIRecoPlanModule(Main)`.
#
# Other parts take objects of your own without a new type: an operator as `:sparseTrafo`, see
# [Use a Custom Sparsifying Transform](@ref custom-transform), and a denoiser in `PnPRegularization` as `:reg`, see
# [Use a Plug-and-Play Prior](@ref plug-and-play).
