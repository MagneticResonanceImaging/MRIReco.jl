export reconstruction

include("RecoParameters.jl")
include("DirectReconstruction.jl")
include("IterativeReconstruction.jl")


"""
    reconstruction(acqData::AcquisitionData, recoParams::Dict)

Performs image reconstruction of an AcquisitionData object.
Parameters are specified in a dictionary.

Reconstruction types are specified by the key `:reco`.
Valid reconstruction names are:
* `"direct"`             - direct Fourier reconstruction (default)
* `"standard"`           - iterative reconstruction for all contrasts, coils & slices independently
* `"multiEcho"`          - iterative joint reconstruction of all echo images
* `"multiCoil"`          - SENSE-type iterative reconstruction
* `"multiCoilMultiEcho"` - SENSE-type iterative reconstruction of all echo images
* `"multiCoilMultiEchoSubspace"` - SENSE-type iterative reconstruction of the coefficient maps of a temporal basis

The dictionary is merged into the defaults of [`defaultRecoParams`](@ref), and `:reconSize` defaults to the encoding
size of `acqData`. `reconstruction` then loads the plan of the reconstruction type with [`MRIRecoPlan`](@ref). For the
iterative reconstructions, a few keys first select the type of a parameter: `:densityWeighting => false` selects
[`UniformWeightingParameters`](@ref), `:weights` selects [`CustomWeightingParameters`](@ref), `:encodingOps` selects
[`CustomEncodingParameters`](@ref), and a `:sparseTrafo` with an operator selects
[`CustomSparsityParameters`](@ref); `"direct"` does not support these keys and warns that it cannot set the
parameter they select. `reconstruction` then sets the parameters of the same name as the other keys, for example
`iterations` of [`LeastSquaresSolverParameter`](@ref) for `:iterations`, and runs the plan.

Other keys that the plan does not have are ignored. If the value of such a key is a parameter object such as
`UniformWeightingParameters()`, an algorithm or a plan, `reconstruction` warns; with
`ENV["JULIA_DEBUG"] = "MRIReco"`, it lists the remaining keys that the plan does not have, `:densityWeighting`
included.
"""
function reconstruction(acqData::AcquisitionData, recoParams::Dict)
  recoParams = copy(recoParams)
  # check dimensionality of encoding
  encodingDims = ndims(trajectory(acqData))
  if encodingDims==3 && numSlices(acqData)>1
    @error "reconstruction of multiple 3d-encoded volumina is not yet supported"
  end

  if !haskey(recoParams, :reconSize)
    recoParams[:reconSize] = encodingSize(acqData)
  end  

  # load reconstruction parameters
  recoParams = merge(defaultRecoParams(), recoParams)

  plan = setupPlanReco!(acqData, recoParams)
  return reconstruct(build(plan), acqData)
end

"""
    reconstruction(acqData::AcquisitionData, recoParams::Dict, filename::String;
                   force=false)

This method is meant to save the result of `reconstruction(acqData, recoParams)` in the file `filename`, or to load
it from there if `force=false` and the file exists. The functions it calls for this are no longer part of MRIReco,
so it throws an `UndefVarError`: right away if `force=false` and the file exists, otherwise after the reconstruction
([#295](https://github.com/MagneticResonanceImaging/MRIReco.jl/issues/295)).
"""
function reconstruction(acqData::AcquisitionData, recoParams::Dict, filename::String;
                        force=false)
  if !force && isfile(filename)
    return recoImage( RecoFileIBI(filename) )
  else
    I = reconstruction(acqData, recoParams)
    saveasRecoFile(filename, I, recoParams)
    return I
  end
end
