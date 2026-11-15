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

The dictionary is merged into [`defaultRecoParams`](@ref), and `:reconSize` defaults to the encoding size of
`acqData`. `reconstruction` then loads the plan of the reconstruction type with [`MRIRecoPlan`](@ref), sets the
parameters of the same name as the other keys, for example `iterations` of [`LeastSquaresSolverParameter`](@ref)
for `:iterations`, and runs the plan. Keys that the plan doesn't have are ignored without a warning; with
`ENV["JULIA_DEBUG"] = "MRIReco"`, `reconstruction` lists them. For the iterative reconstructions, a few keys also
choose the type of a parameter: `:densityWeighting => false`, `:weights`, a `:sparseTrafo` given as an operator and
`:encodingOps`; `"direct"` warns about these keys and ignores them, except `:weights` given as a single vector,
which makes it fail.
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
    reconstruction(acqData::AcquisitionData, recoParams::Dict,filename::String; force=false)

Performs the same image reconstruction as `reconstruction(acqData::AcquisitionData, recoParams::Dict)` and is meant
to save the image in a file with name `filename` or, if `force=false` and the file exists, to load the image from
it. Since the functions it calls to load and save the image are no longer part of MRIReco, it throws an
`UndefVarError`, after the reconstruction if the file doesn't exist
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
