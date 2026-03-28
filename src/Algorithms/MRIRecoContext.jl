"""
    MRIRecoContext

Container for reconstruction context parameters that are used throughout the 
reconstruction pipeline. Uses a ScopedValue for task-local access.

# Fields
* `reconSize::R` - Size of the reconstructed image (NTuple)
* `arrayType::A` - Array type for computation (Array, JLArray, etc.)
* `S::S` - Storage type for operators
* `acqData::Acq` - Acquisition data (optional), provides metadata like trajectory info

# Usage
```julia
# Set the context for reconstruction
ctx = MRIRecoContext((128,128), Float64, Array)

# Or with acquisition data
ctx = MRIRecoContext((128,128), acqData, Array)

# Using with() to set scoped value
with(MRIRECO_CONTEXT => ctx) do
    # reconstruction code
end

# And then access anywhere in the reconstruction pipeline
function foo(sparsity_params)
    rs = ctx_reconSize()  # Get from scoped value
    ad = ctx_acqData()    # Get acquisition data if needed
    # ...
end
```
"""
struct MRIRecoContext{D, arrT, S, Acq <: Union{Nothing, AcquisitionData}}
  reconSize::NTuple{D, Int64}
  arrayType::Type{arrT}
  storageType::Type{S}
  acqData::Acq
  
    function MRIRecoContext(reconSize::NTuple{D, Int64}, arrayType::Type{arrT}, storageType::Type{S}, acqData::Acq) where {D, arrT, S, Acq}
    return new{D, arrT, S, Acq}(reconSize, arrayType, storageType, acqData)
  end
end

function MRIRecoContext(reconSize::NTuple{D, Int64}, acqData::AcquisitionData{T}, arrayType::Type{arrT} = Array) where {D, T, arrT}
  storageType = typeof(arrayType{Complex{T}}(undef, 0))  
  return MRIRecoContext(reconSize, arrayType, storageType, acqData)
end
function MRIRecoContext(reconSize::NTuple{D, Int64}, eltype::Type{T}, arrayType::Type{arrT} = Array) where {D, T <: Real, arrT}
  storageType = typeof(arrayType{Complex{T}}(undef, 0))  
  return MRIRecoContext(reconSize, arrayType, storageType, nothing)
end

const MRIRECO_CONTEXT = ScopedValue{MRIRecoContext}()

export MRIRECO_CONTEXT, MRIRecoContext
export ctx_reconSize, ctx_arrayType, ctx_storageType, ctx_acqData

"""
    ctx_reconSize() -> NTuple

Get the reconstruction size from the current MRI reconstruction context.
"""
ctx_reconSize() = MRIRECO_CONTEXT[].reconSize

"""
    ctx_arrayType() -> Type

Get the array type from the current MRI reconstruction context.
"""
ctx_arrayType() = MRIRECO_CONTEXT[].arrayType

"""
    ctx_storageType() -> Type

Get the storage type from the current MRI reconstruction context.
"""
ctx_storageType() = MRIRECO_CONTEXT[].storageType

"""
    ctx_acqData() -> Union{AcquisitionData, Nothing}

Get the acquisition data from the current MRI reconstruction context.
Returns nothing if not set.
"""
ctx_acqData() = MRIRECO_CONTEXT[].acqData