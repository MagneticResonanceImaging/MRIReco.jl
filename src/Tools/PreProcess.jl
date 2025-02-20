export remove_oversampling, remove_oversampling!

using MRIOperators.FFTW
"""
    raw = remove_oversampling(rawData::RawAcquisitionData)
    remove_oversampling!(rawData::RawAcquisitionData)

Removes oversampling from the profiles in a `RawAcquisitionData` object by creating a deep copy of the input data and applying the `remove_oversampling!` function.

# Arguments
- `rawData::RawAcquisitionData`: The `RawAcquisitionData` object to process.

# Returns
- `RawAcquisitionData`: A new `RawAcquisitionData` object with oversampling removed.

# Examples
```julia
# Assuming `rawData` is a RawAcquisitionData object
processedData = remove_oversampling(rawData)
```
"""
function remove_oversampling(rawData)
  rawData_2 = deepcopy(rawData)
  remove_oversampling!(rawData_2)
  return rawData_2
end

function remove_oversampling!(rawData)
  ratioFOV = rawData.params["encodedFOV"][1] / rawData.params["reconFOV"][1]
    
  # Assuming all profiles have the same size
  sdim = size(rawData.profiles[1].data,1)
  sdim_out = Int(sdim / ratioFOV)
  start = Int((sdim[1] - sdim_out) ÷ 2)

  Threads.@threads for profile in rawData.profiles
      tmp = profile.data
      tmp = ifftshift(ifft(ifftshift(tmp, 1), 1), 1)
      @views tmp = tmp[start:start + sdim_out - 1, :]
      tmp = fftshift(fft(fftshift(tmp, 1), 1), 1)
      
      profile.data = tmp
      profile.head.number_of_samples = sdim_out
      profile.head.center_sample = UInt16(profile.head.center_sample / ratioFOV)
      profile.head.discard_pre = UInt16(profile.head.discard_pre / ratioFOV)
      profile.head.discard_post = UInt16(profile.head.discard_post / ratioFOV)
  end

  rawData.params["encodedSize"][1] = rawData.params["encodedSize"][1] ÷ ratioFOV
  return rawData
end
