export filter_raw_by_flags

"""
    filter_raw_by_flags(rawData::RawAcquisitionData, flags::Vector{String})

Filters the profiles in a `RawAcquisitionData` object based on specified flags.

# Arguments
- `rawData::RawAcquisitionData`: The `RawAcquisitionData` object to filter.
- `flags::Vector{String}`: A vector of flag names to filter the profiles by.

# Returns
- `RawAcquisitionData`: A new `RawAcquisitionData` object containing only the profiles that have all the specified flags set.

# Examples
```julia
# Assuming `rawData` is a RawAcquisitionData object and `flags` is a vector of flag names
filteredData = filter_raw_by_flags(rawData, ["ACQ_IS_PHASECORR", "ACQ_IS_NOISEADJSCAN"])
"""
function filter_raw_by_flags(rawData::RawAcquisitionData, flags::Union{T,Vector{T}}) where T <: Union{Int, AbstractString}
  filteredProfiles = Profile[]
  for profile in rawData.profiles
    if all(flag_is_set(profile, flag) for flag in flags)
      push!(filteredProfiles, profile)
    end
  end
  return RawAcquisitionData(rawData.params, filteredProfiles)
end

