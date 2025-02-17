export FLAGS, create_flag_bitmask, flag_is_set, flag_set!, flag_remove!, flag_remove_all!, flags_of

export b64,
       ACQ_FIRST_IN_ENCODE_STEP1,
       ACQ_LAST_IN_ENCODE_STEP1,
       ACQ_FIRST_IN_ENCODE_STEP2,
       ACQ_LAST_IN_ENCODE_STEP2,
       ACQ_FIRST_IN_AVERAGE,
       ACQ_LAST_IN_AVERAGE,
       ACQ_FIRST_IN_SLICE,
       ACQ_LAST_IN_SLICE,
       ACQ_FIRST_IN_CONTRAST,
       ACQ_LAST_IN_CONTRAST,
       ACQ_FIRST_IN_PHASE,
       ACQ_LAST_IN_PHASE,
       ACQ_FIRST_IN_REPETITION,
       ACQ_LAST_IN_REPETITION,
       ACQ_FIRST_IN_SET,
       ACQ_LAST_IN_SET,
       ACQ_FIRST_IN_SEGMENT,
       ACQ_LAST_IN_SEGMENT,
       ACQ_IS_NOISE_MEASUREMENT,
       ACQ_IS_PARALLEL_CALIBRATION,
       ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING,
       ACQ_IS_REVERSE,
       ACQ_IS_NAVIGATION_DATA,
       ACQ_IS_PHASECORR_DATA,
       ACQ_LAST_IN_MEASUREMENT,
       ACQ_IS_HPFEEDBACK_DATA,
       ACQ_IS_DUMMYSCAN_DATA,
       ACQ_IS_RTFEEDBACK_DATA,
       ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA,
       ACQ_IS_PHASE_STABILIZATION_REFERENCE,
       ACQ_IS_PHASE_STABILIZATION,
       ACQ_COMPRESSION1,
       ACQ_COMPRESSION2,
       ACQ_COMPRESSION3,
       ACQ_COMPRESSION4,
       ACQ_USER1,
       ACQ_USER2,
       ACQ_USER3,
       ACQ_USER4,
       ACQ_USER5,
       ACQ_USER6,
       ACQ_USER7,
       ACQ_USER8


const b64 = UInt64(1)
const ACQ_FIRST_IN_ENCODE_STEP1               = 1b64 << ( 01 - 1 )
const ACQ_LAST_IN_ENCODE_STEP1                = 1b64 << ( 02 - 1 )
const ACQ_FIRST_IN_ENCODE_STEP2               = 1b64 << ( 03 - 1 )
const ACQ_LAST_IN_ENCODE_STEP2                = 1b64 << ( 04 - 1 )
const ACQ_FIRST_IN_AVERAGE                    = 1b64 << ( 05 - 1 )
const ACQ_LAST_IN_AVERAGE                     = 1b64 << ( 06 - 1 )
const ACQ_FIRST_IN_SLICE                      = 1b64 << ( 07 - 1 )
const ACQ_LAST_IN_SLICE                       = 1b64 << ( 08 - 1 )
const ACQ_FIRST_IN_CONTRAST                   = 1b64 << ( 09 - 1 )
const ACQ_LAST_IN_CONTRAST                    = 1b64 << ( 10 - 1 )
const ACQ_FIRST_IN_PHASE                      = 1b64 << ( 11 - 1 )
const ACQ_LAST_IN_PHASE                       = 1b64 << ( 12 - 1 )
const ACQ_FIRST_IN_REPETITION                 = 1b64 << ( 13 - 1 )
const ACQ_LAST_IN_REPETITION                  = 1b64 << ( 14 - 1 )
const ACQ_FIRST_IN_SET                        = 1b64 << ( 15 - 1 )
const ACQ_LAST_IN_SET                         = 1b64 << ( 16 - 1 )
const ACQ_FIRST_IN_SEGMENT                    = 1b64 << ( 17 - 1 )
const ACQ_LAST_IN_SEGMENT                     = 1b64 << ( 18 - 1 )
const ACQ_IS_NOISE_MEASUREMENT                = 1b64 << ( 19 - 1 )
const ACQ_IS_PARALLEL_CALIBRATION             = 1b64 << ( 20 - 1 )
const ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING = 1b64 << ( 21 - 1 )
const ACQ_IS_REVERSE                          = 1b64 << ( 22 - 1 )
const ACQ_IS_NAVIGATION_DATA                  = 1b64 << ( 23 - 1 )
const ACQ_IS_PHASECORR_DATA                   = 1b64 << ( 24 - 1 )
const ACQ_LAST_IN_MEASUREMENT                 = 1b64 << ( 25 - 1 )
const ACQ_IS_HPFEEDBACK_DATA                  = 1b64 << ( 26 - 1 )
const ACQ_IS_DUMMYSCAN_DATA                   = 1b64 << ( 27 - 1 )
const ACQ_IS_RTFEEDBACK_DATA                  = 1b64 << ( 28 - 1 )
const ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA   = 1b64 << ( 29 - 1 )
const ACQ_IS_PHASE_STABILIZATION_REFERENCE    = 1b64 << ( 30 - 1 )
const ACQ_IS_PHASE_STABILIZATION              = 1b64 << ( 31 - 1 )
const ACQ_COMPRESSION1                        = 1b64 << ( 53 - 1 )
const ACQ_COMPRESSION2                        = 1b64 << ( 54 - 1 )
const ACQ_COMPRESSION3                        = 1b64 << ( 55 - 1 )
const ACQ_COMPRESSION4                        = 1b64 << ( 56 - 1 )
const ACQ_USER1                               = 1b64 << ( 57 - 1 )
const ACQ_USER2                               = 1b64 << ( 58 - 1 )
const ACQ_USER3                               = 1b64 << ( 59 - 1 )
const ACQ_USER4                               = 1b64 << ( 60 - 1 )
const ACQ_USER5                               = 1b64 << ( 61 - 1 )
const ACQ_USER6                               = 1b64 << ( 62 - 1 )
const ACQ_USER7                               = 1b64 << ( 63 - 1 )
const ACQ_USER8                               = 1b64 << ( 64 - 1 )

"""
    FLAGS::Dict{String, Int}

A dictionary mapping string keys (representing flag names) to bitmask values.
Flags are used to indicate specific attributes of the corresponding Profile. Example flag names include "ACQ_FIRST_IN_ENCODE_STEP1", "ACQ_LAST_IN_SLICE", etc.
"""
FLAGS = Dict(
    "ACQ_FIRST_IN_ENCODE_STEP1"                => 1,
    "ACQ_LAST_IN_ENCODE_STEP1"                 => 2,
    "ACQ_FIRST_IN_ENCODE_STEP2"                => 3,
    "ACQ_LAST_IN_ENCODE_STEP2"                 => 4,
    "ACQ_FIRST_IN_AVERAGE"                     => 5,
    "ACQ_LAST_IN_AVERAGE"                      => 6,
    "ACQ_FIRST_IN_SLICE"                       => 7,
    "ACQ_LAST_IN_SLICE"                        => 8,
    "ACQ_FIRST_IN_CONTRAST"                    => 9,
    "ACQ_LAST_IN_CONTRAST"                     => 10,
    "ACQ_FIRST_IN_PHASE"                       => 11,
    "ACQ_LAST_IN_PHASE"                        => 12,
    "ACQ_FIRST_IN_REPETITION"                  => 13,
    "ACQ_LAST_IN_REPETITION"                   => 14,
    "ACQ_FIRST_IN_SET"                         => 15,
    "ACQ_LAST_IN_SET"                          => 16,
    "ACQ_FIRST_IN_SEGMENT"                     => 17,
    "ACQ_LAST_IN_SEGMENT"                      => 18,
    "ACQ_IS_NOISE_MEASUREMENT"                 => 19,
    "ACQ_IS_PARALLEL_CALIBRATION"              => 20,
    "ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING"  => 21,
    "ACQ_IS_REVERSE"                           => 22,
    "ACQ_IS_NAVIGATION_DATA"                   => 23,
    "ACQ_IS_PHASECORR_DATA"                    => 24,
    "ACQ_LAST_IN_MEASUREMENT"                  => 25,
    "ACQ_IS_HPFEEDBACK_DATA"                   => 26,
    "ACQ_IS_DUMMYSCAN_DATA"                    => 27,
    "ACQ_IS_RTFEEDBACK_DATA"                   => 28,
    "ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA"    => 29,
    "ACQ_COMPRESSION1"                         => 53,
    "ACQ_COMPRESSION2"                         => 54,
    "ACQ_COMPRESSION3"                         => 55,
    "ACQ_COMPRESSION4"                         => 56,
    "ACQ_USER1"                                => 57,
    "ACQ_USER2"                                => 58,
    "ACQ_USER3"                                => 59,
    "ACQ_USER4"                                => 60,
    "ACQ_USER5"                                => 61,
    "ACQ_USER6"                                => 62,
    "ACQ_USER7"                                => 63,
    "ACQ_USER8"                                => 64
)

function bitshift(A, k)
  k >= 0 ? out = A << k : out = A >> abs(k)
  return out
end

"""
    create_flag_bitmask(flag::AbstractString) -> UInt64
    create_flag_bitmask(flag::Integer) -> UInt64

Creates a bitmask for the given flag.

- If the flag is a string, its integer value is looked up in `FLAGS`.
- If the flag is an integer, it must be positive. The bitmask is created by left-shifting `1` to the corresponding bit position.

# Arguments
- `flag`: A string representing the flag name or a positive integer.

# Returns
- A `UInt64` bitmask where only the bit corresponding to the flag is set.

# Throws
- `KeyError` if the string flag is not found in `FLAGS`.
- `DomainError` if the integer flag is not positive.

# Example
```julia
create_flag_bitmask("ACQ_FIRST_IN_ENCODE_STEP1")  # 0x0000000000000001
create_flag_bitmask(2)                           # 0x0000000000000002
```
"""
create_flag_bitmask(flag::T) where T  = error("Unexpected type for bitmask, expected String or positive Integer, found $T")

create_flag_bitmask(flag::AbstractString) = create_flag_bitmask(FLAGS[flag])
function create_flag_bitmask(flag::Integer)
  (flag > 0 && flag <= 64)|| throw(DomainError(flag, "Bitmask can only be created for  integers from 1 to 64"))
  b = UInt64(flag)
  return bitshift(UInt64(1),  b - 1)
end

"""
    flag_is_set(obj::Profile, flag) -> Bool
    flag_is_set(head::AcquisitionHeader, flag) -> Bool

Checks whether a specific flag is set in the given `Profile` or `AcquisitionHeader` object.

# Arguments
- `obj`: A `Profile` or `AcquisitionHeader` object field containing a `flags` attribute (as a bitmask).
- `flag`: A string (flag name) or an integer (bit position).

# Returns
- `true` if the flag is set in the `obj.head.flags` bitmask; otherwise, `false`.

# Example
```julia
flag_is_set(profile, "ACQ_FIRST_IN_ENCODE_STEP1")  # true or false
flag_is_set(head, "ACQ_FIRST_IN_ENCODE_STEP1")  # true or false
```
"""
function flag_is_set(obj::Profile, flag)
  return flag_is_set(obj.head, flag)
end

function flag_is_set(head::AcquisitionHeader, flag)
  bitmask = create_flag_bitmask(flag)
  return head.flags & bitmask > 0
end

"""
    flag_set!(obj::Profile, flag)
    flag_set!(head::AcquisitionHeader, flag)

Sets the specified flag in the given `Profile` or `AcquisitionHeader` object.

# Arguments
- `obj`: A `Profile` object with a `head` field containing a `flags` attribute (as a bitmask).
- `head`: An `AcquisitionHeader` object containing a `flags` attribute (as a bitmask).
- `flag`: A string (flag name) or an integer (bit position), or a vector of such values.

# Modifies
- `obj.head.flags` or `head.flags` by setting the bit corresponding to the flag.

# Example
```julia
flag_set!(profile, "ACQ_FIRST_IN_SLICE")
flag_set!(profile, ["ACQ_LAST_IN_PHASE","ACQ_USER8"])

head = AcquisitionHeader()
flag_set!(head, "ACQ_FIRST_IN_SLICE")
```
"""
function flag_set!(obj::Profile, flags)
  flag_set!(obj.head, flags)
end

function flag_set!(head::AcquisitionHeader, flags::Union{T,Vector{T}}) where T <: Union{Int, AbstractString}
    broadcast(flags) do flag
      bitmask = create_flag_bitmask(flag)
      head.flags = head.flags | bitmask
    end
end

"""
    flag_remove!(obj::Profile, flag)
    flag_remove!(head::AcquisitionHeader, flag)

Removes (clears) the specified flag in the given `Profile` or `AcquisitionHeader` object.

# Arguments
- `obj`: A `Profile` object with a `head` field containing a `flags` attribute (as a bitmask).
- `head`: An `AcquisitionHeader` object containing a `flags` attribute (as a bitmask).
- `flag`: A string (flag name) or an integer (bit position), or a vector of such values.

# Modifies
- `obj.head.flags` or `head.flags` by clearing the bit corresponding to the flag.

# Example
```julia
flag_remove!(profile, "ACQ_LAST_IN_PHASE")
flag_remove!(profile, ["ACQ_LAST_IN_PHASE","ACQ_USER8"])

head = AcquisitionHeader()
flag_remove!(head, "ACQ_LAST_IN_PHASE")
```
"""
function flag_remove!(obj::Profile, flag)
  flag_remove!(obj.head, flag)
end

function flag_remove!(head::AcquisitionHeader, flags::Union{T,Vector{T}}) where T <: Union{Int, AbstractString}
  broadcast(flags) do flag
    bitmask = create_flag_bitmask(flag)
    head.flags = head.flags & ~bitmask
  end
end

"""
    flag_remove_all!(obj::Profile)
    flag_remove_all!(head::AcquisitionHeader)

Clears all flags in the given `Profile` object or `AcquisitionHeader` object.

# Arguments
- `obj`: A `Profile` object with a `head` field containing a `flags` attribute (as a bitmask).
- `head`: An `AcquisitionHeader` object containing a `flags` attribute (as a bitmask).

# Modifies
- `obj.head.flags` or `head.flags`, setting it to `0` (all flags cleared).

# Example
```julia
flag_remove_all!(profile)
flag_remove_all!(head)
```
"""
function flag_remove_all!(obj::Profile)
    flag_remove_all!(obj.head)
end

function flag_remove_all!(head::AcquisitionHeader)
  head.flags = UInt64(0);
end

"""
    flags_of(obj::Profile) -> Vector{String}
    flags_of(head::AcquisitionHeader) -> Vector{String}

Returns a list of all flags that are set in the given `Profile` object or `AcquisitionHeader` object.

# Arguments
- `obj`: A `Profile` object with a `head` field containing a `flags` attribute (as a bitmask).
- `head`: An `AcquisitionHeader` object containing a `flags` attribute (as a bitmask).

# Returns
- An array of strings representing the names of the flags that are currently set.

# Example
```julia
flags = flags_of(profile)  # ["ACQ_FIRST_IN_ENCODE_STEP1", "ACQ_LAST_IN_SLICE"]
```
"""
function flags_of(obj::Profile)
  return flags_of(obj.head)
end

function flags_of(head::AcquisitionHeader)
  flags_ = String[]
  for f in FLAGS
    flag_is_set(head::AcquisitionHeader, f.first) ? push!(flags_,f.first) : nothing
  end
  return flags_
end

## SANITY Check for head.flags -> should be updated by a UInt64
function Base.setproperty!(header::AcquisitionHeader, name::Symbol, val)
  if name == :flags
    return throw(DomainError("Flags can only be updated with UInt64 found $(typeof(val)) \n \n Use the dedicated interface `h.flags = ACQ_IS_REVERSE + ACQ_IS_NAVIGATION_DATA` or the function flag_set!(h,\"ACQ_IS_REVERSE\")"))
  end
  ty = fieldtype(typeof(header), name)
  tmp = val isa ty ? val : convert(ty, val)
  return setfield!(header, name, tmp)
end
function Base.setproperty!(header::AcquisitionHeader, name::Symbol, val::UInt64)
  ty = fieldtype(typeof(header), name)
  tmp = val isa ty ? val : convert(ty, val)
  return setfield!(header, name, tmp)
end