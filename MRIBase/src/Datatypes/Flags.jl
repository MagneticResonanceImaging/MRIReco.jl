export FLAGS, create_flag_bitmask, flag_is_set, flag_set!, flag_remove!, flag_remove_all!
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

create_flag_bitmask(flag::T) where T  = error("Unexpected type for bitmask, expected String or positive Integer, found $T")
create_flag_bitmask(flag::AbstractString) = create_flag_bitmask(FLAGS[flag])
function create_flag_bitmask(flag::Integer)
  flag > 0 || throw(DomainError(flag, "Bitmask can only be created for positive integers"))
  b = UInt64(flag)
  return bitshift(UInt64(1),  b - 1)
end

function flag_is_set(obj::Profile, flag)
  bitmask = create_flag_bitmask(flag)
  ret = obj.head.flags & bitmask > 0
  return ret
end


function flag_set!(obj::Profile, flag)
  bitmask = create_flag_bitmask(flag)
  obj.head.flags = obj.head.flags | bitmask
end

function flag_remove!(obj::Profile, flag)
  bitmask = create_flag_bitmask(flag)
  obj.head.flags = obj.head.flags & ~bitmask
end

function flag_remove_all!(obj::Profile)
    obj.head.flags = UInt64(0);
end