module HiGHS

export ManagedHiGHS

import HiGHS_jll: libhighs

"""
Wrapper for the HiGHS C library.
"""
module CWrapper

import HiGHS.libhighs
include(joinpath("wrapper", "libhighs_api.jl"))
include(joinpath("wrapper", "libhighs_common.jl"))
include(joinpath("wrapper", "api_helpers.jl"))

end # module CWrapper

include("c_model.jl")
include("options.jl")

include("MOI_wrapper.jl")

end # module HiGHS
