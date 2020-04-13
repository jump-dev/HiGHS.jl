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

module MOIWrapper

import HiGHS.CWrapper
import HiGHS: ManagedHiGHS, reset_model!

include("MOI_wrapper.jl")
end # module MOIWrapper

import .MOIWrapper: Optimizer

end # module HiGHS
