module HiGHS

export ManagedHiGHS

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("HiGHS not properly installed. Please run import Pkg; Pkg.build(\"HiGHS\")")
end

"""
Wrapper for the HiGHS C library.
"""
module CWrapper

import HiGHS.libhighs
include(joinpath("wrapper", "libhighs_api.jl"))
include(joinpath("wrapper", "libhighs_common.jl"))
end # module CWrapper

include("c_model.jl")

module MOIWrapper

import HiGHS.CWrapper
import HiGHS: ManagedHiGHS

include("MOI_wrapper.jl")
end # module MOIWrapper

import .MOIWrapper: Optimizer

end # module HiGHS
