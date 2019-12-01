module HiGHS

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
end

include("c_model.jl")
end # module CWrapper

end # module HiGHS
