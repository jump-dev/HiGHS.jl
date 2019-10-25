module HiGHS

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("HiGHS not properly installed. Please run import Pkg; Pkg.build(\"HiGHS\")")
end

include(joinpath("wrapper", "libhighs_api.jl"))
include(joinpath("wrapper", "libhighs_common.jl"))

end # module
