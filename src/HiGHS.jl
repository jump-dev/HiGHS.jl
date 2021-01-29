module HiGHS

import HiGHS_jll: libhighs

include("wrapper/ctypes.jl")
include("wrapper/libhighs_common.jl")
include("wrapper/libhighs_api.jl")

include("MOI_wrapper.jl")

# HiGHS exports all `Highs_xxx` symbols. If you don't want all of these symbols
# in your environment, then use `import HiGHS` instead of `using HiGHS`.

for sym in names(@__MODULE__, all = true)
    if startswith(string(sym), "Highs_")
        @eval export $sym
    end
end

end # module HiGHS
