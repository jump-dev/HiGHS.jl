# Completes the generated wrapper with default methods
# Automatically generated using Clang.jl

# TODO update these to leverage dispatch
# function Highs_setHighsBoolOptionValue(highs, option, value::Cint)
    # ccall((:Highs_setHighsBoolOptionValue, libhighs), Cint, (Ptr{Cvoid}, Cstring, Cint), highs, option, value)
# end
# 
# function Highs_setHighsIntOptionValue(highs, option, value::Cint)
    # ccall((:Highs_setHighsIntOptionValue, libhighs), Cint, (Ptr{Cvoid}, Cstring, Cint), highs, option, value)
# end
# 
# function Highs_setHighsDoubleOptionValue(highs, option, value::Cdouble)
    # ccall((:Highs_setHighsDoubleOptionValue, libhighs), Cint, (Ptr{Cvoid}, Cstring, Cdouble), highs, option, value)
# end
# 
# function Highs_setHighsStringOptionValue(highs, option, value)
    # ccall((:Highs_setHighsStringOptionValue, libhighs), Cint, (Ptr{Cvoid}, Cstring, Cstring), highs, option, value)
# end
# 
# function Highs_setHighsOptionValue(highs, option, value)
    # ccall((:Highs_setHighsOptionValue, libhighs), Cint, (Ptr{Cvoid}, Cstring, Cstring), highs, option, value)
# end
# 
# function Highs_getHighsBoolOptionValue(highs, option, value)
    # ccall((:Highs_getHighsBoolOptionValue, libhighs), Cint, (Ptr{Cvoid}, Cstring, Ptr{Cint}), highs, option, value)
# end
# 
# function Highs_getHighsIntOptionValue(highs, option, value)
    # ccall((:Highs_getHighsIntOptionValue, libhighs), Cint, (Ptr{Cvoid}, Cstring, Ptr{Cint}), highs, option, value)
# end
# 
# function Highs_getHighsDoubleOptionValue(highs, option, value)
    # ccall((:Highs_getHighsDoubleOptionValue, libhighs), Cint, (Ptr{Cvoid}, Cstring, Ptr{Cdouble}), highs, option, value)
# end
# 
# function Highs_getHighsStringOptionValue(highs, option, value)
    # ccall((:Highs_getHighsStringOptionValue, libhighs), Cint, (Ptr{Cvoid}, Cstring, Cstring), highs, option, value)
# end
# 
# function Highs_getHighsIntInfoValue(highs, info, value)
    # ccall((:Highs_getHighsIntInfoValue, libhighs), Cint, (Ptr{Cvoid}, Cstring, Ptr{Cint}), highs, info, value)
# end
# 
# function Highs_getHighsDoubleInfoValue(highs, info, value)
    # ccall((:Highs_getHighsDoubleInfoValue, libhighs), Cint, (Ptr{Cvoid}, Cstring, Ptr{Cdouble}), highs, info, value)
# end
# 
# function Highs_getSolution(highs, colvalue, coldual, rowvalue, rowdual)
    # ccall((:Highs_getSolution, libhighs), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), highs, colvalue, coldual, rowvalue, rowdual)
# end
# 
# function Highs_getBasis(highs, colstatus, rowstatus)
    # ccall((:Highs_getBasis, libhighs), Cvoid, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}), highs, colstatus, rowstatus)
# end

function Highs_getModelStatus(highs)
    Highs_getModelStatus(highs, Cint(0))
end
