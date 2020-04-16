# Completes the generated wrapper with default methods

# TODO update these to leverage dispatch

"""
    set_option(highs, option, value::T)

Sets the option to the value, dispatches on
the value type T to call appropriate HiGHS function.
"""
function set_option end

function set_option(highs, option, value::Integer)
    Highs_setHighsIntOptionValue(highs, option, Cint(value))
end

function set_option(highs, option, value::Bool)
    Highs_setHighsBoolOptionValue(highs, option, Cint(value))
end

function set_option(highs, option, value::AbstractFloat)
    Highs_setHighsDoubleOptionValue(highs, option, Cdouble(value))
end

function set_option(highs, option, value::String)
    Highs_setHighsStringOptionValue(highs, option, value)
end

"""
    get_option(highs, option, T)

Gets the option of type T on the HiGHS side.
"""
function get_option end

function get_option(highs, option, ::Type{T}) where {T <: Integer}
    value = Ref{Cint}(0)
    Highs_getHighsIntOptionValue(highs, option, value)
    return T(value[])
end

function get_option(highs, option, ::Type{Bool})
    value = Ref{Cint}(0)
    Highs_getHighsBoolOptionValue(highs, option, value)
    return Bool(value[])
end

function get_option(highs, option, ::Type{T}) where {T <: AbstractFloat}
    value = Ref{Cdouble}()
    Highs_getHighsDoubleOptionValue(highs, option, value)
    return T(value[])
end

function get_option(highs, option, ::Type{String})
    v = Vector{Cchar}(undef, 100)
    p = pointer(v)
    Highs_getHighsStringOptionValue(highs, option, p)
    GC.@preserve v s = unsafe_string(p)
    return s
end

# convenience to accomodate/ignore scale argument
function Highs_getModelStatus(highs)
    Highs_getModelStatus(highs, Cint(0))
end
