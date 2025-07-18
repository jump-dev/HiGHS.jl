# Copyright (c) 2019 Mathieu Besançon, Oscar Dowson, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

const CleverDicts = MOI.Utilities.CleverDicts

@enum(
    _RowType,
    _ROW_TYPE_LESSTHAN,
    _ROW_TYPE_GREATERTHAN,
    _ROW_TYPE_INTERVAL,
    _ROW_TYPE_EQUAL_TO,
)

_row_type(::MOI.GreaterThan{Float64}) = _ROW_TYPE_GREATERTHAN
_row_type(::MOI.LessThan{Float64}) = _ROW_TYPE_LESSTHAN
_row_type(::MOI.EqualTo{Float64}) = _ROW_TYPE_EQUAL_TO
_row_type(::MOI.Interval{Float64}) = _ROW_TYPE_INTERVAL

_bounds(s::MOI.GreaterThan{Float64}) = s.lower, Inf
_bounds(s::MOI.LessThan{Float64}) = -Inf, s.upper
_bounds(s::MOI.EqualTo{Float64}) = s.value, s.value
_bounds(s::MOI.Interval{Float64}) = s.lower, s.upper

@enum(
    _BoundEnum,
    _BOUND_NONE,
    _BOUND_LESS_THAN,
    _BOUND_GREATER_THAN,
    _BOUND_LESS_AND_GREATER_THAN,
    _BOUND_INTERVAL,
    _BOUND_EQUAL_TO,
    _BOUND_SEMI_INTEGER,
    _BOUND_SEMI_CONTINUOUS,
)

@enum(
    _TypeEnum,
    _TYPE_CONTINUOUS,
    _TYPE_INTEGER,
    _TYPE_BINARY,
    _TYPE_SEMI_CONTINUOUS,
    _TYPE_SEMI_INTEGER,
)

const _SCALAR_SETS = Union{
    MOI.GreaterThan{Float64},
    MOI.LessThan{Float64},
    MOI.EqualTo{Float64},
    MOI.Interval{Float64},
}

"""
    _VariableInfo

A struct to store information about the variables.
"""
mutable struct _VariableInfo
    # We need to keep the index here because sometimes we call `LinearIndex`.
    index::MOI.VariableIndex
    # The variable name.
    name::String
    # The zero-indexed column in the HiGHS object.
    column::HighsInt
    # Storage to keep track of the variable bounds.
    bound::_BoundEnum
    lower::Float64
    upper::Float64
    # Track integrality
    type::_TypeEnum
    start::Union{Nothing,Float64}
    function _VariableInfo(
        index::MOI.VariableIndex,
        column::HighsInt,
        bound::_BoundEnum = _BOUND_NONE,
    )
        return new(
            index,
            "",
            column,
            bound,
            -Inf,
            Inf,
            _TYPE_CONTINUOUS,
            nothing,
        )
    end
end

function _update_info(info::_VariableInfo, s::MOI.GreaterThan{Float64})
    _throw_if_existing_lower(info, s)
    if info.bound == _BOUND_LESS_THAN
        info.bound = _BOUND_LESS_AND_GREATER_THAN
    else
        info.bound = _BOUND_GREATER_THAN
    end
    info.lower = s.lower
    return
end

function _update_info(info::_VariableInfo, s::MOI.LessThan{Float64})
    _throw_if_existing_upper(info, s)
    if info.bound == _BOUND_GREATER_THAN
        info.bound = _BOUND_LESS_AND_GREATER_THAN
    else
        info.bound = _BOUND_LESS_THAN
    end
    info.upper = s.upper
    return
end

function _update_info(info::_VariableInfo, s::MOI.EqualTo{Float64})
    _throw_if_existing_lower(info, s)
    _throw_if_existing_upper(info, s)
    info.bound = _BOUND_EQUAL_TO
    info.lower = s.value
    info.upper = s.value
    return
end

function _update_info(info::_VariableInfo, s::MOI.Interval{Float64})
    _throw_if_existing_lower(info, s)
    _throw_if_existing_upper(info, s)
    info.bound = _BOUND_INTERVAL
    info.lower = s.lower
    info.upper = s.upper
    return
end

function _variable_info_dict()
    return CleverDicts.CleverDict{MOI.VariableIndex,_VariableInfo}(
        x -> x.value,
        i -> MOI.VariableIndex(i),
    )
end

"""
    _ConstraintInfo

A struct to store information about the affine constraints.
"""
mutable struct _ConstraintInfo
    # The constraint name.
    name::String
    # The zero-indexed row in the HiGHS object.
    row::HighsInt
    # Storage to keep track of the constraint bounds.
    set::_RowType
    lower::Float64
    upper::Float64
end

function _ConstraintInfo(set::_SCALAR_SETS)
    lower, upper = _bounds(set)
    return _ConstraintInfo("", 0, _row_type(set), lower, upper)
end

struct _ConstraintKey
    value::Int64
end

function _constraint_info_dict()
    return CleverDicts.CleverDict{_ConstraintKey,_ConstraintInfo}(
        c -> c.value,
        i -> _ConstraintKey(i),
    )
end

"""
    _set(c::_ConstraintInfo)

Return the set associated with a constraint.
"""
function _set(c::_ConstraintInfo)
    if c.set == _ROW_TYPE_LESSTHAN
        return MOI.LessThan(c.upper)
    elseif c.set == _ROW_TYPE_GREATERTHAN
        return MOI.GreaterThan(c.lower)
    elseif c.set == _ROW_TYPE_INTERVAL
        return MOI.Interval(c.lower, c.upper)
    else
        @assert c.set == _ROW_TYPE_EQUAL_TO
        return MOI.EqualTo(c.lower)
    end
end

@enum(_OptimizeStatus, _OPTIMIZE_NOT_CALLED, _OPTIMIZE_OK, _OPTIMIZE_ERRORED)

"""
    _Solution

A struct to store the vector solution from HiGHS because it doesn't support
accessing them element-wise.
"""
mutable struct _Solution
    status::_OptimizeStatus
    model_status::HighsInt
    colvalue::Vector{Cdouble}
    coldual::Vector{Cdouble}
    colstatus::Vector{HighsInt}
    rowvalue::Vector{Cdouble}
    rowdual::Vector{Cdouble}
    rowstatus::Vector{HighsInt}
    primal_solution_status::HighsInt
    dual_solution_status::HighsInt
    has_primal_ray::Bool
    has_dual_ray::Bool
    function _Solution()
        return new(
            _OPTIMIZE_NOT_CALLED,
            kHighsModelStatusNotset,
            Cdouble[],
            Cdouble[],
            HighsInt[],
            Cdouble[],
            Cdouble[],
            HighsInt[],
            kHighsSolutionStatusNone,
            kHighsSolutionStatusNone,
            false,
            false,
        )
    end
end

function Base.empty!(x::_Solution)
    x.status = _OPTIMIZE_NOT_CALLED
    x.model_status = kHighsModelStatusNotset
    empty!(x.colvalue)
    empty!(x.coldual)
    empty!(x.colstatus)
    empty!(x.rowvalue)
    empty!(x.rowdual)
    empty!(x.rowstatus)
    x.primal_solution_status = kHighsSolutionStatusNone
    x.dual_solution_status = kHighsSolutionStatusNone
    x.has_dual_ray = false
    x.has_primal_ray = false
    return x
end

Base.isempty(x::_Solution) = x.status == _OPTIMIZE_NOT_CALLED

mutable struct _CallbackData
    f::Function
end

function Base.unsafe_convert(::Type{Ptr{Cvoid}}, d::_CallbackData)
    return pointer_from_objref(d)
end

"""
    Optimizer()

Create a new Optimizer object.
"""
mutable struct Optimizer <: MOI.AbstractOptimizer
    # A pointer to the underlying HiGHS optimizer.
    inner::Ptr{Cvoid}

    options::Dict{String,Any}

    # Storage for `MOI.Name`.
    name::String

    # A flag to keep track of MOI.FEASIBILITY_SENSE, since HiGHS only stores
    # MIN_SENSE or MAX_SENSE. This allows us to differentiate between MIN_SENSE
    # and FEASIBILITY_SENSE.
    is_feasibility::Bool
    is_objective_function_set::Bool
    is_objective_sense_set::Bool
    multi_objective::Union{Nothing,MOI.VectorAffineFunction{Float64}}

    # A flag to keep track of whether the objective is linear or quadratic.
    hessian::Union{Nothing,SparseArrays.SparseMatrixCSC{Float64,HighsInt}}

    # HiGHS doesn't have special support for binary variables. Cache them here
    # to modify bounds on solve.
    binaries::Set{_VariableInfo}

    variable_info::typeof(_variable_info_dict())
    affine_constraint_info::typeof(_constraint_info_dict())

    # Mappings from variable and constraint names to their indices. These are
    # lazily built on-demand, so most of the time, they are `nothing`.
    name_to_variable::Union{
        Nothing,
        Dict{String,Union{Nothing,MOI.VariableIndex}},
    }
    name_to_constraint_index::Union{
        Nothing,
        Dict{String,Union{Nothing,MOI.ConstraintIndex}},
    }

    # HiGHS just returns a single solution struct :(
    solution::_Solution
    callback_data::Union{Nothing,_CallbackData}

    compute_infeasibility_certificates::Bool

    function Optimizer()
        model = new(
            C_NULL,
            Dict{String,Any}(),
            "",
            true,
            false,
            false,
            nothing,
            nothing,
            Set{_VariableInfo}(),
            _variable_info_dict(),
            _constraint_info_dict(),
            nothing,
            nothing,
            _Solution(),
            nothing,
            true,
        )
        MOI.empty!(model)
        if Sys.iswindows()
            # Set number of threads by default on Windows. This fixes a hang in
            # some versions of Highs_jll. For discussion, see
            # https://github.com/ERGO-Code/HiGHS/issues/1257
            n_threads = max(1, div(length(Sys.cpu_info()), 2))
            MOI.set(model, MOI.NumberOfThreads(), n_threads)
        end
        finalizer(Highs_destroy, model)
        return model
    end
end

Base.cconvert(::Type{Ptr{Cvoid}}, model::Optimizer) = model
Base.unsafe_convert(::Type{Ptr{Cvoid}}, model::Optimizer) = model.inner

MOI.get(::Optimizer, ::MOI.SolverName) = "HiGHS"

function MOI.get(::Optimizer, ::MOI.SolverVersion)
    X, Y, Z = Highs_versionMajor(), Highs_versionMinor(), Highs_versionPatch()
    return "v$X.$Y.$Z"
end

function _check_ret(ret::HighsInt)
    if ret == kHighsStatusError
        error(
            "Encountered an error in HiGHS (Status $(ret)). Check the log " *
            "for details.",
        )
    end
    return
end

function Base.show(io::IO, model::Optimizer)
    return print(
        io,
        "A HiGHS model with $(Highs_getNumCols(model)) columns and " *
        "$(Highs_getNumRows(model)) rows.",
    )
end

function MOI.empty!(model::Optimizer)
    if model.inner != C_NULL
        Highs_destroy(model)
    end
    model.inner = Highs_create()
    for (key, value) in model.options
        MOI.set(model, MOI.RawOptimizerAttribute(key), value)
    end
    model.name = ""
    model.is_feasibility = true
    model.is_objective_function_set = false
    model.is_objective_sense_set = false
    model.multi_objective = nothing
    model.hessian = nothing
    empty!(model.binaries)
    empty!(model.variable_info)
    empty!(model.affine_constraint_info)
    model.name_to_variable = nothing
    model.name_to_constraint_index = nothing
    empty!(model.solution)
    ret = Highs_changeObjectiveOffset(model, 0.0)
    _check_ret(ret)
    model.callback_data = nothing
    return
end

function MOI.is_empty(model::Optimizer)
    offset = Ref{Cdouble}(0.0)
    ret = Highs_getObjectiveOffset(model, offset)
    _check_ret(ret)
    return Highs_getNumCols(model) == 0 &&
           Highs_getNumRows(model) == 0 &&
           model.is_feasibility &&
           !model.is_objective_function_set &&
           !model.is_objective_sense_set &&
           model.hessian === nothing &&
           isempty(model.binaries) &&
           isempty(model.variable_info) &&
           isempty(model.affine_constraint_info) &&
           model.name_to_variable === nothing &&
           model.name_to_constraint_index === nothing &&
           isempty(model.solution) &&
           iszero(offset[]) &&
           model.callback_data === nothing
end

MOI.get(model::Optimizer, ::MOI.RawSolver) = model

function MOI.get(model::Optimizer, ::MOI.ListOfVariableAttributesSet)
    ret = MOI.AbstractVariableAttribute[MOI.VariableName()]
    if any(info -> info.start !== nothing, values(model.variable_info))
        push!(ret, MOI.VariablePrimalStart())
    end
    return ret
end

function MOI.get(model::Optimizer, ::MOI.ListOfModelAttributesSet)
    attributes = MOI.AbstractModelAttribute[]
    if model.is_objective_sense_set
        push!(attributes, MOI.ObjectiveSense())
    end
    if model.is_objective_function_set
        F = MOI.get(model, MOI.ObjectiveFunctionType())
        push!(attributes, MOI.ObjectiveFunction{F}())
    end
    if !isempty(model.name)
        push!(attributes, MOI.Name())
    end
    return attributes
end

function MOI.get(
    ::Optimizer,
    ::MOI.ListOfConstraintAttributesSet{F,S},
) where {F,S}
    attributes = MOI.AbstractConstraintAttribute[]
    if F != MOI.VariableIndex
        return push!(attributes, MOI.ConstraintName())
    end
    return attributes
end

function MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{F,S}) where {F,S}
    # TODO: this could be more efficient.
    return length(MOI.get(model, MOI.ListOfConstraintIndices{F,S}()))
end

function MOI.get(model::Optimizer, ::MOI.ListOfConstraintTypesPresent)
    constraints = Set{Tuple{Type,Type}}()
    for info in values(model.variable_info)
        if info.bound == _BOUND_NONE
        elseif info.bound == _BOUND_LESS_THAN
            push!(constraints, (MOI.VariableIndex, MOI.LessThan{Float64}))
        elseif info.bound == _BOUND_GREATER_THAN
            push!(constraints, (MOI.VariableIndex, MOI.GreaterThan{Float64}))
        elseif info.bound == _BOUND_LESS_AND_GREATER_THAN
            push!(constraints, (MOI.VariableIndex, MOI.LessThan{Float64}))
            push!(constraints, (MOI.VariableIndex, MOI.GreaterThan{Float64}))
        elseif info.bound == _BOUND_EQUAL_TO
            push!(constraints, (MOI.VariableIndex, MOI.EqualTo{Float64}))
        elseif info.bound == _BOUND_INTERVAL
            push!(constraints, (MOI.VariableIndex, MOI.Interval{Float64}))
        elseif info.bound == _BOUND_SEMI_CONTINUOUS
            push!(constraints, (MOI.VariableIndex, MOI.Semicontinuous{Float64}))
        else
            @assert info.bound == _BOUND_SEMI_INTEGER
            push!(constraints, (MOI.VariableIndex, MOI.Semiinteger{Float64}))
        end
        if info.type == _TYPE_INTEGER
            push!(constraints, (MOI.VariableIndex, MOI.Integer))
        elseif info.type == _TYPE_BINARY
            push!(constraints, (MOI.VariableIndex, MOI.ZeroOne))
        end
    end
    for info in values(model.affine_constraint_info)
        push!(
            constraints,
            (MOI.ScalarAffineFunction{Float64}, typeof(_set(info))),
        )
    end
    return collect(constraints)
end

###
### HiGHS.ComputeInfeasibilityCertificate
###

"""
    HiGHS.ComputeInfeasibilityCertificate()

An `MOI.AbstractOptimizerAttribute` to control whether HiGHS attempts to
compute an infeasibility certificate if the primal is primal or dual infeasible.

The default is `true`. Set to `false` to disable.
"""
struct ComputeInfeasibilityCertificate <: MOI.AbstractOptimizerAttribute end

MOI.supports(::Optimizer, ::ComputeInfeasibilityCertificate) = true

function MOI.get(model::Optimizer, ::ComputeInfeasibilityCertificate)
    return model.compute_infeasibility_certificates
end

function MOI.set(model::Optimizer, ::ComputeInfeasibilityCertificate, v::Bool)
    model.compute_infeasibility_certificates = v
    return
end

###
### MOI.RawOptimizerAttribute
###

function MOI.supports(model::Optimizer, param::MOI.RawOptimizerAttribute)
    typeP = Ref{HighsInt}()
    return Highs_getOptionType(model, param.name, typeP) == 0
end

function _check_option_status(ret::HighsInt)
    if ret != 0
        error("Encountered an error in HiGHS: Check the log for details.")
    end
    return
end

_highs_option_type(::Bool) = kHighsOptionTypeBool
_highs_option_type(::Integer) = kHighsOptionTypeInt
_highs_option_type(::AbstractFloat) = kHighsOptionTypeDouble
_highs_option_type(::String) = kHighsOptionTypeString
_highs_option_type(::Any) = HighsInt(-1)

function _type_for_highs_option(k)
    if k == kHighsOptionTypeBool
        return Bool
    elseif k == kHighsOptionTypeInt
        return Integer
    elseif k == kHighsOptionTypeDouble
        return AbstractFloat
    else
        @assert k == 3
        return String
    end
end

function _set_option(model::Optimizer, option::String, value::Bool)
    return Highs_setBoolOptionValue(model, option, HighsInt(value))
end

function _set_option(model::Optimizer, option::String, value::Integer)
    return Highs_setIntOptionValue(model, option, HighsInt(value))
end

function _set_option(model::Optimizer, option::String, value::AbstractFloat)
    return Highs_setDoubleOptionValue(model, option, Cdouble(value))
end

function _set_option(model::Optimizer, option::String, value::String)
    return Highs_setStringOptionValue(model, option, value)
end

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    typeP = Ref{HighsInt}()
    ret = Highs_getOptionType(model, param.name, typeP)
    if ret != 0
        throw(MOI.UnsupportedAttribute(param))
    elseif typeP[] != _highs_option_type(value)
        throw(
            MOI.SetAttributeNotAllowed(
                param,
                "\n\nInvalid value `$(value)::$(typeof(value))` for option " *
                "\"$(param.name)\", expected a value of type " *
                "`$(_type_for_highs_option(typeP[]))`.\n\n",
            ),
        )
    end
    model.options[param.name] = value
    ret = _set_option(model, param.name, value)
    return _check_option_status(ret)
end

### MOI.get

function _get_bool_option(model::Optimizer, option::String)
    value = Ref{HighsInt}(0)
    ret = Highs_getBoolOptionValue(model, option, value)
    _check_option_status(ret)
    return Bool(value[])
end

function _get_int_option(model::Optimizer, option::String)
    value = Ref{HighsInt}(0)
    ret = Highs_getIntOptionValue(model, option, value)
    _check_option_status(ret)
    return value[]
end

function _get_double_option(model::Optimizer, option::String)
    value = Ref{Cdouble}()
    ret = Highs_getDoubleOptionValue(model, option, value)
    _check_option_status(ret)
    return value[]
end

function _get_string_option(model::Optimizer, option::String)
    buffer = Vector{Cchar}(undef, kHighsMaximumStringLength)
    bufferP = pointer(buffer)
    GC.@preserve buffer begin
        ret = Highs_getStringOptionValue(model, option, bufferP)
        _check_option_status(ret)
        return unsafe_string(bufferP)
    end
end

function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    typeP = Ref{HighsInt}()
    ret = Highs_getOptionType(model, param.name, typeP)
    if ret != 0
        throw(MOI.UnsupportedAttribute(param))
    elseif typeP[] == 0
        return _get_bool_option(model, param.name)::Bool
    elseif typeP[] == 1
        return _get_int_option(model, param.name)::HighsInt
    elseif typeP[] == 2
        return _get_double_option(model, param.name)::Cdouble
    else
        @assert typeP[] == 3
        return _get_string_option(model, param.name)::String
    end
end

###
### MOI.TimeLimitSec
###

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, ::Nothing)
    return MOI.set(model, MOI.RawOptimizerAttribute("time_limit"), Inf)
end

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, limit::Real)
    return MOI.set(
        model,
        MOI.RawOptimizerAttribute("time_limit"),
        Float64(limit),
    )
end

function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    value = MOI.get(model, MOI.RawOptimizerAttribute("time_limit"))
    return value == Inf ? nothing : value
end

###
### MOI.NodeLimit
###

MOI.supports(::Optimizer, ::MOI.NodeLimit) = true

function MOI.set(model::Optimizer, ::MOI.NodeLimit, ::Nothing)
    attr = MOI.RawOptimizerAttribute("mip_max_nodes")
    MOI.set(model, attr, typemax(Cint))
    return
end

function MOI.set(model::Optimizer, ::MOI.NodeLimit, limit::Real)
    attr = MOI.RawOptimizerAttribute("mip_max_nodes")
    MOI.set(model, attr, Int(limit))
    return
end

function MOI.get(model::Optimizer, ::MOI.NodeLimit)
    value = MOI.get(model, MOI.RawOptimizerAttribute("mip_max_nodes"))
    return value == typemax(Cint) ? nothing : Int(value)
end

###
### MOI.AbsoluteGapTolerance
###

MOI.supports(::Optimizer, ::MOI.AbsoluteGapTolerance) = true

function MOI.set(model::Optimizer, ::MOI.AbsoluteGapTolerance, limit::Real)
    return MOI.set(
        model,
        MOI.RawOptimizerAttribute("mip_abs_gap"),
        Float64(limit),
    )
end

function MOI.get(model::Optimizer, ::MOI.AbsoluteGapTolerance)
    return MOI.get(model, MOI.RawOptimizerAttribute("mip_abs_gap"))
end

###
### MOI.RelativeGapTolerance
###

MOI.supports(::Optimizer, ::MOI.RelativeGapTolerance) = true

function MOI.set(model::Optimizer, ::MOI.RelativeGapTolerance, limit::Real)
    return MOI.set(
        model,
        MOI.RawOptimizerAttribute("mip_rel_gap"),
        Float64(limit),
    )
end

function MOI.get(model::Optimizer, ::MOI.RelativeGapTolerance)
    return MOI.get(model, MOI.RawOptimizerAttribute("mip_rel_gap"))
end

###
### MOI.ObjectiveLimit
###

MOI.supports(::Optimizer, ::MOI.ObjectiveLimit) = true

function MOI.set(model::Optimizer, ::MOI.ObjectiveLimit, limit::Real)
    attr = MOI.RawOptimizerAttribute("objective_target")
    MOI.set(model, attr, Float64(limit))
    return
end

function MOI.set(model::Optimizer, ::MOI.ObjectiveLimit, ::Nothing)
    MOI.set(model, MOI.ObjectiveLimit(), -Inf)
    return
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveLimit)
    value = MOI.get(model, MOI.RawOptimizerAttribute("objective_target"))
    return value == -Inf ? nothing : value
end

###
### MOI.Silent
###

MOI.supports(::Optimizer, ::MOI.Silent) = true

function MOI.get(model::Optimizer, ::MOI.Silent)
    return !MOI.get(model, MOI.RawOptimizerAttribute("output_flag"))
end

function MOI.set(model::Optimizer, ::MOI.Silent, flag::Bool)
    MOI.set(model, MOI.RawOptimizerAttribute("output_flag"), !flag)
    return
end

###
### MOI.Name
###

MOI.supports(::Optimizer, ::MOI.Name) = true

MOI.get(model::Optimizer, ::MOI.Name) = model.name

MOI.set(model::Optimizer, ::MOI.Name, name::String) = (model.name = name)

###
### MOI.NumberOfThreads
###

MOI.supports(::Optimizer, ::MOI.NumberOfThreads) = true

function MOI.get(model::Optimizer, ::MOI.NumberOfThreads)::Int
    return MOI.get(model, MOI.RawOptimizerAttribute("threads"))
end

function MOI.set(model::Optimizer, ::MOI.NumberOfThreads, threads::Int)
    MOI.set(model, MOI.RawOptimizerAttribute("threads"), threads)
    return
end

###
### Variables
###

function MOI.get(model::Optimizer, ::MOI.NumberOfVariables)
    return length(model.variable_info)
end

function MOI.get(model::Optimizer, ::MOI.ListOfVariableIndices)
    return sort!(collect(keys(model.variable_info)), by = x -> x.value)
end

function _info(model::Optimizer, key::MOI.VariableIndex)
    info = get(model.variable_info, key, nothing)
    if info === nothing
        return throw(MOI.InvalidIndex(key))
    end
    return info
end

"""
    column(model::Optimizer, x::MOI.VariableIndex)

Return the 0-indexed column associated with `x` in `model`.
"""
column(model::Optimizer, x::MOI.VariableIndex) = _info(model, x).column

function MOI.add_variable(model::Optimizer)
    # Initialize `_VariableInfo` with a dummy `VariableIndex` and a column,
    # because we need `add_item` to tell us what the `VariableIndex` is.
    index = CleverDicts.add_item(
        model.variable_info,
        _VariableInfo(MOI.VariableIndex(0), HighsInt(0)),
    )
    info = _info(model, index)
    # Now, set `.index` and `.column`.
    info.index = index
    info.column = HighsInt(length(model.variable_info) - 1)
    ret = Highs_addCol(model, 0.0, -Inf, Inf, 0, C_NULL, C_NULL)
    _check_ret(ret)
    return index
end

function MOI.supports_add_constrained_variable(
    ::Optimizer,
    ::Type{
        <:Union{
            MOI.GreaterThan{Float64},
            MOI.LessThan{Float64},
            MOI.EqualTo{Float64},
            MOI.Interval{Float64},
        },
    },
)
    return true
end

function MOI.add_constrained_variable(
    model::Optimizer,
    set::S,
) where {
    S<:Union{
        MOI.GreaterThan{Float64},
        MOI.LessThan{Float64},
        MOI.EqualTo{Float64},
        MOI.Interval{Float64},
    },
}
    # Initialize `_VariableInfo` with a dummy `VariableIndex` and a column,
    # because we need `add_item` to tell us what the `VariableIndex` is.
    index = CleverDicts.add_item(
        model.variable_info,
        _VariableInfo(MOI.VariableIndex(0), HighsInt(0)),
    )
    info = _info(model, index)
    # Now, set `.index` and `.column`.
    info.index = index
    info.column = HighsInt(length(model.variable_info) - 1)
    lower, upper = _bounds(set)
    ret = Highs_addCol(model, 0.0, lower, upper, 0, C_NULL, C_NULL)
    _check_ret(ret)
    _update_info(info, set)
    return index, MOI.ConstraintIndex{MOI.VariableIndex,S}(index.value)
end

function MOI.supports_add_constrained_variable(
    ::Optimizer,
    ::Type{Tuple{MOI.GreaterThan{Float64},MOI.LessThan{Float64}}},
)
    return true
end

function MOI.add_constrained_variable(
    model::Optimizer,
    set::Tuple{MOI.GreaterThan{Float64},MOI.LessThan{Float64}},
)
    # Initialize `_VariableInfo` with a dummy `VariableIndex` and a column,
    # because we need `add_item` to tell us what the `VariableIndex` is.
    index = CleverDicts.add_item(
        model.variable_info,
        _VariableInfo(MOI.VariableIndex(0), HighsInt(0)),
    )
    info = _info(model, index)
    # Now, set `.index` and `.column`.
    info.index = index
    info.column = HighsInt(length(model.variable_info) - 1)
    l, u = set[1].lower, set[2].upper
    ret = Highs_addCol(model, 0.0, l, u, 0, C_NULL, C_NULL)
    _check_ret(ret)
    _update_info(info, set[1])
    _update_info(info, set[2])
    c_gt = MOI.ConstraintIndex{MOI.VariableIndex,(typeof(set[1]))}(index.value)
    c_lt = MOI.ConstraintIndex{MOI.VariableIndex,(typeof(set[2]))}(index.value)
    return index, (c_gt, c_lt)
end

function MOI.is_valid(model::Optimizer, v::MOI.VariableIndex)
    return haskey(model.variable_info, v)
end

function MOI.delete(model::Optimizer, v::MOI.VariableIndex)
    col = column(model, v)
    ret = Highs_deleteColsByRange(model, col, col)
    _check_ret(ret)
    delete!(model.variable_info, v)
    for other_info in values(model.variable_info)
        if other_info.column > col
            other_info.column -= 1
        end
    end
    if model.multi_objective !== nothing
        model.multi_objective =
            MOI.Utilities.filter_variables(!=(v), model.multi_objective)
    end
    model.name_to_variable = nothing
    model.name_to_constraint_index = nothing
    return
end

#
# Variable names
#

function MOI.supports(
    ::Optimizer,
    ::MOI.VariableName,
    ::Type{MOI.VariableIndex},
)
    return true
end

function MOI.get(model::Optimizer, ::Type{MOI.VariableIndex}, name::String)
    if model.name_to_variable === nothing
        _rebuild_name_to_variable(model)
    end
    variable = get(model.name_to_variable, name, missing)
    if variable === nothing
        error("Duplicate name detected: $(name)")
    elseif variable === missing
        return nothing  # Variable with this name does not exist
    end
    return variable::MOI.VariableIndex
end

function MOI.get(model::Optimizer, ::MOI.VariableName, v::MOI.VariableIndex)
    return _info(model, v).name
end

function MOI.set(
    model::Optimizer,
    ::MOI.VariableName,
    v::MOI.VariableIndex,
    name::String,
)
    info = _info(model, v)
    if info.name == name
        return  # The name is unchanged
    end
    info.name = name
    model.name_to_variable = nothing
    if isempty(name)
        # HiGHS does not support passing "" or C_NULL for names. So we default
        # to a fairly standard C$i for the column name.
        name = string("C", info.column)
    end
    ret = Highs_passColName(model, info.column, name)
    _check_ret(ret)
    return
end

function _rebuild_name_to_variable(model::Optimizer)
    model.name_to_variable = Dict{String,Union{Nothing,MOI.VariableIndex}}()
    for (index, info) in model.variable_info
        if isempty(info.name)
            continue
        end
        if haskey(model.name_to_variable, info.name)
            model.name_to_variable[info.name] = nothing
        else
            model.name_to_variable[info.name] = index
        end
    end
    return
end

#
# MOI.VariablePrimalStart
#

function MOI.supports(
    ::Optimizer,
    ::MOI.VariablePrimalStart,
    ::Type{MOI.VariableIndex},
)
    return true
end

function MOI.set(
    model::Optimizer,
    ::MOI.VariablePrimalStart,
    x::MOI.VariableIndex,
    start::Union{Nothing,Float64},
)
    info = _info(model, x)
    info.start = start
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.VariablePrimalStart,
    x::MOI.VariableIndex,
)
    return _info(model, x).start
end

###
### Objectives
###

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveSense,
    sense::MOI.OptimizationSense,
)
    x = sense == MOI.MAX_SENSE ? kHighsObjSenseMaximize : kHighsObjSenseMinimize
    ret = Highs_changeObjectiveSense(model, x)
    _check_ret(ret)
    if sense == MOI.FEASIBILITY_SENSE
        model.is_feasibility = true
        n = MOI.get(model, MOI.NumberOfVariables())
        ret = Highs_changeColsCostByRange(
            model,
            HighsInt(0),
            HighsInt(n - 1),
            zeros(Cdouble, n),
        )
        _check_ret(ret)
        ret = Highs_changeObjectiveOffset(model, 0.0)
        _check_ret(ret)
    else
        model.is_feasibility = false
    end
    model.is_objective_sense_set = true
    return
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveSense)
    if model.is_feasibility
        return MOI.FEASIBILITY_SENSE
    end
    senseP = Ref{HighsInt}()
    ret = Highs_getObjectiveSense(model, senseP)
    _check_ret(ret)
    return senseP[] == kHighsObjSenseMinimize ? MOI.MIN_SENSE : MOI.MAX_SENSE
end

function MOI.supports(
    ::Optimizer,
    ::Union{
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
        MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}},
    },
)
    return true
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveFunctionType)
    if model.multi_objective !== nothing
        return MOI.VectorAffineFunction{Float64}
    elseif model.hessian === nothing
        return MOI.ScalarAffineFunction{Float64}
    else
        return MOI.ScalarQuadraticFunction{Float64}
    end
end

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
    f::MOI.ScalarAffineFunction{Float64},
)
    if model.multi_objective !== nothing
        ret = Highs_clearLinearObjectives(model)
        _check_ret(ret)
        model.multi_objective = nothing
    end
    num_vars = HighsInt(length(model.variable_info))
    obj = zeros(Float64, num_vars)
    for term in f.terms
        col = column(model, term.variable)
        obj[col+1] += term.coefficient
    end
    ret = Highs_changeColsCostByRange(model, HighsInt(0), num_vars - 1, obj)
    _check_ret(ret)
    ret = Highs_changeObjectiveOffset(model, f.constant)
    _check_ret(ret)
    model.is_objective_function_set = true
    if model.hessian !== nothing
        start = zeros(HighsInt, num_vars)
        ret = Highs_passHessian(
            model,
            num_vars,
            0,
            kHighsHessianFormatTriangular,
            start,
            C_NULL,
            C_NULL,
        )
        _check_ret(ret)
        model.hessian = nothing
    end
    return
end

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}},
    f::MOI.ScalarQuadraticFunction{Float64},
)
    if model.multi_objective !== nothing
        ret = Highs_clearLinearObjectives(model)
        _check_ret(ret)
        model.multi_objective = nothing
    end
    numcol = length(model.variable_info)
    obj = zeros(Float64, numcol)
    for term in f.affine_terms
        obj[column(model, term.variable)+1] += term.coefficient
    end
    ret = Highs_changeColsCostByMask(model, ones(HighsInt, numcol), obj)
    _check_ret(ret)
    ret = Highs_changeObjectiveOffset(model, f.constant)
    _check_ret(ret)
    I, J, V = HighsInt[], HighsInt[], Cdouble[]
    for term in f.quadratic_terms
        if iszero(term.coefficient)
            continue
        end
        i = model.variable_info[term.variable_1].column
        j = model.variable_info[term.variable_2].column
        row, col = (i > j) ? (i, j) : (j, i)
        push!(I, row + 1)
        push!(J, col + 1)
        push!(V, term.coefficient)
    end
    Q = SparseArrays.sparse(I, J, V, numcol, numcol)
    ret = Highs_passHessian(
        model,
        numcol,
        length(Q.nzval),
        kHighsHessianFormatTriangular,
        Q.colptr .- HighsInt(1),
        Q.rowval .- HighsInt(1),
        Q.nzval,
    )
    _check_ret(ret)
    model.hessian = Q
    model.is_objective_function_set = true
    return
end

function _get_objective_function(model::Optimizer, ::Nothing)
    numcol = MOI.get(model, MOI.NumberOfVariables())
    offset = Ref{Cdouble}(0.0)
    ret = Highs_getObjectiveOffset(model, offset)
    _check_ret(ret)
    if numcol == 0
        return MOI.ScalarAffineFunction{Float64}(
            MOI.ScalarAffineTerm{Float64}[],
            offset[],
        )
    end
    num_colsP, nnzP = Ref{HighsInt}(0), Ref{HighsInt}(0)
    costs = Vector{Cdouble}(undef, numcol)
    ret = Highs_getColsByRange(
        model,
        0,
        numcol - 1,
        num_colsP,
        costs,
        C_NULL,
        C_NULL,
        nnzP,
        C_NULL,
        C_NULL,
        C_NULL,
    )
    return MOI.ScalarAffineFunction{Float64}(
        MOI.ScalarAffineTerm{Float64}[
            MOI.ScalarAffineTerm(costs[info.column+1], index) for
            (index, info) in model.variable_info if
            !iszero(costs[info.column+1])
        ],
        offset[],
    )
end

function _get_objective_function(model::Optimizer, Q)
    f = _get_objective_function(model, nothing)
    return MOI.ScalarQuadraticFunction(
        MOI.ScalarQuadraticTerm{Float64}[
            MOI.ScalarQuadraticTerm(
                Q.nzval[j],
                model.variable_info[CleverDicts.LinearIndex(col)].index,
                model.variable_info[CleverDicts.LinearIndex(Q.rowval[j])].index,
            ) for col in 1:Q.m for j in Q.colptr[col]:(Q.colptr[col+1]-1)
        ],
        f.terms,
        f.constant,
    )
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveFunction{F}) where {F}
    return convert(F, _get_objective_function(model, model.hessian))
end

function MOI.modify(
    model::Optimizer,
    ::MOI.ObjectiveFunction{
        <:Union{
            MOI.ScalarAffineFunction{Float64},
            MOI.ScalarQuadraticFunction{Float64},
        },
    },
    chg::MOI.ScalarConstantChange{Float64},
)
    ret = Highs_changeObjectiveOffset(model, chg.new_constant)
    _check_ret(ret)
    model.is_objective_function_set = true
    return
end

function MOI.modify(
    model::Optimizer,
    ::MOI.ObjectiveFunction{
        <:Union{
            MOI.ScalarAffineFunction{Float64},
            MOI.ScalarQuadraticFunction{Float64},
        },
    },
    chg::MOI.ScalarCoefficientChange{Float64},
)
    ret = Highs_changeColCost(
        model,
        column(model, chg.variable),
        chg.new_coefficient,
    )
    _check_ret(ret)
    model.is_objective_function_set = true
    return
end

function MOI.modify(
    model::Optimizer,
    attr::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}},
    chg::MOI.ScalarQuadraticCoefficientChange{Float64},
)
    if model.hessian === nothing
        f = MOI.get(
            model,
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        )
        term = MOI.ScalarQuadraticTerm(
            chg.new_coefficient,
            chg.variable_1,
            chg.variable_2,
        )
        new_f = MOI.ScalarQuadraticFunction([term], f.terms, f.constant)
        MOI.set(model, attr, new_f)
        return
    end
    i, j = column(model, chg.variable_1), column(model, chg.variable_2)
    if i < j
        j, i = i, j
    end
    model.hessian[i+1, j+1] = chg.new_coefficient
    ret = Highs_passHessian(
        model,
        size(model.hessian, 1),
        length(model.hessian.nzval),
        kHighsHessianFormatTriangular,
        model.hessian.colptr .- HighsInt(1),
        model.hessian.rowval .- HighsInt(1),
        model.hessian.nzval,
    )
    _check_ret(ret)
    model.is_objective_function_set = true
    return
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{MOI.VectorAffineFunction{Float64}},
)
    return true
end

function MOI.get(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.VectorAffineFunction{Float64}},
)
    return model.multi_objective
end

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.VectorAffineFunction{Float64}},
    f::MOI.VectorAffineFunction{Float64},
)
    num_vars = HighsInt(length(model.variable_info))
    O = MOI.output_dimension(f)
    obj_coefs = zeros(Float64, O * num_vars)
    for term in f.terms
        col = column(model, term.scalar_term.variable) + 1
        obj_coefs[O*(term.output_index-1)+col] += term.scalar_term.coefficient
    end
    # senseP will be 1 if MIN and -1 if MAX
    senseP = Ref{HighsInt}()
    ret = Highs_getObjectiveSense(model, senseP)
    _check_ret(ret)
    ret = Highs_passLinearObjectives(
        model,
        O,                           # num_linear_objective
        fill(Float64(senseP[]), O),  # weight: set to -1 if maximizing
        f.constants,                 # offset
        obj_coefs,                   # coefficients
        zeros(Float64, O),           # abs_tolerance
        zeros(Float64, O),           # rel_tolerance
        ones(Cint, O),               # priority
    )
    _check_ret(ret)
    model.multi_objective = f
    model.is_objective_function_set = true
    return
end

###
### VariableIndex-in-Set constraints.
###

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VariableIndex},
    ::Type{<:_SCALAR_SETS},
)
    return true
end

function _bound_enums(::Type{MOI.LessThan{Float64}})
    return (_BOUND_LESS_THAN, _BOUND_LESS_AND_GREATER_THAN)
end

function _bound_enums(::Type{MOI.GreaterThan{Float64}})
    return (_BOUND_GREATER_THAN, _BOUND_LESS_AND_GREATER_THAN)
end

_bound_enums(::Type{MOI.Interval{Float64}}) = (_BOUND_INTERVAL,)

_bound_enums(::Type{MOI.EqualTo{Float64}}) = (_BOUND_EQUAL_TO,)

function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintIndices{MOI.VariableIndex,S},
) where {S<:_SCALAR_SETS}
    indices = MOI.ConstraintIndex{MOI.VariableIndex,S}[
        MOI.ConstraintIndex{MOI.VariableIndex,S}(key.value) for
        (key, info) in model.variable_info if info.bound in _bound_enums(S)
    ]
    return sort!(indices, by = x -> x.value)
end

function _info(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,<:Any},
)
    var_index = MOI.VariableIndex(c.value)
    info = get(model.variable_info, var_index, nothing)
    if info === nothing
        return throw(MOI.InvalidIndex(c))
    end
    return info
end

"""
    column(
        model::Optimizer,
        c::MOI.ConstraintIndex{MOI.VariableIndex,<:Any},
    )

Return the 0-indexed column associated with the variable bounds `c` in `model`.
"""
function column(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,<:Any},
)
    return _info(model, c).column
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}},
)
    info = get(model.variable_info, MOI.VariableIndex(c.value), nothing)
    if info === nothing
        return false
    end
    return info.bound == _BOUND_LESS_THAN ||
           info.bound == _BOUND_LESS_AND_GREATER_THAN
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}},
)
    info = get(model.variable_info, MOI.VariableIndex(c.value), nothing)
    if info === nothing
        return false
    end
    return info.bound == _BOUND_GREATER_THAN ||
           info.bound == _BOUND_LESS_AND_GREATER_THAN
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Interval{Float64}},
)
    info = get(model.variable_info, MOI.VariableIndex(c.value), nothing)
    if info === nothing
        return false
    end
    return info.bound == _BOUND_INTERVAL
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.EqualTo{Float64}},
)
    info = get(model.variable_info, MOI.VariableIndex(c.value), nothing)
    if info === nothing
        return false
    end
    return info.bound == _BOUND_EQUAL_TO
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.VariableIndex,<:Any},
)
    MOI.throw_if_not_valid(model, c)
    return MOI.VariableIndex(c.value)
end

function MOI.set(
    ::Optimizer,
    ::MOI.ConstraintFunction,
    ::MOI.ConstraintIndex{MOI.VariableIndex,<:Any},
    ::MOI.VariableIndex,
)
    return throw(MOI.SettingVariableIndexNotAllowed())
end

function _throw_if_existing_lower(
    info::_VariableInfo,
    ::S,
) where {S<:MOI.AbstractSet}
    if info.bound == _BOUND_LESS_AND_GREATER_THAN
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Float64},S}(info.index))
    elseif info.bound == _BOUND_GREATER_THAN
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Float64},S}(info.index))
    elseif info.bound == _BOUND_INTERVAL
        throw(MOI.LowerBoundAlreadySet{MOI.Interval{Float64},S}(info.index))
    elseif info.bound == _BOUND_EQUAL_TO
        throw(MOI.LowerBoundAlreadySet{MOI.EqualTo{Float64},S}(info.index))
    elseif info.bound == _BOUND_SEMI_CONTINUOUS
        T = MOI.Semicontinuous{Float64}
        throw(MOI.LowerBoundAlreadySet{T,S}(info.index))
    elseif info.bound == _BOUND_SEMI_INTEGER
        throw(MOI.LowerBoundAlreadySet{MOI.Semiinteger{Float64},S}(info.index))
    end
    return
end

function _throw_if_existing_upper(
    info::_VariableInfo,
    ::S,
) where {S<:MOI.AbstractSet}
    if info.bound == _BOUND_LESS_AND_GREATER_THAN
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Float64},S}(info.index))
    elseif info.bound == _BOUND_LESS_THAN
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Float64},S}(info.index))
    elseif info.bound == _BOUND_INTERVAL
        throw(MOI.UpperBoundAlreadySet{MOI.Interval{Float64},S}(info.index))
    elseif info.bound == _BOUND_EQUAL_TO
        throw(MOI.UpperBoundAlreadySet{MOI.EqualTo{Float64},S}(info.index))
    elseif info.bound == _BOUND_SEMI_CONTINUOUS
        T = MOI.Semicontinuous{Float64}
        throw(MOI.UpperBoundAlreadySet{T,S}(info.index))
    elseif info.bound == _BOUND_SEMI_INTEGER
        throw(MOI.UpperBoundAlreadySet{MOI.Semiinteger{Float64},S}(info.index))
    end
    return
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.VariableIndex,
    s::S,
) where {S<:_SCALAR_SETS}
    info = _info(model, f)
    _update_info(info, s)
    index = MOI.ConstraintIndex{MOI.VariableIndex,typeof(s)}(f.value)
    MOI.set(model, MOI.ConstraintSet(), index, s)
    return index
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    ret = Highs_changeColBounds(model, info.column, info.lower, Inf)
    _check_ret(ret)
    info.upper = Inf
    if info.bound == _BOUND_LESS_AND_GREATER_THAN
        info.bound = _BOUND_GREATER_THAN
    else
        info.bound = _BOUND_NONE
    end
    model.name_to_constraint_index = nothing
    return
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    ret = Highs_changeColBounds(model, info.column, -Inf, info.upper)
    _check_ret(ret)
    info.lower = -Inf
    if info.bound == _BOUND_LESS_AND_GREATER_THAN
        info.bound = _BOUND_LESS_THAN
    else
        info.bound = _BOUND_NONE
    end
    model.name_to_constraint_index = nothing
    return
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Interval{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    ret = Highs_changeColBounds(model, info.column, -Inf, Inf)
    _check_ret(ret)
    info.lower, info.upper = -Inf, Inf
    info.bound = _BOUND_NONE
    model.name_to_constraint_index = nothing
    return
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.EqualTo{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    ret = Highs_changeColBounds(model, info.column, -Inf, Inf)
    _check_ret(ret)
    info.lower, info.upper = -Inf, Inf
    info.bound = _BOUND_NONE
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    return MOI.GreaterThan(_info(model, c).lower)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    return MOI.LessThan(_info(model, c).upper)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.EqualTo{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    return MOI.EqualTo(_info(model, c).lower)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Interval{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    return MOI.Interval(info.lower, info.upper)
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex,S},
    s::S,
) where {S<:_SCALAR_SETS}
    MOI.throw_if_not_valid(model, c)
    lower, upper = _bounds(s)
    info = _info(model, c)
    if S == MOI.LessThan{Float64}
        ret = Highs_changeColBounds(model, info.column, info.lower, upper)
        _check_ret(ret)
        info.upper = upper
    elseif S == MOI.GreaterThan{Float64}
        ret = Highs_changeColBounds(model, info.column, lower, info.upper)
        _check_ret(ret)
        info.lower = lower
    else
        ret = Highs_changeColBounds(model, info.column, lower, upper)
        _check_ret(ret)
        info.lower = lower
        info.upper = upper
    end
    return
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::Vector{MOI.ConstraintIndex{MOI.VariableIndex,S}},
    s::Vector{S},
) where {S<:_SCALAR_SETS}
    if length(c) != length(s)
        msg = "number of constraints does not match number of sets"
        throw(DimensionMismatch(msg))
    end
    for ci in c
        MOI.throw_if_not_valid(model, ci)
    end
    N = length(c)
    columns, lower, upper = zeros(Cint, N), zeros(Cdouble, N), zeros(Cdouble, N)
    for (i, (ci, si)) in enumerate(zip(c, s))
        info = _info(model, ci)
        columns[i] = info.column
        l, u = _bounds(si)
        if S == MOI.LessThan{Float64}
            info.upper = u
            lower[i], upper[i] = info.lower, u
        elseif S == MOI.GreaterThan{Float64}
            info.lower = l
            lower[i], upper[i] = l, info.upper
        else
            lower[i], upper[i] = l, u
            info.lower, info.upper = l, u
        end
    end
    ret = Highs_changeColsBoundsBySet(model, N, columns, lower, upper)
    _check_ret(ret)
    return
end

###
### ScalarAffineFunction-in-Set
###

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarAffineFunction{Float64}},
    ::Type{<:_SCALAR_SETS},
)
    return true
end

function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64},S},
) where {S<:_SCALAR_SETS}
    indices = MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S}[
        MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S}(key.value)
        for (key, info) in model.affine_constraint_info if _set(info) isa S
    ]
    return sort!(indices; by = x -> x.value)
end

function _info(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:_SCALAR_SETS},
)
    key = _ConstraintKey(c.value)
    info = get(model.affine_constraint_info, key, nothing)
    if info === nothing
        return throw(MOI.InvalidIndex(c))
    end
    return info
end

function row(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:_SCALAR_SETS},
)
    return _info(model, c).row
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
) where {S<:_SCALAR_SETS}
    key = _ConstraintKey(c.value)
    info = get(model.affine_constraint_info, key, nothing)
    if info === nothing
        return false
    end
    return _set(info) isa S
end

function _indices_and_coefficients(
    indices::Vector{HighsInt},
    coefficients::Vector{Float64},
    model::Optimizer,
    f::MOI.ScalarAffineFunction{Float64},
)
    i = 1
    for term in f.terms
        indices[i] = column(model, term.variable)
        coefficients[i] = term.coefficient
        i += 1
    end
    return indices, coefficients
end

function _indices_and_coefficients(
    model::Optimizer,
    f::MOI.ScalarAffineFunction{Float64},
)
    f_canon = MOI.Utilities.canonical(f)
    nnz = length(f_canon.terms)
    indices, coefficients = zeros(HighsInt, nnz), zeros(Cdouble, nnz)
    _indices_and_coefficients(indices, coefficients, model, f_canon)
    return indices, coefficients
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.ScalarAffineFunction{Float64},
    s::_SCALAR_SETS,
)
    if !iszero(f.constant)
        throw(
            MOI.ScalarFunctionConstantNotZero{Float64,typeof(f),typeof(s)}(
                f.constant,
            ),
        )
    end
    key = CleverDicts.add_item(model.affine_constraint_info, _ConstraintInfo(s))
    model.affine_constraint_info[key].row =
        HighsInt(length(model.affine_constraint_info) - 1)
    indices, coefficients = _indices_and_coefficients(model, f)
    lower, upper = _bounds(s)
    ret = Highs_addRow(
        model,
        lower,
        upper,
        length(indices),
        indices,
        coefficients,
    )
    _check_ret(ret)
    return MOI.ConstraintIndex{typeof(f),typeof(s)}(key.value)
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:_SCALAR_SETS},
)
    r = row(model, c)
    ret = Highs_deleteRowsByRange(model, r, r)
    _check_ret(ret)
    for info in values(model.affine_constraint_info)
        if info.row > r
            info.row -= 1
        end
    end
    key = _ConstraintKey(c.value)
    delete!(model.affine_constraint_info, key)
    model.name_to_constraint_index = nothing
    return
end

function MOI.delete(
    model::Optimizer,
    c::Vector{MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S}},
) where {S<:_SCALAR_SETS}
    rows = sort(row.(model, c))
    ret = Highs_deleteRowsBySet(model, length(rows), rows)
    _check_ret(ret)
    for info in values(model.affine_constraint_info)
        i = findlast(Base.Fix2(<, info.row), rows)
        if i !== nothing
            info.row -= i
        end
    end
    for ci in c
        delete!(model.affine_constraint_info, _ConstraintKey(ci.value))
    end
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
) where {S<:_SCALAR_SETS}
    return _set(_info(model, c))
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
    s::S,
) where {S<:_SCALAR_SETS}
    lower, upper = _bounds(s)
    info = _info(model, c)
    ret = Highs_changeRowBounds(model, info.row, lower, upper)
    _check_ret(ret)
    info.lower, info.upper = lower, upper
    return
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::Vector{MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S}},
    s::Vector{S},
) where {S<:_SCALAR_SETS}
    if length(c) != length(s)
        msg = "number of constraints does not match number of sets"
        throw(DimensionMismatch(msg))
    end
    for ci in c
        MOI.throw_if_not_valid(model, ci)
    end
    N = length(c)
    rows, lower, upper = zeros(Cint, N), zeros(Cdouble, N), zeros(Cdouble, N)
    for (i, (ci, si)) in enumerate(zip(c, s))
        info = _info(model, ci)
        rows[i] = info.row
        lower[i], upper[i] = info.lower, info.upper = _bounds(si)
    end
    ret = Highs_changeRowsBoundsBySet(model, N, rows, lower, upper)
    _check_ret(ret)
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
) where {S<:_SCALAR_SETS}
    r = row(model, c)
    num_row = Ref{HighsInt}(0)
    num_nz = Ref{HighsInt}(0)
    matrix_start = Ref{HighsInt}(0)
    lower, upper = Ref{Cdouble}(), Ref{Cdouble}()
    _ = Highs_getRowsByRange(
        model,
        r,
        r,
        num_row,
        lower,
        upper,
        num_nz,
        matrix_start,
        C_NULL,
        C_NULL,
    )
    matrix_index = Vector{HighsInt}(undef, num_nz[])
    matrix_value = Vector{Cdouble}(undef, num_nz[])
    ret = Highs_getRowsByRange(
        model,
        r,
        r,
        num_row,
        lower,
        upper,
        num_nz,
        matrix_start,
        matrix_index,
        matrix_value,
    )
    _check_ret(ret)
    return MOI.ScalarAffineFunction(
        MOI.ScalarAffineTerm{Float64}[
            MOI.ScalarAffineTerm(
                val,
                model.variable_info[CleverDicts.LinearIndex(col + 1)].index,
            ) for (col, val) in zip(matrix_index, matrix_value) if !iszero(val)
        ],
        0.0,
    )
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ConstraintName,
    ::Type{
        <:MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:_SCALAR_SETS},
    },
)
    return true
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:Any},
)
    return _info(model, c).name
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:Any},
    name::String,
)
    info = _info(model, c)
    if info.name == name
        return  # The name is unchanged
    end
    info.name = name
    model.name_to_constraint_index = nothing
    if isempty(name)
        # HiGHS does not support passing "" or C_NULL for names. So we default
        # to a fairly standard R$i for the crowolumn name.
        name = string("R", info.row)
    end
    ret = Highs_passRowName(model, info.row, name)
    _check_ret(ret)
    return
end

function MOI.get(model::Optimizer, ::Type{MOI.ConstraintIndex}, name::String)
    if model.name_to_constraint_index === nothing
        _rebuild_name_to_constraint_index(model)
    end
    constraint = get(model.name_to_constraint_index, name, missing)
    if constraint === nothing
        error("Duplicate constraint name detected: $(name)")
    elseif constraint === missing
        return nothing  # No constraint exists with this name
    end
    return constraint
end

function MOI.get(
    model::Optimizer,
    C::Type{MOI.ConstraintIndex{F,S}},
    name::String,
) where {F,S}
    index = MOI.get(model, MOI.ConstraintIndex, name)
    if index isa C
        return index::MOI.ConstraintIndex{F,S}
    end
    return nothing
end

function _rebuild_name_to_constraint_index(model::Optimizer)
    model.name_to_constraint_index =
        Dict{String,Union{Nothing,MOI.ConstraintIndex}}()
    for (key, info) in model.affine_constraint_info
        if isempty(info.name)
            continue
        end
        S = typeof(_set(info))
        _set_name_to_constraint_index(
            model,
            info.name,
            MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S}(key.value),
        )
    end
    return
end

function _set_name_to_constraint_index(
    model::Optimizer,
    name::String,
    index::MOI.ConstraintIndex,
)
    if haskey(model.name_to_constraint_index, name)
        model.name_to_constraint_index[name] = nothing
    else
        model.name_to_constraint_index[name] = index
    end
    return
end

function MOI.modify(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
    chg::MOI.ScalarCoefficientChange{Float64},
) where {S<:_SCALAR_SETS}
    ret = Highs_changeCoeff(
        model,
        row(model, c),
        column(model, chg.variable),
        chg.new_coefficient,
    )
    _check_ret(ret)
    return
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{F,S},
    f::F,
) where {F<:MOI.ScalarAffineFunction{Float64},S<:_SCALAR_SETS}
    if !iszero(f.constant)
        throw(MOI.ScalarFunctionConstantNotZero{Float64,F,S}(f.constant))
    end
    # This loop is a bit slow because we query the existing function first, and
    # then set each coefficient one-by-one. It's probably okay until someone
    # complains but it will require some upstream changes to introduce a faster
    # way of modifying a list of coefficients.
    old = MOI.get(model, MOI.ConstraintFunction(), c)
    terms = Dict(x.variable => 0.0 for x in old.terms)
    for term in f.terms
        terms[term.variable] = term.coefficient
    end
    r = row(model, c)
    for (k, v) in terms
        ret = Highs_changeCoeff(model, r, column(model, k), v)
        _check_ret(ret)
    end
    return
end

###
### Optimize methods.
###

"""
    _store_solution(model::Optimizer, ret::HighsInt)

Get the solution from a run of HiGHS.

`ret` is the `kHighsStatus` constant from `Highs_run`.
"""
function _store_solution(model::Optimizer, ret::HighsInt)
    x = model.solution
    x.status = ret == kHighsStatusError ? _OPTIMIZE_ERRORED : _OPTIMIZE_OK
    x.primal_solution_status = kHighsSolutionStatusNone
    x.dual_solution_status = kHighsSolutionStatusNone
    x.has_dual_ray = false
    x.has_primal_ray = false
    numCols = Highs_getNumCols(model)
    numRows = Highs_getNumRows(model)
    resize!(x.colvalue, numCols)
    resize!(x.coldual, numCols)
    resize!(x.rowvalue, numRows)
    resize!(x.rowdual, numRows)
    x.model_status = Highs_getModelStatus(model)
    statusP = Ref{HighsInt}()
    Highs_getIntInfoValue(model, "primal_solution_status", statusP)
    x.primal_solution_status = statusP[]
    Highs_getIntInfoValue(model, "dual_solution_status", statusP)
    x.dual_solution_status = statusP[]
    if x.primal_solution_status != kHighsSolutionStatusNone
        Highs_getSolution(model, x.colvalue, x.coldual, x.rowvalue, x.rowdual)
        if model.hessian === nothing
            # No basis is present in a QP.
            resize!(x.colstatus, numCols)
            resize!(x.rowstatus, numRows)
            Highs_getBasis(model, x.colstatus, x.rowstatus)
        end
    end
    if !MOI.get(model, ComputeInfeasibilityCertificate())
        return
    end
    if x.model_status == kHighsModelStatusInfeasible
        ret = Highs_getDualRay(model, statusP, x.rowdual)
        # Don't `_check_ret(ret)` here, just bail is there isn't a dual ray.
        x.has_dual_ray = (ret == kHighsStatusOk) && (statusP[] == 1)
        if x.has_dual_ray
            _compute_farkas_variable_dual(model, x.coldual)
            x.dual_solution_status = kHighsSolutionStatusNone
        end
    elseif x.model_status == kHighsModelStatusUnbounded
        ret = Highs_getPrimalRay(model, statusP, x.colvalue)
        # Don't `_check_ret(ret)` here, just bail is there isn't a dual ray.
        x.has_primal_ray = (ret == kHighsStatusOk) && (statusP[] == 1)
        if x.has_primal_ray
            x.primal_solution_status = kHighsSolutionStatusNone
        end
    end
    return
end

"""
    _compute_farkas_variable_dual(model::Optimizer, dual::Vector{Cdouble})

Return a Farkas dual associated with the variable bounds of `col`. Given a dual
ray:

    ā * x = λ' * A * x <= λ' * b = -β + sum(āᵢ * Uᵢ | āᵢ < 0) + sum(āᵢ * Lᵢ | āᵢ > 0)

The Farkas dual of the variable is ā, and it applies to the upper bound if ā < 0,
and it applies to the lower bound if ā > 0.
"""
function _compute_farkas_variable_dual(model::Optimizer, dual::Vector{Cdouble})
    dual_objective_value = 0.0
    # Column components of the dual objective value
    num_nz, num_cols = Ref{HighsInt}(0), Ref{HighsInt}(0)
    m = Highs_getNumRows(model)
    n = Highs_getNumCols(model)
    ret = Highs_getColsByRange(
        model,
        0,
        n - 1,
        num_cols,
        C_NULL,     # c_obj
        C_NULL,     # lower
        C_NULL,     # upper
        num_nz,
        C_NULL,     # matrix_start,
        C_NULL,     # matrix_index,
        C_NULL,     # matrix_value,
    )
    _check_ret(ret)
    matrix_start = zeros(HighsInt, n)
    matrix_index = zeros(HighsInt, num_nz[])
    matrix_value = zeros(Cdouble, num_nz[])
    lower, upper = zeros(n), zeros(n)
    push!(matrix_start, num_nz[])
    ret = Highs_getColsByRange(
        model,
        0,
        n - 1,
        num_cols,
        C_NULL,
        lower,
        upper,
        num_nz,
        matrix_start,
        matrix_index,
        matrix_value,
    )
    _check_ret(ret)
    for (i, (li, ui)) in enumerate(zip(lower, upper))
        dual[i] = 0.0
        for k in (matrix_start[i]+1):matrix_start[i+1]
            row = matrix_index[k]
            dual[i] -= model.solution.rowdual[row+1] * matrix_value[k]
        end
    end
    return
end

function _set_variable_primal_start(model::Optimizer)
    if all(info -> info.start === nothing, values(model.variable_info))
        return
    end
    start = zeros(Cdouble, length(model.variable_info))
    for (x, info) in model.variable_info
        # For the default start, pick the lower bound if it exists, otherwise
        # the minimum of the upper bound and zero.
        default = isfinite(info.lower) ? info.lower : min(info.upper, 0.0)
        start[info.column+1] = something(info.start, default)
    end
    ret = Highs_setSolution(model, start, C_NULL, C_NULL, C_NULL)
    _check_ret(ret)
    return
end

function MOI.optimize!(model::Optimizer)
    for info in model.binaries
        Highs_changeColBounds(
            model,
            info.column,
            max(0.0, info.lower),
            min(1.0, info.upper),
        )
    end
    _set_variable_primal_start(model)
    if Highs_versionMajor() == 1 && Highs_versionMinor() >= 6
        # Not available in earlier versions
        ret = Highs_zeroAllClocks(model)
        _check_ret(ret)
    end
    has_default_callback = model.callback_data === nothing
    if has_default_callback
        # From the docstring of disable_sigint, "External functions that do not
        # call julia code or julia runtime automatically disable sigint during
        # their execution." We don't want this though! We want to be able to
        # SIGINT HiGHS, and then catch it as an interrupt. As a hack, until
        # Julia introduces an interruptible ccall --- which it likely won't
        # https://github.com/JuliaLang/julia/issues/2622 --- set a null
        # callback.
        cb_types = [
            # See Issue 291. The SimplexInterrupt callback in HiGHS@1.11 is very
            # slow.
            # kHighsCallbackSimplexInterrupt,
            kHighsCallbackIpmInterrupt,
            kHighsCallbackMipInterrupt,
        ]
        MOI.set(model, CallbackFunction(cb_types), (args...) -> Cint(0))
    end
    # if `Highs_run` implicitly uses memory or other resources owned by `model`, preserve it
    GC.@preserve model begin
        # allow Julia to GC while this thread is in `Highs_run`
        gc_state = ccall(:jl_gc_safe_enter, Int8, ())
        # We disable sigint here so that it can be called only when we are in a
        # try-catch of our CallbackFunction.
        ret = disable_sigint(() -> Highs_run(model))
        # leave GC-safe region, waiting for GC to complete if it's running
        ccall(:jl_gc_safe_leave, Cvoid, (Int8,), gc_state)
    end
    if has_default_callback
        MOI.set(model, CallbackFunction(), nothing)
    end
    _store_solution(model, ret)
    # TODO(odow): resetting the bounds here invalidates previously stored
    #             solutions.
    # for info in model.binaries
    #     Highs_changeColBounds(model, info.column, info.lower, info.upper)
    # end
    return
end

const _TerminationStatusMap = Dict(
    kHighsModelStatusNotset => (MOI.OTHER_ERROR, "kHighsModelStatusNotset"),
    kHighsModelStatusLoadError =>
        (MOI.OTHER_ERROR, "kHighsModelStatusLoadError"),
    kHighsModelStatusModelError =>
        (MOI.INVALID_MODEL, "kHighsModelStatusModelError"),
    kHighsModelStatusPresolveError =>
        (MOI.OTHER_ERROR, "kHighsModelStatusPresolveError"),
    kHighsModelStatusSolveError =>
        (MOI.OTHER_ERROR, "kHighsModelStatusSolveError"),
    kHighsModelStatusPostsolveError =>
        (MOI.OTHER_ERROR, "kHighsModelStatusPostsolveError"),
    kHighsModelStatusModelEmpty =>
        (MOI.INVALID_MODEL, "kHighsModelStatusModelEmpty"),
    kHighsModelStatusOptimal => (MOI.OPTIMAL, "kHighsModelStatusOptimal"),
    kHighsModelStatusInfeasible =>
        (MOI.INFEASIBLE, "kHighsModelStatusInfeasible"),
    kHighsModelStatusUnboundedOrInfeasible => (
        MOI.INFEASIBLE_OR_UNBOUNDED,
        "kHighsModelStatusUnboundedOrInfeasible",
    ),
    kHighsModelStatusUnbounded =>
        (MOI.DUAL_INFEASIBLE, "kHighsModelStatusUnbounded"),
    kHighsModelStatusObjectiveBound =>
        (MOI.OBJECTIVE_LIMIT, "kHighsModelStatusObjectiveBound"),
    kHighsModelStatusObjectiveTarget =>
        (MOI.OBJECTIVE_LIMIT, "kHighsModelStatusObjectiveTarget"),
    kHighsModelStatusTimeLimit =>
        (MOI.TIME_LIMIT, "kHighsModelStatusTimeLimit"),
    kHighsModelStatusIterationLimit =>
        (MOI.ITERATION_LIMIT, "kHighsModelStatusIterationLimit"),
    kHighsModelStatusUnknown =>
        (MOI.OTHER_ERROR, "kHighsModelStatusUnknown"),
    kHighsModelStatusSolutionLimit =>
        (MOI.SOLUTION_LIMIT, "kHighsModelStatusSolutionLimit"),
    kHighsModelStatusInterrupt =>
        (MOI.INTERRUPTED, "kHighsModelStatusInterrupt"),
)

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if model.solution.status == _OPTIMIZE_NOT_CALLED
        return MOI.OPTIMIZE_NOT_CALLED
    elseif model.solution.status == _OPTIMIZE_ERRORED
        return MOI.OTHER_ERROR
    end
    return _TerminationStatusMap[model.solution.model_status][1]
end

function MOI.get(model::Optimizer, ::MOI.ResultCount)
    x = model.solution
    if x.primal_solution_status != kHighsSolutionStatusNone
        return 1
    elseif x.has_dual_ray || x.has_primal_ray
        return 1
    else
        return 0
    end
end

function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    if model.solution.status == _OPTIMIZE_NOT_CALLED
        return "OPTIMIZE_NOT_CALLED"
    elseif model.solution.status == _OPTIMIZE_ERRORED
        return "There was an error calling optimize!"
    end
    return _TerminationStatusMap[model.solution.model_status][2]
end

function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    sol = model.solution
    if attr.result_index != 1
        return MOI.NO_SOLUTION
    elseif sol.primal_solution_status == kHighsSolutionStatusFeasible
        return MOI.FEASIBLE_POINT
    elseif sol.primal_solution_status == kHighsSolutionStatusInfeasible
        return MOI.INFEASIBLE_POINT
    elseif sol.has_primal_ray
        return MOI.INFEASIBILITY_CERTIFICATE
    end
    return MOI.NO_SOLUTION
end

function MOI.get(model::Optimizer, attr::MOI.DualStatus)
    sol = model.solution
    if attr.result_index != 1
        return MOI.NO_SOLUTION
    elseif sol.dual_solution_status == kHighsSolutionStatusFeasible
        return MOI.FEASIBLE_POINT
    elseif sol.dual_solution_status == kHighsSolutionStatusInfeasible
        return MOI.INFEASIBLE_POINT
    elseif model.solution.has_dual_ray
        return MOI.INFEASIBILITY_CERTIFICATE
    end
    return MOI.NO_SOLUTION
end

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    if model.solution.has_primal_ray
        return MOI.Utilities.get_fallback(model, attr)
    end
    return Highs_getObjectiveValue(model)
end

function MOI.get(model::Optimizer, attr::MOI.DualObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    dual_objective_value = 0.0
    if model.solution.has_dual_ray
        # Do nothing
    elseif model.solution.dual_solution_status == kHighsSolutionStatusNone
        # For MIPs, we cannot compute a dual objective value
        return NaN
    else
        offset = Ref{Cdouble}()
        ret = Highs_getObjectiveOffset(model, offset)
        _check_ret(ret)
        dual_objective_value += offset[]
    end
    # Column components of the dual objective value
    set = Cint[
        i - 1 for (i, d) in enumerate(model.solution.coldual) if !iszero(d)
    ]
    n = length(set)
    lower, upper = zeros(n), zeros(n)
    ret = Highs_getColsBySet(
        model,
        n,
        set,
        Ref{Cint}(),  # num_col
        C_NULL,
        lower,
        upper,
        Ref{Cint}(),  # num_nnz
        C_NULL,  # matrix_start,
        C_NULL,  # matrix_index,
        C_NULL,  # matrix_value,
    )
    _check_ret(ret)
    sense = _sense_corrector(model)
    for (li, i, ui) in zip(lower, set, upper)
        di = model.solution.coldual[i+1]
        if model.solution.has_dual_ray
            dual_objective_value += sense * _active_bound(li, NaN, ui, di) * di
        else
            xi = model.solution.colvalue[i+1]
            dual_objective_value += _active_bound(li, xi, ui, di) * di
        end
    end
    # Row components of the dual objective value
    set = Cint[
        i - 1 for (i, d) in enumerate(model.solution.rowdual) if !iszero(d)
    ]
    n = length(set)
    resize!(lower, n)
    resize!(upper, n)
    ret = Highs_getRowsBySet(
        model,
        n,
        set,
        Ref{Cint}(),  # num_row
        lower,
        upper,
        Ref{Cint}(),  # num_nnz
        C_NULL,  # matrix_start,
        C_NULL,  # matrix_index,
        C_NULL,  # matrix_value,
    )
    _check_ret(ret)
    for (li, i, ui) in zip(lower, set, upper)
        di = model.solution.rowdual[i+1]
        if model.solution.has_dual_ray
            dual_objective_value += sense * _active_bound(li, NaN, ui, di) * di
        else
            ri = model.solution.rowvalue[i+1]
            dual_objective_value += _active_bound(li, ri, ui, di) * di
        end
    end
    return dual_objective_value
end

function _active_bound(l, x, u, d)
    if isfinite(l) && isfinite(u)
        if isfinite(x)
            # Pick the bound that is closest to the primal value
            return ifelse(abs(x - l) < abs(x - u), l, u)
        else
            # Pick the bound depending on the sign of the dual value
            return ifelse(d >= 0, l, u)
        end
    elseif isfinite(l)
        return l
    elseif isfinite(u)
        return u
    else
        return 0.0
    end
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveBound)
    if model.solution.dual_solution_status != kHighsSolutionStatusNone
        return MOI.get(model, MOI.DualObjectiveValue())
    end
    p = Ref{Cdouble}()
    Highs_getDoubleInfoValue(model, "mip_dual_bound", p)
    return p[]
end

function MOI.get(model::Optimizer, ::MOI.SolveTimeSec)
    return Highs_getRunTime(model)
end

function MOI.get(model::Optimizer, ::MOI.RelativeGap)
    if model.solution.dual_solution_status != kHighsSolutionStatusNone
        primal = MOI.get(model, MOI.ObjectiveValue())
        dual = MOI.get(model, MOI.DualObjectiveValue())
        return abs(primal - dual) / max(1, abs(primal))
    end
    p = Ref{Cdouble}(0.0)
    ret = Highs_getDoubleInfoValue(model, "mip_gap", p)
    _check_ret(ret)
    return p[]
end

function MOI.get(model::Optimizer, ::MOI.SimplexIterations)
    p = Ref{HighsInt}(0)
    ret = Highs_getIntInfoValue(model, "simplex_iteration_count", p)
    _check_ret(ret)
    return Int64(p[])
end

function MOI.get(model::Optimizer, ::MOI.BarrierIterations)
    p = Ref{HighsInt}(0)
    ret = Highs_getIntInfoValue(model, "ipm_iteration_count", p)
    _check_ret(ret)
    return Int64(p[])
end

function MOI.get(model::Optimizer, ::MOI.NodeCount)
    p = Ref{Int64}(0)
    ret = Highs_getInt64InfoValue(model, "mip_node_count", p)
    _check_ret(ret)
    return p[]
end

function MOI.get(
    model::Optimizer,
    attr::MOI.VariablePrimal,
    x::MOI.VariableIndex,
)
    MOI.check_result_index_bounds(model, attr)
    return model.solution.colvalue[column(model, x)+1]
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.VariableIndex,<:_SCALAR_SETS},
)
    MOI.check_result_index_bounds(model, attr)
    return MOI.get(model, MOI.VariablePrimal(), MOI.VariableIndex(c.value))
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:_SCALAR_SETS},
)
    MOI.check_result_index_bounds(model, attr)
    if model.solution.has_primal_ray
        return MOI.Utilities.get_fallback(model, attr, c)
    end
    return model.solution.rowvalue[row(model, c)+1]
end

"""
    _sense_corrector(model::Optimizer)

Return `-1` if `MAX_SENSE`, and `1` is `MIN_SENSE`. Useful for correcting
solution attributes which are sense-dependent.
"""
function _sense_corrector(model::Optimizer)
    senseP = Ref{HighsInt}()
    ret = Highs_getObjectiveSense(model, senseP)
    _check_ret(ret)
    return senseP[]
end

"""
    _signed_dual(dual::Float64, ::Type{Set})

A heuristic for determining whether the dual of an interval constraint applies
to the lower or upper bound. It can be wrong by at most the solver's tolerance.
"""
_signed_dual(dual::Float64, ::Type{MOI.LessThan{Float64}}) = min(dual, 0.0)
_signed_dual(dual::Float64, ::Type{MOI.GreaterThan{Float64}}) = max(dual, 0.0)
_signed_dual(dual::Float64, ::Any) = dual

"""
    _signed_dual(dual::Float64, ::Type{Set}, status::HighsInt)

Determine whether the dual of an interval constraint applies to the lower or
upper bound using the basis status reported by HiGHS.
"""
function _signed_dual(
    dual::Float64,
    ::Type{MOI.LessThan{Float64}},
    status::HighsInt,
)
    return status == kHighsBasisStatusUpper ? dual : 0.0
end

function _signed_dual(
    dual::Float64,
    ::Type{MOI.GreaterThan{Float64}},
    status::HighsInt,
)
    return status == kHighsBasisStatusLower ? dual : 0.0
end

function _signed_dual(dual::Float64, ::Type{<:MOI.LessThan}, ::Nothing)
    return min(dual, 0.0)
end

function _signed_dual(dual::Float64, ::Type{<:MOI.GreaterThan}, ::Nothing)
    return max(dual, 0.0)
end

_signed_dual(dual::Float64, ::Any, ::Any) = dual

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.VariableIndex,S},
) where {S<:_SCALAR_SETS}
    MOI.check_result_index_bounds(model, attr)
    col = column(model, c)
    if model.solution.has_dual_ray[] == 1
        return _signed_dual(model.solution.coldual[col+1], S)
    end
    dual = _sense_corrector(model) * model.solution.coldual[col+1]
    stat = get(model.solution.colstatus, col + 1, nothing)
    return _signed_dual(dual, S, stat)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
) where {S<:_SCALAR_SETS}
    MOI.check_result_index_bounds(model, attr)
    r = row(model, c) + 1
    dual = model.solution.rowdual[r]
    if model.solution.has_dual_ray[] == 1
        return _signed_dual(dual, S)
    end
    stat = get(model.solution.rowstatus, r, nothing)
    return _signed_dual(_sense_corrector(model) * dual, S, stat)
end

###
### MOI.ConstraintBasisStatus
###

const _BASIS_STATUS_MAP = Dict(
    kHighsBasisStatusLower => MOI.NONBASIC_AT_LOWER,
    kHighsBasisStatusBasic => MOI.BASIC,
    kHighsBasisStatusUpper => MOI.NONBASIC_AT_UPPER,
    kHighsBasisStatusZero => MOI.SUPER_BASIC,
    kHighsBasisStatusNonbasic => MOI.NONBASIC,
)

_replace_basis_status(b::MOI.BasisStatusCode, ::Type{<:MOI.Interval}) = b

function _replace_basis_status(basis::MOI.BasisStatusCode, ::Type{S}) where {S}
    if basis == MOI.NONBASIC_AT_LOWER || basis == MOI.NONBASIC_AT_UPPER
        return MOI.NONBASIC
    end
    return basis
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintBasisStatus,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
) where {S<:_SCALAR_SETS}
    MOI.check_result_index_bounds(model, attr)
    stat = get(model.solution.rowstatus, row(model, c) + 1, nothing)
    if stat === nothing
        throw(MOI.GetAttributeNotAllowed(attr, "no basis is present"))
    end
    return _replace_basis_status(_BASIS_STATUS_MAP[stat], S)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.VariableBasisStatus,
    x::MOI.VariableIndex,
)
    MOI.check_result_index_bounds(model, attr)
    stat = get(model.solution.colstatus, column(model, x) + 1, nothing)
    if stat === nothing
        throw(MOI.GetAttributeNotAllowed(attr, "no basis is present"))
    end
    return _BASIS_STATUS_MAP[stat]
end

###
### Integrality
###

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VariableIndex},
    ::Type{<:Union{MOI.ZeroOne,MOI.Integer}},
)
    return true
end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Integer},
)
    info = get(model.variable_info, MOI.VariableIndex(ci.value), nothing)
    if info === nothing
        return false
    end
    return info.type == _TYPE_INTEGER
end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.ZeroOne},
)
    info = get(model.variable_info, MOI.VariableIndex(ci.value), nothing)
    if info === nothing
        return false
    end
    return info.type == _TYPE_BINARY
end

function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.Integer},
)
    indices = MOI.ConstraintIndex{MOI.VariableIndex,MOI.Integer}[
        MOI.ConstraintIndex{MOI.VariableIndex,MOI.Integer}(key.value) for
        (key, info) in model.variable_info if info.type == _TYPE_INTEGER
    ]
    return sort!(indices, by = x -> x.value)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.ZeroOne},
)
    indices = MOI.ConstraintIndex{MOI.VariableIndex,MOI.ZeroOne}[
        MOI.ConstraintIndex{MOI.VariableIndex,MOI.ZeroOne}(key.value) for
        (key, info) in model.variable_info if info.type == _TYPE_BINARY
    ]
    return sort!(indices, by = x -> x.value)
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.VariableIndex,
    ::MOI.Integer,
)
    info = _info(model, f)
    info.type = _TYPE_INTEGER
    ci = MOI.ConstraintIndex{MOI.VariableIndex,MOI.Integer}(f.value)
    Highs_changeColIntegrality(model, info.column, kHighsVarTypeInteger)
    return ci
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.VariableIndex,
    ::MOI.ZeroOne,
)
    info = _info(model, f)
    info.type = _TYPE_BINARY
    ci = MOI.ConstraintIndex{MOI.VariableIndex,MOI.ZeroOne}(f.value)
    Highs_changeColIntegrality(model, info.column, kHighsVarTypeInteger)
    push!(model.binaries, info)
    return ci
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Integer},
)
    MOI.throw_if_not_valid(model, ci)
    return MOI.Integer()
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.ZeroOne},
)
    MOI.throw_if_not_valid(model, ci)
    return MOI.ZeroOne()
end

function MOI.delete(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.Integer},
)
    MOI.throw_if_not_valid(model, ci)
    info = _info(model, ci)
    info.type = _TYPE_CONTINUOUS
    Highs_changeColIntegrality(model, info.column, kHighsVarTypeContinuous)
    return
end

function MOI.delete(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,MOI.ZeroOne},
)
    MOI.throw_if_not_valid(model, ci)
    info = _info(model, ci)
    info.type = _TYPE_CONTINUOUS
    Highs_changeColIntegrality(model, info.column, kHighsVarTypeContinuous)
    # optimize! sets bounds to [0, 1] if not present, so if we delete the
    # ZeroOne constraint after a call to optimize!, then we need to reset the
    # bounds.
    Highs_changeColBounds(model, info.column, info.lower, info.upper)
    delete!(model.binaries, info)
    return
end

###
### Semicontinuous and Semiinteger
###

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VariableIndex},
    ::Type{S},
) where {S<:Union{MOI.Semicontinuous{Float64},MOI.Semiinteger}}
    return true
end

_bound_type(::Type{<:MOI.Semicontinuous}) = _BOUND_SEMI_CONTINUOUS
_bound_type(::Type{<:MOI.Semiinteger}) = _BOUND_SEMI_INTEGER

_type_type(::Type{<:MOI.Semicontinuous}) = _TYPE_SEMI_CONTINUOUS
_type_type(::Type{<:MOI.Semiinteger}) = _TYPE_SEMI_INTEGER

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,S},
) where {S<:Union{MOI.Semicontinuous{Float64},MOI.Semiinteger{Float64}}}
    info = get(model.variable_info, MOI.VariableIndex(c.value), nothing)
    if info === nothing
        return false
    end
    return info.bound == _bound_type(S)
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.VariableIndex,
    s::S,
) where {S<:Union{MOI.Semicontinuous{Float64},MOI.Semiinteger{Float64}}}
    info = _info(model, f)
    _throw_if_existing_lower(info, s)
    _throw_if_existing_upper(info, s)
    info.bound, info.type = _bound_type(S), _type_type(S)
    info.lower, info.upper = s.lower, s.upper
    index = MOI.ConstraintIndex{MOI.VariableIndex,S}(f.value)
    MOI.set(model, MOI.ConstraintSet(), index, s)
    return index
end

function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintIndices{MOI.VariableIndex,S},
) where {S<:Union{MOI.Semicontinuous{Float64},MOI.Semiinteger{Float64}}}
    indices = MOI.ConstraintIndex{MOI.VariableIndex,S}[
        MOI.ConstraintIndex{MOI.VariableIndex,S}(key.value) for
        (key, info) in model.variable_info if info.type == _type_type(S)
    ]
    return sort!(indices, by = x -> x.value)
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.VariableIndex,S},
) where {S<:Union{MOI.Semicontinuous{Float64},MOI.Semiinteger{Float64}}}
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    ret = Highs_changeColBounds(model, info.column, -Inf, Inf)
    _check_ret(ret)
    info.lower, info.upper = -Inf, Inf
    info.bound, info.type = _BOUND_NONE, _TYPE_CONTINUOUS
    ret =
        Highs_changeColIntegrality(model, info.column, kHighsVarTypeContinuous)
    _check_ret(ret)
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex,S},
) where {S<:Union{MOI.Semicontinuous{Float64},MOI.Semiinteger{Float64}}}
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    return S(info.lower, info.upper)
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex,S},
    set::S,
) where {S<:Union{MOI.Semicontinuous{Float64},MOI.Semiinteger{Float64}}}
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    ret = Highs_changeColBounds(model, info.column, set.lower, set.upper)
    _check_ret(ret)
    var_type = if S == MOI.Semicontinuous{Float64}
        kHighsVarTypeSemiContinuous
    else
        kHighsVarTypeSemiInteger
    end
    ret = Highs_changeColIntegrality(model, info.column, var_type)
    _check_ret(ret)
    return
end

###
### MOI.copy_to
###

function _add_bounds(::Vector{Float64}, ub, i, s::MOI.LessThan{Float64})
    ub[i] = s.upper
    return
end

function _add_bounds(lb, ::Vector{Float64}, i, s::MOI.GreaterThan{Float64})
    lb[i] = s.lower
    return
end

function _add_bounds(lb, ub, i, s::MOI.EqualTo{Float64})
    lb[i], ub[i] = s.value, s.value
    return
end

function _add_bounds(lb, ub, i, s::MOI.Interval{Float64})
    lb[i], ub[i] = s.lower, s.upper
    return
end

_constraints(src, F, S) = MOI.get(src, MOI.ListOfConstraintIndices{F,S}())

function _extract_bound_data(
    dest::Optimizer,
    src::MOI.ModelLike,
    mapping,
    collower::Vector{Float64},
    colupper::Vector{Float64},
    ::Type{S},
) where {S}
    for c_index in _constraints(src, MOI.VariableIndex, S)
        f = MOI.get(src, MOI.ConstraintFunction(), c_index)
        s = MOI.get(src, MOI.ConstraintSet(), c_index)
        new_f = mapping[f]
        info = _info(dest, new_f)
        _add_bounds(collower, colupper, info.column + 1, s)
        _update_info(info, s)
        mapping[c_index] = MOI.ConstraintIndex{MOI.VariableIndex,S}(new_f.value)
    end
    return
end

function _copy_to_columns(dest::Optimizer, src::MOI.ModelLike, mapping)
    x_src = MOI.get(src, MOI.ListOfVariableIndices())
    numcols = HighsInt(length(x_src))
    for (i, x) in enumerate(x_src)
        index = CleverDicts.add_item(
            dest.variable_info,
            _VariableInfo(MOI.VariableIndex(0), HighsInt(0)),
        )
        info = _info(dest, index)
        if MOI.supports(src, MOI.VariableName(), MOI.VariableIndex)
            info.name = MOI.get(src, MOI.VariableName(), x)
        end
        info.index = index
        info.column = HighsInt(i - 1)
        mapping[x] = index
    end
    return numcols
end

_add_sizehint!(vec, n) = sizehint!(vec, length(vec) + n)

function _extract_row_data(
    dest::Optimizer,
    src::MOI.ModelLike,
    mapping,
    rowlower::Vector{Float64},
    rowupper::Vector{Float64},
    I::Vector{HighsInt},
    J::Vector{HighsInt},
    V::Vector{Float64},
    ::Type{S},
) where {S}
    F = MOI.ScalarAffineFunction{Float64}
    row = length(rowlower) + 1
    list = _constraints(src, MOI.ScalarAffineFunction{Float64}, S)
    numrows = length(list)
    _add_sizehint!(rowlower, numrows)
    _add_sizehint!(rowupper, numrows)
    n_terms = 0
    fs = Array{F}(undef, numrows)
    for (i, c_index) in enumerate(list)
        f = MOI.get(src, MOI.ConstraintFunction(), c_index)
        fs[i] = f
        set = MOI.get(src, MOI.ConstraintSet(), c_index)
        l, u = _bounds(set)
        push!(rowlower, l - f.constant)
        push!(rowupper, u - f.constant)
        n_terms += length(f.terms)
        key = CleverDicts.add_item(
            dest.affine_constraint_info,
            _ConstraintInfo(set),
        )
        info = dest.affine_constraint_info[key]
        info.row = HighsInt(length(dest.affine_constraint_info) - 1)
        if MOI.supports(src, MOI.ConstraintName(), MOI.ConstraintIndex{F,S})
            info.name = MOI.get(src, MOI.ConstraintName(), c_index)
        end
        mapping[c_index] = MOI.ConstraintIndex{F,S}(key.value)
    end
    _add_sizehint!(I, n_terms)
    _add_sizehint!(J, n_terms)
    _add_sizehint!(V, n_terms)
    for f in fs
        for term in f.terms
            push!(I, row)
            push!(J, HighsInt(mapping[term.variable].value))
            push!(V, term.coefficient)
        end
        row += 1
    end
    return
end

function _check_input_data(dest::Optimizer, src::MOI.ModelLike)
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        if !MOI.supports_constraint(dest, F, S)
            throw(
                MOI.UnsupportedConstraint{F,S}(
                    "HiGHS does not support constraints of type $F-in-$S.",
                ),
            )
        end
        for attr in MOI.get(src, MOI.ListOfConstraintAttributesSet{F,S}())
            if attr in (
                MOI.ConstraintName(),
                MOI.ConstraintPrimalStart(),
                MOI.ConstraintDualStart(),
            )
                continue
            else
                throw(MOI.UnsupportedAttribute(attr))
            end
        end
    end
    for attr in MOI.get(src, MOI.ListOfModelAttributesSet())
        if attr in (
            MOI.Name(),
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
            MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(),
            MOI.ObjectiveSense(),
        )
            continue
        else
            throw(MOI.UnsupportedAttribute(attr))
        end
    end
    for attr in MOI.get(src, MOI.ListOfVariableAttributesSet())
        if attr in (MOI.VariableName(), MOI.VariablePrimalStart())
            continue
        else
            throw(MOI.UnsupportedAttribute(attr))
        end
    end
    return
end

MOI.supports_incremental_interface(::Optimizer) = true

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    MOI.empty!(dest)
    _check_input_data(dest, src)
    mapping = MOI.Utilities.IndexMap()
    numcol = _copy_to_columns(dest, src, mapping)
    collower, colupper = fill(-Inf, numcol), fill(Inf, numcol)
    rowlower, rowupper = Float64[], Float64[]
    I, J, V = HighsInt[], HighsInt[], Float64[]
    for S in (
        MOI.GreaterThan{Float64},
        MOI.LessThan{Float64},
        MOI.EqualTo{Float64},
        MOI.Interval{Float64},
    )
        _extract_bound_data(dest, src, mapping, collower, colupper, S)
        _extract_row_data(dest, src, mapping, rowlower, rowupper, I, J, V, S)
    end
    numrow = HighsInt(length(rowlower))
    A = SparseArrays.sparse(I, J, V, numrow, numcol)
    # Extract integrality constraints
    integrality = fill(kHighsVarTypeContinuous, numcol)
    has_integrality = false
    for ci in _constraints(src, MOI.VariableIndex, MOI.ZeroOne)
        info = _info(dest, ci)
        info.type = _TYPE_BINARY
        integrality[info.column+1] = kHighsVarTypeInteger
        new_x = mapping[MOI.VariableIndex(ci.value)]
        mapping[ci] = typeof(ci)(new_x.value)
        push!(dest.binaries, _info(dest, mapping[ci]))
        has_integrality = true
    end
    for ci in _constraints(src, MOI.VariableIndex, MOI.Integer)
        info = _info(dest, ci)
        info.type = _TYPE_INTEGER
        integrality[info.column+1] = kHighsVarTypeInteger
        new_x = mapping[MOI.VariableIndex(ci.value)]
        mapping[ci] = typeof(ci)(new_x.value)
        has_integrality = true
    end
    # Semicontinuous sets
    for ci in _constraints(src, MOI.VariableIndex, MOI.Semicontinuous{Float64})
        x = MOI.VariableIndex(ci.value)
        info = _info(dest, x)
        column = info.column + 1
        set = MOI.get(src, MOI.ConstraintSet(), ci)
        collower[column], colupper[column] = set.lower, set.upper
        info.bound, info.type = _BOUND_SEMI_CONTINUOUS, _TYPE_SEMI_CONTINUOUS
        info.lower, info.upper = set.lower, set.upper
        integrality[column] = kHighsVarTypeSemiContinuous
        mapping[ci] = typeof(ci)(mapping[x].value)
        has_integrality = true
    end
    # Semiinteger sets
    for ci in _constraints(src, MOI.VariableIndex, MOI.Semiinteger{Float64})
        x = MOI.VariableIndex(ci.value)
        info = _info(dest, x)
        column = info.column + 1
        set = MOI.get(src, MOI.ConstraintSet(), ci)
        collower[column], colupper[column] = set.lower, set.upper
        info.bound, info.type = _BOUND_SEMI_INTEGER, _TYPE_SEMI_INTEGER
        info.lower, info.upper = set.lower, set.upper
        integrality[column] = kHighsVarTypeSemiInteger
        mapping[ci] = typeof(ci)(mapping[x].value)
        has_integrality = true
    end
    # Build the model
    is_max = MOI.get(src, MOI.ObjectiveSense()) == MOI.MAX_SENSE
    dest.is_feasibility = false
    dest.is_objective_sense_set = true
    ret = Highs_passModel(
        dest,
        numcol,
        numrow,
        length(A.nzval),
        0,          # q_num_nz,
        kHighsMatrixFormatColwise,   # a_format,
        0,          # q_format,
        is_max ? kHighsObjSenseMaximize : kHighsObjSenseMinimize,
        0.0,                 # The objective function is handled separately at
        fill(0.0, numcol),  # the end.
        collower,
        colupper,
        rowlower,
        rowupper,
        A.colptr .- HighsInt(1),
        A.rowval .- HighsInt(1),
        A.nzval,
        C_NULL,  # qstart,  # The objective function is handled separately at
        C_NULL,  # qindex,  # the end.
        C_NULL,  # qvalue,
        has_integrality ? integrality : C_NULL,
    )
    _check_ret(ret)
    # Handle the objective function at the end.
    F = MOI.get(src, MOI.ObjectiveFunctionType())
    f_obj = MOI.get(src, MOI.ObjectiveFunction{F}())
    MOI.set(
        dest,
        MOI.ObjectiveFunction{F}(),
        MOI.Utilities.map_indices(mapping, f_obj),
    )
    # VariablePrimalStart
    start_attr = MOI.VariablePrimalStart()
    if start_attr in MOI.get(src, MOI.ListOfVariableAttributesSet())
        for x in MOI.get(src, MOI.ListOfVariableIndices())
            start = MOI.get(src, start_attr, x)
            MOI.set(dest, start_attr, mapping[x], start)
        end
    end
    return mapping
end

"""
    CallbackFunction(
        callback_types::Vector{Cint} = Cint[
            # kHighsCallbackLogging,
            kHighsCallbackSimplexInterrupt,
            kHighsCallbackIpmInterrupt,
            kHighsCallbackMipSolution,
            kHighsCallbackMipImprovingSolution,
            # kHighsCallbackMipLogging,
            kHighsCallbackMipInterrupt,
            kHighsCallbackMipGetCutPool,
            kHighsCallbackMipDefineLazyConstraints,
        ],
    ) <: MOI.AbstractModelAttribute

An `MOI.AbstractModelAttribute` for setting the callback function in HiGHS.

`callback_types` is a vector of places at which HiGHS will call the callback
function. By default, this is everything except the `Logging` callbacks.

## Callback signature

The function you provide as a callback must have the signature:
```julia
function user_callback(
    callback_type::Cint,
    message::Ptr{Cchar},
    data_out::HighsCallbackDataOut,
)::Cint
    return user_interrupt
end
```

The input arguments are interpreted as follows:

 * `callback_type` is a `Cint` corresponding to the location of where the
   callback is being called from
 * `message` is a pointer with the logging string (if the callback type is a
   `Logging` callback). Get a Julia `String` with `unsafe_string(message)`.
 * `data_out` is a `HiGHS.HighsCallbackDataOut` struct. This has various
   information fields, only some of which are filled.

The return argument is a `Cint`. If it is non-zero and the callback type is an
`Interrupt`, then HiGHS will interupt the solve.
"""
struct CallbackFunction <: MOI.AbstractModelAttribute
    callback_types::Vector{Cint}
end

function CallbackFunction()
    # By default, turn everything on except logging
    return CallbackFunction([
        # kHighsCallbackLogging,
        kHighsCallbackSimplexInterrupt,
        kHighsCallbackIpmInterrupt,
        kHighsCallbackMipSolution,
        kHighsCallbackMipImprovingSolution,
        # kHighsCallbackMipLogging,
        kHighsCallbackMipInterrupt,
        kHighsCallbackMipGetCutPool,
        kHighsCallbackMipDefineLazyConstraints,
    ])
end

# HiGHS v1.10 and HiGHS v1.11 introduced additional fields. We maintain this for
# backwards compatibility.
function HighsCallbackDataIn(terminate)
    return HighsCallbackDataIn(terminate, C_NULL, C_NULL, 0, zero(HighsInt))
end

function _cfn_user_callback(
    callback_type::Cint,
    message::Ptr{Cchar},
    p_data_out::Ptr{HighsCallbackDataOut},
    p_data_in::Ptr{HighsCallbackDataIn},
    p_user_data::Ptr{Cvoid},
)
    user_data = unsafe_pointer_to_objref(p_user_data)::_CallbackData
    data_out = unsafe_load(p_data_out)::HighsCallbackDataOut
    if callback_type in (
        kHighsCallbackSimplexInterrupt,
        kHighsCallbackIpmInterrupt,
        kHighsCallbackMipInterrupt,
    )
        @assert p_data_in !== C_NULL
        try
            reenable_sigint() do
                terminate = user_data.f(callback_type, message, data_out)
                unsafe_store!(p_data_in, HighsCallbackDataIn(terminate))
                return
            end
        catch err
            unsafe_store!(p_data_in, HighsCallbackDataIn(Cint(1)))
            if !(err isa InterruptException)
                rethrow(err)
            end
        end
    else
        # Ignore what the user says about terminating
        _ = user_data.f(callback_type, message, data_out)
    end
    return
end

function MOI.set(model::Optimizer, attr::CallbackFunction, f::Function)
    callback_cfn = @cfunction(
        _cfn_user_callback,
        Cvoid,
        (
            Cint,
            Ptr{Cchar},
            Ptr{HighsCallbackDataOut},
            Ptr{HighsCallbackDataIn},
            Ptr{Cvoid},
        ),
    )
    model.callback_data = _CallbackData(f)
    ret = Highs_setCallback(model, callback_cfn, model.callback_data)
    _check_ret(ret)
    for callback_type in attr.callback_types
        ret = Highs_startCallback(model, callback_type)
        _check_ret(ret)
    end
    return
end

function MOI.set(model::Optimizer, attr::CallbackFunction, ::Nothing)
    model.callback_data = nothing
    ret = Highs_setCallback(model, C_NULL, C_NULL)
    _check_ret(ret)
    for callback_type in attr.callback_types
        ret = Highs_stopCallback(model, callback_type)
        _check_ret(ret)
    end
    return
end

function MOI.write_to_file(model::Optimizer, filename::String)
    ret = Highs_writeModel(model, filename)
    _check_ret(ret)
    return
end

# These enums are deprecated. Use the `kHighsXXX` constants defined in
# libhighs.jl instead.

@enum(HighsBasisStatus, kLower = 0, kBasic, kUpper, kZero, kNonbasic)
@enum(HighsHessianFormat, kNoneHessian = 0, kTriangular, kSquare)
@enum(HighsMatrixFormat, kColwise = 1, kRowwise, kRowwisePartitioned)
@enum(
    HighsModelStatus,
    kNotset = 0,
    kLoadError,
    kModelError,
    kPresolveError,
    kSolveError,
    kPostsolveError,
    kModelEmpty,
    kOptimal,
    kInfeasible,
    kUnboundedOrInfeasible,
    kUnbounded,
    kObjectiveBound,
    kObjectiveTarget,
    kTimeLimit,
    kIterationLimit,
    kUnknown,
)
@enum(HighsObjSense, kMinimize = 1, kMaximize = -1)
@enum(HighsVartype, kContinuous = 0, kInteger = 1, kImplicitInteger = 2)
@enum(HighsStatus, HighsStatuskError = -1, HighsStatuskOk, HighsStatuskWarning)
