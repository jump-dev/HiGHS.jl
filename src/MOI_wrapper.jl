import MathOptInterface
import SparseArrays

const MOI = MathOptInterface
const CleverDicts = MOI.Utilities.CleverDicts

"""
    HighsBasisStatus

https://github.com/ERGO-Code/HiGHS/blob/4a5dd7499522f1fa730a31c59bba419b2bcc6839/src/lp_data/HConst.h#L151-L159
"""
@enum(HighsBasisStatus, kLower = 0, kBasic, kUpper, kZero, kNonbasic)

"""
    HighsMatrixFormat

https://github.com/ERGO-Code/HiGHS/blob/4a5dd7499522f1fa730a31c59bba419b2bcc6839/src/lp_data/HConst.h#L88
"""
@enum(HighsMatrixFormat, kNone = 0, kColwise, kRowwise)

"""
    HighsModelStatus

https://github.com/ERGO-Code/HiGHS/blob/4a5dd7499522f1fa730a31c59bba419b2bcc6839/src/lp_data/HConst.h#L126-L148
"""
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

"""
    HighsObjSense

https://github.com/ERGO-Code/HiGHS/blob/4a5dd7499522f1fa730a31c59bba419b2bcc6839/src/lp_data/HConst.h#L86
"""
@enum(HighsObjSense, kMinimize = 1, kMaximize = -1)

"""
    HighsVartype

    https://github.com/ERGO-Code/HiGHS/blob/4a5dd7499522f1fa730a31c59bba419b2bcc6839/src/lp_data/HConst.h#L69-L73
"""
@enum(HighsVartype, kContinuous = 0, kInteger = 1, kImplicitInteger = 2)

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
)

@enum(_TypeEnum, _TYPE_CONTINUOUS, _TYPE_INTEGER, _TYPE_BINARY,)

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
    column::Cint
    # Storage to keep track of the variable bounds.
    bound::_BoundEnum
    lower::Float64
    upper::Float64
    # Track integrality
    type::_TypeEnum
    # We can perform an optimization and only store two strings for the
    # constraint names because, at most, there can be two SingleVariable
    # constraints, e.g., LessThan, GreaterThan.
    lessthan_name::String
    greaterthan_interval_or_equalto_name::String
    integrality_name::String

    function _VariableInfo(
        index::MOI.VariableIndex,
        column::Cint,
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
            "",
            "",
            "",
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
        x::MOI.VariableIndex -> x.value,
        x::Int64 -> MOI.VariableIndex(x),
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
    row::Cint
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
        x::_ConstraintKey -> x.value,
        x::Int64 -> _ConstraintKey(x),
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
    colvalue::Vector{Cdouble}
    coldual::Vector{Cdouble}
    colstatus::Vector{Cint}
    rowvalue::Vector{Cdouble}
    rowdual::Vector{Cdouble}
    rowstatus::Vector{Cint}
    has_primal_solution::Bool
    has_dual_solution::Bool
    has_primal_ray::Bool
    has_dual_ray::Bool
    function _Solution()
        return new(
            _OPTIMIZE_NOT_CALLED,
            Cdouble[],
            Cdouble[],
            Cint[],
            Cdouble[],
            Cdouble[],
            Cint[],
            false,
            false,
            false,
            false,
        )
    end
end

function Base.empty!(x::_Solution)
    x.status = _OPTIMIZE_NOT_CALLED
    empty!(x.colvalue)
    empty!(x.coldual)
    empty!(x.colstatus)
    empty!(x.rowvalue)
    empty!(x.rowdual)
    empty!(x.rowstatus)
    x.has_primal_solution = false
    x.has_dual_solution = false
    x.has_dual_ray = false
    x.has_primal_ray = false
    return x
end

Base.isempty(x::_Solution) = x.status == _OPTIMIZE_NOT_CALLED

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

    """
        Optimizer()

    Create a new Optimizer object.
    """
    function Optimizer()
        model = new(
            C_NULL,
            Dict{String,Any}(),
            "",
            true,
            Set{_VariableInfo}(),
            _variable_info_dict(),
            _constraint_info_dict(),
            nothing,
            nothing,
            _Solution(),
        )
        MOI.empty!(model)
        finalizer(Highs_destroy, model)
        return model
    end
end

Base.cconvert(::Type{Ptr{Cvoid}}, model::Optimizer) = model
Base.unsafe_convert(::Type{Ptr{Cvoid}}, model::Optimizer) = model.inner

function _check_ret(ret::Cint)
    if ret != Cint(0)
        error(
            "Encountered an error in HiGHS (Status $(ref)). Check the log " *
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
        MOI.set(model, MOI.RawParameter(key), value)
    end
    model.is_feasibility = true
    empty!(model.binaries)
    empty!(model.variable_info)
    empty!(model.affine_constraint_info)
    model.name_to_variable = nothing
    model.name_to_constraint_index = nothing
    empty!(model.solution)
    return
end

function MOI.is_empty(model::Optimizer)
    offset = Ref{Cdouble}(0.0)
    ret = Highs_getObjectiveOffset(model, offset)
    _check_ret(ret)
    return Highs_getNumCols(model) == 0 &&
           Highs_getNumRows(model) == 0 &&
           model.is_feasibility &&
           isempty(model.binaries) &&
           isempty(model.variable_info) &&
           isempty(model.affine_constraint_info) &&
           model.name_to_variable === nothing &&
           model.name_to_constraint_index === nothing &&
           isempty(model.solution) &&
           iszero(offset[])
end

MOI.get(::Optimizer, ::MOI.SolverName) = "HiGHS"

MOI.get(model::Optimizer, ::MOI.RawSolver) = model

function MOI.get(::Optimizer, ::MOI.ListOfVariableAttributesSet)
    return MOI.AbstractVariableAttribute[MOI.VariableName()]
end

function MOI.get(model::Optimizer, ::MOI.ListOfModelAttributesSet)
    attributes = [
        MOI.ObjectiveSense(),
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
    ]
    if MOI.get(model, MOI.Name()) != ""
        push!(attributes, MOI.Name())
    end
    return attributes
end

function MOI.get(::Optimizer, ::MOI.ListOfConstraintAttributesSet)
    return MOI.AbstractConstraintAttribute[MOI.ConstraintName()]
end

function MOI.get(model::Optimizer, ::MOI.NumberOfConstraints{F,S}) where {F,S}
    # TODO: this could be more efficient.
    return length(MOI.get(model, MOI.ListOfConstraintIndices{F,S}()))
end

function MOI.get(model::Optimizer, ::MOI.ListOfConstraints)
    constraints = Set{Tuple{DataType,DataType}}()
    for info in values(model.variable_info)
        if info.bound == _BOUND_NONE
        elseif info.bound == _BOUND_LESS_THAN
            push!(constraints, (MOI.SingleVariable, MOI.LessThan{Float64}))
        elseif info.bound == _BOUND_GREATER_THAN
            push!(constraints, (MOI.SingleVariable, MOI.GreaterThan{Float64}))
        elseif info.bound == _BOUND_LESS_AND_GREATER_THAN
            push!(constraints, (MOI.SingleVariable, MOI.LessThan{Float64}))
            push!(constraints, (MOI.SingleVariable, MOI.GreaterThan{Float64}))
        elseif info.bound == _BOUND_EQUAL_TO
            push!(constraints, (MOI.SingleVariable, MOI.EqualTo{Float64}))
        else
            @assert info.bound == _BOUND_INTERVAL
            push!(constraints, (MOI.SingleVariable, MOI.Interval{Float64}))
        end
        if info.type == _TYPE_INTEGER
            push!(constraints, (MOI.SingleVariable, MOI.Integer))
        elseif info.type == _TYPE_BINARY
            push!(constraints, (MOI.SingleVariable, MOI.ZeroOne))
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
### MOI.RawParameter
###

function MOI.supports(model::Optimizer, param::MOI.RawParameter)
    if !(param.name isa String)
        return false
    end
    typeP = Ref{Cint}()
    return Highs_getOptionType(model, param.name, typeP) == 0
end

function _check_option_status(ret::Cint)
    if ret != 0
        error("Encountered an error in HiGHS: Check the log for details.")
    end
    return
end

function _set_option(model::Optimizer, option::String, value::Bool)
    return Highs_setBoolOptionValue(model, option, Cint(value))
end

function _set_option(model::Optimizer, option::String, value::Integer)
    return Highs_setIntOptionValue(model, option, Cint(value))
end

function _set_option(model::Optimizer, option::String, value::AbstractFloat)
    return Highs_setDoubleOptionValue(model, option, Cdouble(value))
end

function _set_option(model::Optimizer, option::String, value::String)
    return Highs_setStringOptionValue(model, option, value)
end

function MOI.set(model::Optimizer, param::MOI.RawParameter, value)
    if !MOI.supports(model, param)
        throw(MOI.UnsupportedAttribute(param))
    end
    model.options[param.name] = value
    ret = _set_option(model, param.name, value)
    return _check_option_status(ret)
end

### MOI.get

function _get_bool_option(model::Optimizer, option::String)
    value = Ref{Cint}(0)
    ret = Highs_getBoolOptionValue(model, option, value)
    _check_option_status(ret)
    return Bool(value[])
end

function _get_int_option(model::Optimizer, option::String)
    value = Ref{Cint}(0)
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
    buffer = Vector{Cchar}(undef, 1024)
    bufferP = pointer(buffer)
    GC.@preserve buffer begin
        ret = Highs_getStringOptionValue(model, option, bufferP)
        _check_option_status(ret)
        return unsafe_string(bufferP)
    end
end

function MOI.get(model::Optimizer, param::MOI.RawParameter)
    if !(param.name isa String)
        throw(MOI.UnsupportedAttribute(param))
    end
    typeP = Ref{Cint}()
    ret = Highs_getOptionType(model, param.name, typeP)
    if ret != 0
        throw(MOI.UnsupportedAttribute(param))
    elseif typeP[] == 0
        return _get_bool_option(model, param.name)::Bool
    elseif typeP[] == 1
        return _get_int_option(model, param.name)::Cint
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
    return MOI.set(model, MOI.RawParameter("time_limit"), Inf)
end

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, limit::Real)
    return MOI.set(model, MOI.RawParameter("time_limit"), Float64(limit))
end

function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    return MOI.get(model, MOI.RawParameter("time_limit"))
end

###
### MOI.Silent
###

MOI.supports(::Optimizer, ::MOI.Silent) = true

function MOI.get(model::Optimizer, ::MOI.Silent)
    return !MOI.get(model, MOI.RawParameter("output_flag"))
end

function MOI.set(model::Optimizer, ::MOI.Silent, flag::Bool)
    MOI.set(model, MOI.RawParameter("output_flag"), !flag)
    return
end

###
### MOI.Name
###

MOI.supports(::Optimizer, ::MOI.Name) = true

MOI.get(model::Optimizer, ::MOI.Name) = model.name

MOI.set(model::Optimizer, ::MOI.Name, name::String) = (model.name = name)

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
    if haskey(model.variable_info, key)
        return model.variable_info[key]
    end
    return throw(MOI.InvalidIndex(key))
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
        _VariableInfo(MOI.VariableIndex(0), Cint(0)),
    )
    info = _info(model, index)
    # Now, set `.index` and `.column`.
    info.index = index
    info.column = Cint(length(model.variable_info) - 1)
    ret = Highs_addCol(model, 0.0, -Inf, Inf, 0, C_NULL, C_NULL)
    _check_ret(ret)
    return index
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
    if haskey(model.name_to_variable, name)
        variable = model.name_to_variable[name]
        if variable === nothing
            error("Duplicate name detected: $(name)")
        end
        return variable
    end
    return nothing
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
    info.name = name
    model.name_to_variable = nothing
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

###
### Objectives
###

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveSense,
    sense::MOI.OptimizationSense,
)
    x = sense == MOI.MAX_SENSE ? Cint(kMaximize) : Cint(kMinimize)
    ret = Highs_changeObjectiveSense(model, x)
    _check_ret(ret)
    if sense == MOI.FEASIBILITY_SENSE
        model.is_feasibility = true
        # TODO(odow): cache the mask.
        n = MOI.get(model, MOI.NumberOfVariables())
        ret =
            Highs_changeColsCostByMask(model, ones(Cint, n), zeros(Cdouble, n))
        _check_ret(ret)
        ret = Highs_changeObjectiveOffset(model, 0.0)
        _check_ret(ret)
    else
        model.is_feasibility = false
    end
    return
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveSense)
    if model.is_feasibility
        return MOI.FEASIBILITY_SENSE
    end
    senseP = Ref{Cint}()
    ret = Highs_getObjectiveSense(model, senseP)
    _check_ret(ret)
    return senseP[] == Cint(kMinimize) ? MOI.MIN_SENSE : MOI.MAX_SENSE
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
)
    return true
end

function MOI.get(::Optimizer, ::MOI.ObjectiveFunctionType)
    return MOI.ScalarAffineFunction{Float64}
end

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
    f::MOI.ScalarAffineFunction{Float64},
)
    num_vars = length(model.variable_info)
    obj = zeros(Float64, num_vars)
    for term in f.terms
        col = column(model, term.variable_index)
        obj[col+1] += term.coefficient
    end
    # TODO(odow): cache the mask.
    mask = ones(Cint, num_vars)
    ret = Highs_changeColsCostByMask(model, mask, obj)
    _check_ret(ret)
    ret = Highs_changeObjectiveOffset(model, f.constant)
    _check_ret(ret)
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
)
    ncols = MOI.get(model, MOI.NumberOfVariables())
    offset = Ref{Cdouble}(0.0)
    ret = Highs_getObjectiveOffset(model, offset)
    _check_ret(ret)
    if ncols == 0
        return MOI.ScalarAffineFunction{Float64}(
            MOI.ScalarAffineTerm{Float64}[],
            offset[],
        )
    end
    num_colsP, nnzP = Ref{Cint}(0), Ref{Cint}(0)
    costs = Vector{Cdouble}(undef, ncols)
    ret = Highs_getColsByRange(
        model,
        0,
        ncols - 1,
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

function MOI.get(model::Optimizer, ::MOI.ObjectiveFunction{F}) where {F}
    obj = MOI.get(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
    )
    return convert(F, obj)
end

function MOI.modify(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
    chg::MOI.ScalarConstantChange{Float64},
)
    ret = Highs_changeObjectiveOffset(model, chg.new_constant)
    _check_ret(ret)
    return
end

function MOI.modify(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
    chg::MOI.ScalarCoefficientChange{Float64},
)
    ret = Highs_changeColCost(
        model,
        column(model, chg.variable),
        chg.new_coefficient,
    )
    _check_ret(ret)
    return
end

###
### SingleVariable-in-Set constraints.
###

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.SingleVariable},
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
    ::MOI.ListOfConstraintIndices{MOI.SingleVariable,S},
) where {S<:_SCALAR_SETS}
    indices = MOI.ConstraintIndex{MOI.SingleVariable,S}[
        MOI.ConstraintIndex{MOI.SingleVariable,S}(key.value) for
        (key, info) in model.variable_info if info.bound in _bound_enums(S)
    ]
    return sort!(indices, by = x -> x.value)
end

function _info(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable,<:Any},
)
    var_index = MOI.VariableIndex(c.value)
    if haskey(model.variable_info, var_index)
        return _info(model, var_index)
    end
    return throw(MOI.InvalidIndex(c))
end

"""
    column(
        model::Optimizer,
        c::MOI.ConstraintIndex{MOI.SingleVariable,<:Any},
    )

Return the 0-indexed column associated with the variable bounds `c` in `model`.
"""
function column(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable,<:Any},
)
    return _info(model, c).column
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable,MOI.LessThan{Float64}},
)
    if haskey(model.variable_info, MOI.VariableIndex(c.value))
        info = _info(model, c)
        return info.bound == _BOUND_LESS_THAN ||
               info.bound == _BOUND_LESS_AND_GREATER_THAN
    end
    return false
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable,MOI.GreaterThan{Float64}},
)
    if haskey(model.variable_info, MOI.VariableIndex(c.value))
        info = _info(model, c)
        return info.bound == _BOUND_GREATER_THAN ||
               info.bound == _BOUND_LESS_AND_GREATER_THAN
    end
    return false
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable,MOI.Interval{Float64}},
)
    return haskey(model.variable_info, MOI.VariableIndex(c.value)) &&
           _info(model, c).bound == _BOUND_INTERVAL
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable,MOI.EqualTo{Float64}},
)
    return haskey(model.variable_info, MOI.VariableIndex(c.value)) &&
           _info(model, c).bound == _BOUND_EQUAL_TO
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.SingleVariable,<:Any},
)
    MOI.throw_if_not_valid(model, c)
    return MOI.SingleVariable(MOI.VariableIndex(c.value))
end

function MOI.set(
    ::Optimizer,
    ::MOI.ConstraintFunction,
    ::MOI.ConstraintIndex{MOI.SingleVariable,<:Any},
    ::MOI.SingleVariable,
)
    return throw(MOI.SettingSingleVariableFunctionNotAllowed())
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
    end
    return
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.SingleVariable,
    s::S,
) where {S<:_SCALAR_SETS}
    info = _info(model, f.variable)
    _update_info(info, s)
    index = MOI.ConstraintIndex{MOI.SingleVariable,typeof(s)}(f.variable.value)
    MOI.set(model, MOI.ConstraintSet(), index, s)
    return index
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable,MOI.LessThan{Float64}},
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
    info.lessthan_name = ""
    model.name_to_constraint_index = nothing
    return
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable,MOI.GreaterThan{Float64}},
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
    info.greaterthan_interval_or_equalto_name = ""
    model.name_to_constraint_index = nothing
    return
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable,MOI.Interval{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    ret = Highs_changeColBounds(model, info.column, -Inf, Inf)
    _check_ret(ret)
    info.lower, info.upper = -Inf, Inf
    info.bound = _BOUND_NONE
    info.greaterthan_interval_or_equalto_name = ""
    model.name_to_constraint_index = nothing
    return
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.SingleVariable,MOI.EqualTo{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    ret = Highs_changeColBounds(model, info.column, -Inf, Inf)
    _check_ret(ret)
    info.lower, info.upper = -Inf, Inf
    info.bound = _BOUND_NONE
    info.greaterthan_interval_or_equalto_name = ""
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable,MOI.GreaterThan{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    return MOI.GreaterThan(_info(model, c).lower)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable,MOI.LessThan{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    return MOI.LessThan(_info(model, c).upper)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable,MOI.EqualTo{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    return MOI.EqualTo(_info(model, c).lower)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable,MOI.Interval{Float64}},
)
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    return MOI.Interval(info.lower, info.upper)
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable,S},
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

function MOI.supports(
    ::Optimizer,
    ::MOI.ConstraintName,
    ::Type{<:MOI.ConstraintIndex},
)
    return true
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.SingleVariable,S},
) where {S<:_SCALAR_SETS}
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    if S <: MOI.LessThan
        return info.lessthan_name
    else
        return info.greaterthan_interval_or_equalto_name
    end
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.SingleVariable,S},
    name::String,
) where {S<:_SCALAR_SETS}
    MOI.throw_if_not_valid(model, c)
    info = _info(model, c)
    if S <: MOI.LessThan
        info.lessthan_name = name
    else
        info.greaterthan_interval_or_equalto_name = name
    end
    model.name_to_constraint_index = nothing
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
    if haskey(model.affine_constraint_info, key)
        return model.affine_constraint_info[key]
    end
    return throw(MOI.InvalidIndex(c))
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
    indices::Vector{Cint},
    coefficients::Vector{Float64},
    model::Optimizer,
    f::MOI.ScalarAffineFunction{Float64},
)
    i = 1
    for term in f.terms
        indices[i] = column(model, term.variable_index)
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
    indices, coefficients = zeros(Cint, nnz), zeros(Cdouble, nnz)
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
        Cint(length(model.affine_constraint_info) - 1)
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

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
) where {S<:_SCALAR_SETS}
    r = row(model, c)
    num_row = Ref{Cint}(0)
    num_nz = Ref{Cint}(0)
    matrix_start = Ref{Cint}(0)
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
    matrix_index = Vector{Cint}(undef, num_nz[])
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
    info.name = name
    model.name_to_constraint_index = nothing
    return
end

function MOI.get(model::Optimizer, ::Type{MOI.ConstraintIndex}, name::String)
    if model.name_to_constraint_index === nothing
        _rebuild_name_to_constraint_index(model)
    end
    if haskey(model.name_to_constraint_index, name)
        constr = model.name_to_constraint_index[name]
        if constr === nothing
            error("Duplicate constraint name detected: $(name)")
        end
        return constr
    end
    return nothing
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
    for (key, info) in model.variable_info
        if !isempty(info.lessthan_name)
            _set_name_to_constraint_index(
                model,
                info.lessthan_name,
                MOI.ConstraintIndex{MOI.SingleVariable,MOI.LessThan{Float64}}(
                    key.value,
                ),
            )
        end
        if !isempty(info.greaterthan_interval_or_equalto_name)
            S = if info.bound == _BOUND_GREATER_THAN
                MOI.GreaterThan{Float64}
            elseif info.bound == _BOUND_LESS_AND_GREATER_THAN
                MOI.GreaterThan{Float64}
            elseif info.bound == _BOUND_EQUAL_TO
                MOI.EqualTo{Float64}
            else
                @assert info.bound == _BOUND_INTERVAL
                MOI.Interval{Float64}
            end
            _set_name_to_constraint_index(
                model,
                info.greaterthan_interval_or_equalto_name,
                MOI.ConstraintIndex{MOI.SingleVariable,S}(key.value),
            )
        end
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
    terms = Dict(x.variable_index => 0.0 for x in old.terms)
    for term in f.terms
        terms[term.variable_index] = term.coefficient
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
    _store_solution(model::Optimizer, ret::Cint)

Get the solution from a run of HiGHS.

`ret` is the `HighsStatus` enum {OK = 0, Warning, Error} from `Highs_run`.
"""
function _store_solution(model::Optimizer, ret::Cint)
    x = model.solution
    x.status = ret == 0 ? _OPTIMIZE_OK : _OPTIMIZE_ERRORED
    x.has_primal_solution = false
    x.has_dual_solution = false
    x.has_dual_ray = false
    x.has_primal_ray = false
    numCols = Highs_getNumCols(model)
    numRows = Highs_getNumRows(model)
    resize!(x.colvalue, numCols)
    resize!(x.coldual, numCols)
    resize!(x.colstatus, numCols)
    resize!(x.rowvalue, numRows)
    resize!(x.rowdual, numRows)
    resize!(x.rowstatus, numRows)
    status = HighsModelStatus(Highs_getModelStatus(model))
    statusP = Ref{Cint}()
    if status == kInfeasible
        ret = Highs_getDualRay(model, statusP, x.rowdual)
        _check_ret(ret)
        x.has_dual_ray = statusP[] == 1
    elseif status == kUnbounded
        ret = Highs_getPrimalRay(model, statusP, x.colvalue)
        _check_ret(ret)
        x.has_primal_ray = statusP[] == 1
    else
        Highs_getIntInfoValue(model, "primal_solution_status", statusP)
        x.has_primal_solution = statusP[] == 2
        Highs_getIntInfoValue(model, "dual_solution_status", statusP)
        x.has_dual_solution = statusP[] == 2
        if x.has_primal_solution
            Highs_getSolution(
                model,
                x.colvalue,
                x.coldual,
                x.rowvalue,
                x.rowdual,
            )
            Highs_getBasis(model, x.colstatus, x.rowstatus)
        end
    end
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
    ret = Highs_run(model)
    _store_solution(model, ret)
    for info in model.binaries
        Highs_changeColBounds(model, info.column, info.lower, info.upper)
    end
    return
end

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if model.solution.status == _OPTIMIZE_NOT_CALLED
        return MOI.OPTIMIZE_NOT_CALLED
    elseif model.solution.status == _OPTIMIZE_ERRORED
        return MOI.OTHER_ERROR
    end
    status = HighsModelStatus(Highs_getModelStatus(model))
    if status == kNotset
        return MOI.OTHER_ERROR
    elseif status == kLoadError
        return MOI.OTHER_ERROR
    elseif status == kModelError
        return MOI.INVALID_MODEL
    elseif status == kPresolveError
        return MOI.OTHER_ERROR
    elseif status == kSolveError
        return MOI.OTHER_ERROR
    elseif status == kPostsolveError
        return MOI.OTHER_ERROR
    elseif status == kModelEmpty
        return MOI.INVALID_MODEL
    elseif status == kOptimal
        return MOI.OPTIMAL
    elseif status == kInfeasible
        return MOI.INFEASIBLE
    elseif status == kUnboundedOrInfeasible
        return MOI.INFEASIBLE_OR_UNBOUNDED
    elseif status == kUnbounded
        return MOI.DUAL_INFEASIBLE
    elseif status == kObjectiveBound
        return MOI.OBJECTIVE_LIMIT
    elseif status == kObjectiveTarget
        return MOI.OBJECTIVE_LIMIT
    elseif status == kTimeLimit
        return MOI.TIME_LIMIT
    elseif status == kIterationLimit
        return MOI.ITERATION_LIMIT
    else
        @assert status == kUnknown
        return MOI.OTHER_ERROR
    end
end

function MOI.get(model::Optimizer, ::MOI.ResultCount)
    x = model.solution
    if x.has_primal_solution || x.has_dual_ray || x.has_primal_ray
        return 1
    else
        return 0
    end
end

function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    status = HighsModelStatus(Highs_getModelStatus(model))
    return string(status)
end

function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    if attr.N != 1
        return MOI.NO_SOLUTION
    elseif model.solution.has_primal_solution
        return MOI.FEASIBLE_POINT
    elseif model.solution.has_primal_ray
        return MOI.INFEASIBILITY_CERTIFICATE
    end
    return MOI.NO_SOLUTION
end

function MOI.get(model::Optimizer, attr::MOI.DualStatus)
    if attr.N != 1
        return MOI.NO_SOLUTION
    elseif model.solution.has_dual_solution
        return MOI.FEASIBLE_POINT
    elseif model.solution.has_dual_ray
        return MOI.INFEASIBILITY_CERTIFICATE
    end
    return MOI.NO_SOLUTION
end

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    return Highs_getObjectiveValue(model)
end

function MOI.get(model::Optimizer, attr::MOI.DualObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    return MOI.Utilities.get_fallback(model, attr, Float64)
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveBound)
    # TODO(odow): there is a bug in HiGHS where it reports incorrect values for
    # mip_dual_bound in the LP case. As a work-around, just return the most
    # optimistic of the primal and dual values.
    p = Ref{Cdouble}()
    Highs_getDoubleInfoValue(model, "mip_dual_bound", p)
    sense = _sense_corrector(model)
    primal = Highs_getObjectiveValue(model)
    return sense == -1 ? max(primal, -p[]) : min(primal, p[])
end

function MOI.get(model::Optimizer, ::MOI.SolveTime)
    return Highs_getRunTime(model)
end

function MOI.get(model::Optimizer, ::MOI.SimplexIterations)
    return Highs_getSimplexIterationCount(model)
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
    c::MOI.ConstraintIndex{MOI.SingleVariable,<:_SCALAR_SETS},
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
    return model.solution.rowvalue[row(model, c)+1]
end

"""
    _sense_corrector(model::Optimizer)

Return `-1` if `MAX_SENSE`, and `1` is `MIN_SENSE`. Useful for correcting
solution attributes which are sense-dependent.
"""
function _sense_corrector(model::Optimizer)
    senseP = Ref{Cint}()
    ret = Highs_getObjectiveSense(model, senseP)
    _check_ret(ret)
    return senseP[]
end

"""
    _farkas_variable_dual(model::Optimizer, col::Cint)

Return a Farkas dual associated with the variable bounds of `col`. Given a dual
ray:

    ā * x = λ' * A * x <= λ' * b = -β + sum(āᵢ * Uᵢ | āᵢ < 0) + sum(āᵢ * Lᵢ | āᵢ > 0)

The Farkas dual of the variable is ā, and it applies to the upper bound if ā < 0,
and it applies to the lower bound if ā > 0.
"""
function _farkas_variable_dual(model::Optimizer, col::Cint)
    num_nz, num_cols = Ref{Cint}(0), Ref{Cint}(0)
    # TODO(odow): how does getColsByRangeWork???
    m = Highs_getNumRows(model)
    matrix_start = zeros(Cint, 2)
    matrix_index = Vector{Cint}(undef, m)
    matrix_value = Vector{Cdouble}(undef, m)
    ret = Highs_getColsByRange(
        model,
        col,
        col,
        num_cols,
        C_NULL,
        C_NULL,
        C_NULL,
        num_nz,
        matrix_start,
        matrix_index,
        matrix_value,
    )
    _check_ret(ret)
    dual = 0.0
    for i in 1:num_nz[]
        dual += -model.solution.rowdual[matrix_index[i]+1] * matrix_value[i]
    end
    return dual
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
    _signed_dual(dual::Float64, ::Type{Set}, status::HighsBasisStatus)

Determine whether the dual of an interval constraint applies to the lower or
upper bound using the basis status reported by HiGHS.
"""
function _signed_dual(
    dual::Float64,
    ::Type{MOI.LessThan{Float64}},
    status::HighsBasisStatus,
)
    return status == kUpper ? dual : 0.0
end

function _signed_dual(
    dual::Float64,
    ::Type{MOI.GreaterThan{Float64}},
    status::HighsBasisStatus,
)
    return status == kLower ? dual : 0.0
end

_signed_dual(dual::Float64, ::Any, ::HighsBasisStatus) = dual

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable,S},
) where {S<:_SCALAR_SETS}
    MOI.check_result_index_bounds(model, attr)
    col = column(model, c)
    if model.solution.has_dual_ray[] == 1
        return _signed_dual(_farkas_variable_dual(model, col), S)
    end
    dual = _sense_corrector(model) * model.solution.coldual[col+1]
    stat = HighsBasisStatus(model.solution.colstatus[col+1])
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
    stat = HighsBasisStatus(model.solution.rowstatus[r])
    return _signed_dual(_sense_corrector(model) * dual, S, stat)
end

###
### MOI.ConstraintBasisStatus
###

# HiGHS only reports a single basis status for each ranged constraint. Therefore
# if the status is kLower or kUpper, we must distinguish whether or not it is
# refering to the given constraint.
function _nonbasic_status(code::HighsBasisStatus, ::Type{<:MOI.Interval})
    return code == kLower ? MOI.NONBASIC_AT_LOWER : MOI.NONBASIC_AT_UPPER
end
function _nonbasic_status(code::HighsBasisStatus, ::Type{<:MOI.LessThan})
    return code == kLower ? MOI.BASIC : MOI.NONBASIC
end
function _nonbasic_status(code::HighsBasisStatus, ::Type{<:MOI.GreaterThan})
    return code == kLower ? MOI.NONBASIC : MOI.BASIC
end
_nonbasic_status(::HighsBasisStatus, ::Type{<:MOI.EqualTo}) = MOI.NONBASIC

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintBasisStatus,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
) where {S<:_SCALAR_SETS}
    MOI.check_result_index_bounds(model, attr)
    stat = HighsBasisStatus(model.solution.rowstatus[row(model, c)+1])
    if stat == kLower
        return _nonbasic_status(stat, S)
    elseif stat == kBasic
        return MOI.BASIC
    elseif stat == kUpper
        return _nonbasic_status(stat, S)
    elseif stat == kZero || stat == kNonbasic
        return MOI.NONBASIC
    else
        @assert stat == SUPER
        return MOI.SUPER_BASIC
    end
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintBasisStatus,
    c::MOI.ConstraintIndex{MOI.SingleVariable,S},
) where {S<:_SCALAR_SETS}
    MOI.check_result_index_bounds(model, attr)
    stat = HighsBasisStatus(model.solution.colstatus[column(model, c)+1])
    if stat == kLower
        return _nonbasic_status(stat, S)
    elseif stat == kBasic
        return MOI.BASIC
    elseif stat == kUpper
        return _nonbasic_status(stat, S)
    elseif stat == kZero || stat == kNonbasic
        return MOI.NONBASIC
    else
        @assert stat == SUPER
        return MOI.SUPER_BASIC
    end
end

###
### Integrality
###

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.SingleVariable},
    ::Type{<:Union{MOI.ZeroOne,MOI.Integer}},
)
    return true
end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.SingleVariable,MOI.Integer},
)
    return _info(model, ci).type == _TYPE_INTEGER
end

function MOI.is_valid(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.SingleVariable,MOI.ZeroOne},
)
    return _info(model, ci).type == _TYPE_BINARY
end

function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintIndices{MOI.SingleVariable,MOI.Integer},
)
    indices = MOI.ConstraintIndex{MOI.SingleVariable,MOI.Integer}[
        MOI.ConstraintIndex{MOI.SingleVariable,MOI.Integer}(key.value) for
        (key, info) in model.variable_info if info.type == _TYPE_INTEGER
    ]
    return sort!(indices, by = x -> x.value)
end

function MOI.get(
    model::Optimizer,
    ::MOI.ListOfConstraintIndices{MOI.SingleVariable,MOI.ZeroOne},
)
    indices = MOI.ConstraintIndex{MOI.SingleVariable,MOI.ZeroOne}[
        MOI.ConstraintIndex{MOI.SingleVariable,MOI.ZeroOne}(key.value) for
        (key, info) in model.variable_info if info.type == _TYPE_BINARY
    ]
    return sort!(indices, by = x -> x.value)
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.SingleVariable,
    ::MOI.Integer,
)
    info = _info(model, f.variable)
    info.type = _TYPE_INTEGER
    ci = MOI.ConstraintIndex{MOI.SingleVariable,MOI.Integer}(f.variable.value)
    Highs_changeColIntegrality(model, info.column, kInteger)
    return ci
end

function MOI.add_constraint(
    model::Optimizer,
    f::MOI.SingleVariable,
    ::MOI.ZeroOne,
)
    info = _info(model, f.variable)
    info.type = _TYPE_BINARY
    ci = MOI.ConstraintIndex{MOI.SingleVariable,MOI.ZeroOne}(f.variable.value)
    Highs_changeColIntegrality(model, info.column, kInteger)
    push!(model.binaries, info)
    return ci
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.SingleVariable,MOI.Integer},
)
    MOI.throw_if_not_valid(model, ci)
    return MOI.Integer()
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.SingleVariable,MOI.ZeroOne},
)
    MOI.throw_if_not_valid(model, ci)
    return MOI.ZeroOne()
end

function MOI.delete(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.SingleVariable,MOI.Integer},
)
    MOI.throw_if_not_valid(model, ci)
    info = _info(model, ci)
    info.type = _TYPE_CONTINUOUS
    Highs_changeColIntegrality(model, info.column, kContinuous)
    return
end

function MOI.delete(
    model::Optimizer,
    ci::MOI.ConstraintIndex{MOI.SingleVariable,MOI.ZeroOne},
)
    MOI.throw_if_not_valid(model, ci)
    info = _info(model, ci)
    info.type = _TYPE_CONTINUOUS
    Highs_changeColIntegrality(model, info.column, kContinuous)
    delete!(model.binaries, info)
    return
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintName,
    ci::MOI.ConstraintIndex{
        MOI.SingleVariable,
        <:Union{MOI.Integer,MOI.ZeroOne},
    },
    name::String,
)
    MOI.throw_if_not_valid(model, ci)
    info = _info(model, ci)
    info.integrality_name = name
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintName,
    ci::MOI.ConstraintIndex{
        MOI.SingleVariable,
        <:Union{MOI.Integer,MOI.ZeroOne},
    },
)
    MOI.throw_if_not_valid(model, ci)
    info = _info(model, ci)
    return info.integrality_name
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

function _extract_bound_data(
    dest::Optimizer,
    src::MOI.ModelLike,
    mapping,
    collower::Vector{Float64},
    colupper::Vector{Float64},
    ::Type{S},
) where {S}
    for c_index in
        MOI.get(src, MOI.ListOfConstraintIndices{MOI.SingleVariable,S}())
        f = MOI.get(src, MOI.ConstraintFunction(), c_index)
        s = MOI.get(src, MOI.ConstraintSet(), c_index)
        new_f = mapping.varmap[f.variable]
        info = _info(dest, new_f)
        _add_bounds(collower, colupper, info.column + 1, s)
        _update_info(info, s)
        mapping.conmap[c_index] =
            MOI.ConstraintIndex{MOI.SingleVariable,S}(new_f.value)
    end
    return
end

function _copy_to_columns(dest::Optimizer, src::MOI.ModelLike, mapping)
    x_src = MOI.get(src, MOI.ListOfVariableIndices())
    numcols = Cint(length(x_src))
    for i in 1:numcols
        index = CleverDicts.add_item(
            dest.variable_info,
            _VariableInfo(MOI.VariableIndex(0), Cint(0)),
        )
        info = _info(dest, index)
        info.name = MOI.get(dest, MOI.VariableName(), x_src[i])
        info.index = index
        info.column = Cint(i - 1)
        mapping.varmap[x_src[i]] = index
    end
    fobj =
        MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    c = fill(0.0, numcols)
    for term in fobj.terms
        i = mapping.varmap[term.variable_index].value
        c[i] += term.coefficient
    end
    ret = Highs_changeObjectiveOffset(dest, fobj.constant)
    _check_ret(ret)
    return numcols, c
end

add_sizehint!(vec, n) = sizehint!(vec, length(vec) + n)

function _extract_row_data(
    dest::Optimizer,
    src::MOI.ModelLike,
    mapping,
    rowlower::Vector{Float64},
    rowupper::Vector{Float64},
    I::Vector{Cint},
    J::Vector{Cint},
    V::Vector{Float64},
    ::Type{S},
) where {S}
    row = length(I) == 0 ? 1 : I[end] + 1
    list = MOI.get(
        src,
        MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64},S}(),
    )
    numrows = length(list)
    add_sizehint!(rowlower, numrows)
    add_sizehint!(rowupper, numrows)
    n_terms = 0
    fs = Array{MOI.ScalarAffineFunction{Float64}}(undef, numrows)
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
        dest.affine_constraint_info[key].row =
            Cint(length(dest.affine_constraint_info) - 1)
        mapping.conmap[c_index] =
            MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S}(key.value)
    end
    add_sizehint!(I, n_terms)
    add_sizehint!(J, n_terms)
    add_sizehint!(V, n_terms)
    for (i, c_index) in enumerate(list)
        for term in fs[i].terms
            push!(I, row)
            push!(J, Cint(mapping.varmap[term.variable_index].value))
            push!(V, term.coefficient)
        end
        row += 1
    end
    return
end

function _check_input_data(dest::Optimizer, src::MOI.ModelLike)
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        if !MOI.supports_constraint(dest, F, S)
            throw(
                MOI.UnsupportedConstraint{F,S}(
                    "HiGHS does not support constraints of type $F-in-$S.",
                ),
            )
        end
    end
    fobj_type = MOI.get(src, MOI.ObjectiveFunctionType())
    if !MOI.supports(dest, MOI.ObjectiveFunction{fobj_type}())
        throw(MOI.UnsupportedAttribute(MOI.ObjectiveFunction(fobj_type)))
    end
    return
end

MOI.Utilities.supports_default_copy_to(::Optimizer, ::Bool) = true

function MOI.copy_to(
    dest::Optimizer,
    src::MOI.ModelLike;
    copy_names::Bool = false,
    kwargs...,
)
    if copy_names
        return MOI.Utilities.automatic_copy_to(
            dest,
            src;
            copy_names = true,
            kwargs...,
        )
    end
    @assert MOI.is_empty(dest)
    _check_input_data(dest, src)
    mapping = MOI.Utilities.IndexMap()
    numcol, colcost = _copy_to_columns(dest, src, mapping)
    collower, colupper = fill(-Inf, numcol), fill(Inf, numcol)
    rowlower, rowupper = Float64[], Float64[]
    I, J, V = Cint[], Cint[], Float64[]
    for S in (
        MOI.GreaterThan{Float64},
        MOI.LessThan{Float64},
        MOI.EqualTo{Float64},
        MOI.Interval{Float64},
    )
        _extract_bound_data(dest, src, mapping, collower, colupper, S)
        _extract_row_data(dest, src, mapping, rowlower, rowupper, I, J, V, S)
    end
    numrow = Cint(length(rowlower))
    A = SparseArrays.sparse(I, J, V, numrow, numcol)
    integrality = fill(kContinuous, numcol)
    for ci in MOI.get(
        src,
        MOI.ListOfConstraintIndices{MOI.SingleVariable,MOI.ZeroOne}(),
    )
        push!(dest.binaries, mapping[ci])
        integrality[_info(dest, ci).column+1] = kInteger
    end
    for ci in MOI.get(
        src,
        MOI.ListOfConstraintIndices{MOI.SingleVariable,MOI.Integer}(),
    )
        integrality[_info(dest, ci).column+1] = kInteger
    end
    Highs_passMip(
        dest,
        numcol,
        numrow,
        length(V),
        kColwise,  # The A matrix is given is column-wise.
        MOI.get(src, MOI.ObjectiveSense()) == MOI.MAX_SENSE ? kMaximize :
        kMinimize,
        0.0,
        colcost,
        collower,
        colupper,
        rowlower,
        rowupper,
        A.colptr .- Cint(1),
        A.rowval .- Cint(1),
        A.nzval,
        integrality,
    )
    return mapping
end
