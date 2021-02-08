import MathOptInterface

const MOI = MathOptInterface
const CleverDicts = MOI.Utilities.CleverDicts

_bounds(s::MOI.GreaterThan{Float64}) = s.lower, Inf
_bounds(s::MOI.LessThan{Float64}) = -Inf, s.upper
_bounds(s::MOI.EqualTo{Float64}) = s.value, s.value
_bounds(s::MOI.Interval{Float64}) = s.lower, s.upper

@enum(
    _RowType,
    _ROW_TYPE_LESSTHAN,
    _ROW_TYPE_GREATERTHAN,
    _ROW_TYPE_INTERVAL,
    _ROW_TYPE_EQUAL_TO,
)

@enum(
    _BoundEnum,
    _BOUND_NONE,
    _BOUND_LESS_THAN,
    _BOUND_GREATER_THAN,
    _BOUND_LESS_AND_GREATER_THAN,
    _BOUND_INTERVAL,
    _BOUND_EQUAL_TO,
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
    column::Cint
    # Storage to keep track of the variable bounds.
    bound::_BoundEnum
    lower::Float64
    upper::Float64
    # We can perform an optimization and only store two strings for the
    # constraint names because, at most, there can be two SingleVariable
    # constraints, e.g., LessThan, GreaterThan.
    lessthan_name::String
    greaterthan_interval_or_equalto_name::String

    function _VariableInfo(
        index::MOI.VariableIndex, column::Cint, bound::_BoundEnum = _BOUND_NONE,
    )
        return new(index, "", column, bound, -Inf, Inf, "", "")
    end
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

function _ConstraintInfo(set::MOI.LessThan{Float64})
    _ConstraintInfo("", 0, _ROW_TYPE_LESSTHAN, -Inf, set.upper)
end

function _ConstraintInfo(set::MOI.GreaterThan{Float64})
    _ConstraintInfo("", 0, _ROW_TYPE_GREATERTHAN, set.lower, Inf)
end

function _ConstraintInfo(set::MOI.Interval{Float64})
    _ConstraintInfo("", 0, _ROW_TYPE_INTERVAL, set.lower, set.upper)
end

function _ConstraintInfo(set::MOI.EqualTo{Float64})
    _ConstraintInfo("", 0, _ROW_TYPE_EQUAL_TO, set.value, set.value)
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

"""
    _Solution

A struct to store the vector solution from HiGHS because it doesn't support
accessing them element-wise.
"""
struct _Solution
    colvalue::Vector{Cdouble}
    coldual::Vector{Cdouble}
    rowvalue::Vector{Cdouble}
    rowdual::Vector{Cdouble}
    _Solution() = new(Cdouble[], Cdouble[], Cdouble[], Cdouble[])
end

mutable struct Optimizer <: MOI.AbstractOptimizer
    # A pointer to the underlying HiGHS optimizer.
    inner::Ptr{Cvoid}

    # Storage for `MOI.Name`.
    name::String

    # A flag to keep track of MOI.Silent, which over-rides the print_level
    # parameter.
    silent::Bool

    # A flag to keep track of MOI.FEASIBILITY_SENSE, since HiGHS only stores
    # MIN_SENSE or MAX_SENSE. This allows us to differentiate between MIN_SENSE
    # and FEASIBILITY_SENSE.
    is_feasibility::Bool

    # HiGHS doesn't support constants in the objective function.
    objective_constant::Float64

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
    optimize_called::Bool
    solution::_Solution

    """
        Optimizer()

    Create a new Optimizer object.
    """
    function Optimizer()
        model = new(
            C_NULL,
            "",
            false,
            false,
            0.0,
            _variable_info_dict(),
            _constraint_info_dict(),
            nothing,
            nothing,
            false,
            _Solution(),
        )
        MOI.empty!(model)
        finalizer(Highs_destroy, model)
        return model
    end
end

Base.cconvert(::Type{Ptr{Cvoid}}, model::Optimizer) = model
Base.unsafe_convert(::Type{Ptr{Cvoid}}, model::Optimizer) = model.inner

### HiGHS has "interesting" error returns. Sometimes it uses non-zero to
### indicate an error, and sometimes it uses true (1) to indicate success and
### false (0) to indicate an error. Annoyingly inconsistent.

function _check_ret(::Optimizer, ret::Cint)
    if ret != Cint(0)
        error(
            "Encountered error code $(ret) in HiGHS. Check the log for details."
        )
    end
    return
end

function _check_bool_ret(::Optimizer, ret::Cint)
    if ret != Cint(1)
        error("Encountered an error in HiGHS. Check the log for details.")
    end
    return
end

function Base.show(io::IO, ::Optimizer)
    print(
        io,
        "A HiGHS model with $(Highs_getNumCols(model)) columns and " *
        "$(Highs_getNumRows(model)) rows."
    )
end

function MOI.empty!(model::Optimizer)
    if model.inner != C_NULL
        Highs_destroy(model)
    end
    model.inner = Highs_create()
    if model.inner == C_NULL
        error("Unable to create an internal model via the C API.")
    end
    model.objective_constant = 0.0
    model.is_feasibility = true
    empty!(model.variable_info)
    empty!(model.affine_constraint_info)
    model.name_to_variable = nothing
    model.name_to_constraint_index = nothing
    model.optimize_called = false
    return
end

function MOI.is_empty(model::Optimizer)
    return Highs_getNumCols(model) == 0 &&
        Highs_getNumRows(model) == 0 &&
        iszero(model.objective_constant) &&
        model.is_feasibility &&
        isempty(model.variable_info) &&
        isempty(model.affine_constraint_info) &&
        model.name_to_variable === nothing &&
        model.name_to_constraint_index === nothing &&
        model.optimize_called == false
end

MOI.get(::Optimizer, ::MOI.SolverName) = "HiGHS"

MOI.get(model::Optimizer, ::MOI.RawSolver) = model

MOI.Utilities.supports_default_copy_to(::Optimizer, ::Bool) = true

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; kws...)
    return MOI.Utilities.automatic_copy_to(dest, src; kws...)
end

function MOI.get(::Optimizer, ::MOI.ListOfVariableAttributesSet)
    return MOI.AbstractVariableAttribute[MOI.VariableName()]
end

function MOI.get(model::Optimizer, ::MOI.ListOfModelAttributesSet)
    attributes = [
        MOI.ObjectiveSense(),
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}()
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

"""
    _OPTIONS

A dictionary mapping the `String` name of HiGHS options to the expected input
type.
"""
const _OPTIONS = Dict{String,DataType}(
    "presolve" => String,
    "solver" => String,
    "parallel" => String,
    "simplex_strategy" => Cint,
    "simplex_iteration_limit" => Cint,
    "highs_min_threads" => Cint,
    "message_level" => Cint,
    "time_limit" => Cdouble,
)

function MOI.supports(::Optimizer, param::MOI.RawParameter)
    return haskey(_OPTIONS, param.name)
end

function MOI.set(model::Optimizer, param::MOI.RawParameter, value::Integer)
    ret = Highs_setHighsIntOptionValue(model, param.name, Cint(value))
    return _check_ret(model, ret)
end

function MOI.set(model::Optimizer, param::MOI.RawParameter, value::Bool)
    ret = Highs_setHighsBoolOptionValue(model, param.name, Cint(value))
    return _check_ret(model, ret)
end

function MOI.set(model::Optimizer, param::MOI.RawParameter, value::AbstractFloat)
    ret = Highs_setHighsDoubleOptionValue(model, param.name, Cdouble(value))
    return _check_ret(model, ret)
end

function MOI.set(model::Optimizer, param::MOI.RawParameter, value::String)
    ret = Highs_setHighsStringOptionValue(model, param.name, value)
    return _check_ret(model, ret)
end

function _get_option(model::Optimizer, option::String, ::Type{Cint})
    value = Ref{Cint}(0)
    Highs_getHighsIntOptionValue(model, option, value)
    return value[]
end

function _get_option(model::Optimizer, option::String, ::Type{Bool})
    value = Ref{Cint}(0)
    Highs_getHighsBoolOptionValue(model, option, value)
    return Bool(value[])
end

function _get_option(model::Optimizer, option::String, ::Type{Cdouble})
    value = Ref{Cdouble}()
    Highs_getHighsDoubleOptionValue(model, option, value)
    return value[]
end

function _get_option(model::Optimizer, option::String, ::Type{String})
    buffer = Vector{Cchar}(undef, 100)
    bufferP = pointer(buffer)
    GC.@preserve buffer begin
        Highs_getHighsStringOptionValue(model, option, bufferP)
        return unsafe_string(bufferP)
    end
end

function MOI.get(o::Optimizer, param::MOI.RawParameter)
    param_type = get(_OPTIONS, param.name, nothing)
    if param_type === nothing
        throw(ArgumentError("Parameter $(param.name) is not supported"))
    end
    return _get_option(o, param.name, param_type)
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

MOI.get(model::Optimizer, ::MOI.Silent) = model.silent

function MOI.set(model::Optimizer, ::MOI.Silent, flag::Bool)
    model.silent = flag
    # TODO(odow): what is the message_level default?
    MOI.set(model, MOI.RawParameter("message_level"), flag ? 1 : 0)
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
    throw(MOI.InvalidIndex(key))
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
        model.variable_info, _VariableInfo(MOI.VariableIndex(0), Cint(0))
    )
    info = _info(model, index)
    # Now, set `.index` and `.column`.
    info.index = index
    info.column = Cint(length(model.variable_info) - 1)
    ret = Highs_addCol(model, 0.0, -Inf, Inf, 0, C_NULL, C_NULL)
    _check_bool_ret(model, ret)
    return index
end

function MOI.is_valid(model::Optimizer, v::MOI.VariableIndex)
    return haskey(model.variable_info, v)
end

function MOI.delete(model::Optimizer, v::MOI.VariableIndex)
    col = column(model, v)
    ret = Highs_deleteColsByRange(model, col, col)
    _check_bool_ret(model, ret)
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
    ::Optimizer, ::MOI.VariableName, ::Type{MOI.VariableIndex}
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
    model::Optimizer, ::MOI.VariableName, v::MOI.VariableIndex, name::String
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
    model::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense
)
    x = sense == MOI.MAX_SENSE ? Cint(-1) : Cint(1)
    ret = Highs_changeObjectiveSense(model, x)
    _check_bool_ret(model, ret)
    if sense == MOI.FEASIBILITY_SENSE
        model.is_feasibility = true
        # TODO(odow): cache the mask.
        n = MOI.get(model, MOI.NumberOfVariables())
        ret = Highs_changeColsCostByMask(
            model, ones(Cint, n), zeros(Cdouble, n)
        )
        _check_bool_ret(model, ret)
        model.objective_constant = 0.0
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
    _check_bool_ret(model, ret)
    return senseP[] == 1 ? MOI.MIN_SENSE : MOI.MAX_SENSE
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
        obj[col + 1] += term.coefficient
    end
    # TODO(odow): cache the mask.
    mask = ones(Cint, num_vars)
    ret = Highs_changeColsCostByMask(model, mask, obj)
    _check_bool_ret(model, ret)
    model.objective_constant = f.constant
    return
end

function MOI.get(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
)
    ncols = MOI.get(model, MOI.NumberOfVariables())
    if ncols == 0
        return MOI.ScalarAffineFunction{Float64}(
            MOI.ScalarAffineTerm{Float64}[], model.objective_constant
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
            MOI.ScalarAffineTerm(costs[info.column + 1], index)
            for (index, info) in model.variable_info
            if !iszero(costs[info.column + 1])
        ],
        model.objective_constant,
    )
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveFunction{F}) where {F}
    obj = MOI.get(
        model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}()
    )
    return convert(F, obj)
end

function MOI.modify(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
    chg::MOI.ScalarConstantChange{Float64}
)
    model.objective_constant = chg.new_constant
    return
end

function MOI.modify(
    model::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
    chg::MOI.ScalarCoefficientChange{Float64}
)
    ret = Highs_changeColCost(
        model, column(model, chg.variable), chg.new_coefficient,
    )
    _check_bool_ret(model, ret)
    return
end

###
### SingleVariable-in-Set constraints.
###

function MOI.supports_constraint(
    ::Optimizer, ::Type{MOI.SingleVariable}, ::Type{<:_SCALAR_SETS}
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
    model::Optimizer, ::MOI.ListOfConstraintIndices{MOI.SingleVariable,S}
) where {S<:_SCALAR_SETS}
    indices = MOI.ConstraintIndex{MOI.SingleVariable, S}[
        MOI.ConstraintIndex{MOI.SingleVariable,S}(key.value)
        for (key, info) in model.variable_info if info.bound in _bound_enums(S)
    ]
    return sort!(indices, by = x -> x.value)
end

function _info(
    model::Optimizer, c::MOI.ConstraintIndex{MOI.SingleVariable,<:Any}
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
    bound::_BoundEnum, ::Type{S}, variable::MOI.VariableIndex,
) where {S<:MOI.AbstractSet}
    if bound == _BOUND_LESS_AND_GREATER_THAN
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Float64},S}(variable))
    elseif bound == _BOUND_GREATER_THAN
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Float64},S}(variable))
    elseif bound == _BOUND_INTERVAL
        throw(MOI.LowerBoundAlreadySet{MOI.Interval{Float64},S}(variable))
    elseif bound == _BOUND_EQUAL_TO
        throw(MOI.LowerBoundAlreadySet{MOI.EqualTo{Float64},S}(variable))
    end
    return
end

function _throw_if_existing_upper(
    bound::_BoundEnum, ::Type{S}, variable::MOI.VariableIndex
) where {S<:MOI.AbstractSet}
    if bound == _BOUND_LESS_AND_GREATER_THAN
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Float64},S}(variable))
    elseif bound == _BOUND_LESS_THAN
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Float64},S}(variable))
    elseif bound == _BOUND_INTERVAL
        throw(MOI.UpperBoundAlreadySet{MOI.Interval{Float64},S}(variable))
    elseif bound == _BOUND_EQUAL_TO
        throw(MOI.UpperBoundAlreadySet{MOI.EqualTo{Float64},S}(variable))
    end
    return
end

function MOI.add_constraint(
    model::Optimizer, f::MOI.SingleVariable, s::S
) where {S <: _SCALAR_SETS}
    info = _info(model, f.variable)
    if S <: MOI.LessThan{Float64}
        _throw_if_existing_upper(info.bound, S, f.variable)
        if info.bound == _BOUND_GREATER_THAN
            info.bound = _BOUND_LESS_AND_GREATER_THAN
        else
            info.bound = _BOUND_LESS_THAN
        end
        info.upper = s.upper
    elseif S <: MOI.GreaterThan{Float64}
        _throw_if_existing_lower(info.bound, S, f.variable)
        if info.bound == _BOUND_LESS_THAN
            info.bound = _BOUND_LESS_AND_GREATER_THAN
        else
            info.bound = _BOUND_GREATER_THAN
        end
        info.lower = s.lower
    elseif S <: MOI.EqualTo{Float64}
        _throw_if_existing_lower(info.bound, S, f.variable)
        _throw_if_existing_upper(info.bound, S, f.variable)
        info.bound = _BOUND_EQUAL_TO
        info.lower = s.value
        info.upper = s.value
    else
        @assert S <: MOI.Interval{Float64}
        _throw_if_existing_lower(info.bound, S, f.variable)
        _throw_if_existing_upper(info.bound, S, f.variable)
        info.bound = _BOUND_INTERVAL
        info.lower = s.lower
        info.upper = s.upper
    end
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
    _check_bool_ret(model, ret)
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
    _check_bool_ret(model, ret)
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
    _check_bool_ret(model, ret)
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
    _check_bool_ret(model, ret)
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
        _check_bool_ret(model, ret)
        info.upper = upper
    elseif S == MOI.GreaterThan{Float64}
        ret = Highs_changeColBounds(model, info.column, lower, info.upper)
        _check_bool_ret(model, ret)
        info.lower = lower
    else
        ret = Highs_changeColBounds(model, info.column, lower, upper)
        _check_bool_ret(model, ret)
        info.lower = lower
        info.upper = upper
    end
    return
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ConstraintName,
    ::Type{<:MOI.ConstraintIndex{MOI.SingleVariable,<:_SCALAR_SETS}}
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
    throw(MOI.InvalidIndex(c))
end

function row(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:_SCALAR_SETS},
)
    return _info(model, c).row
end

function MOI.is_valid(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S}
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
    model::Optimizer, f::MOI.ScalarAffineFunction{Float64}
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
        throw(MOI.ScalarFunctionConstantNotZero{Float64,typeof(f),typeof(s)}(f.constant))
    end
    key = CleverDicts.add_item(
        model.affine_constraint_info, _ConstraintInfo(s)
    )
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
    _check_bool_ret(model, ret)
    return MOI.ConstraintIndex{typeof(f),typeof(s)}(key.value)
end

function MOI.delete(
    model::Optimizer,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:_SCALAR_SETS}
)
    r = row(model, c)
    ret = Highs_deleteRowsByRange(model, r, r)
    _check_bool_ret(model, ret)
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
    _check_bool_ret(model, ret)
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
    _check_bool_ret(model, ret)
    return MOI.ScalarAffineFunction(
        MOI.ScalarAffineTerm{Float64}[
            MOI.ScalarAffineTerm(
                val, model.variable_info[CleverDicts.LinearIndex(col + 1)].index
            )
            for (col, val) in zip(matrix_index, matrix_value) if !iszero(val)
        ],
        0.0,
    )
end

function MOI.supports(
    ::Optimizer,
    ::MOI.ConstraintName,
    ::Type{<:MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:_SCALAR_SETS}}
)
    return true
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:Any},
)
    return _info(model, c).name
end

function MOI.set(
    model::Optimizer,
    ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:Any},
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
    model::Optimizer, C::Type{MOI.ConstraintIndex{F,S}}, name::String
) where {F,S}
    index = MOI.get(model, MOI.ConstraintIndex, name)
    if index isa C
        return index::MOI.ConstraintIndex{F,S}
    end
    return nothing
end

function _rebuild_name_to_constraint_index(model::Optimizer)
    model.name_to_constraint_index = Dict{String,Union{Nothing,MOI.ConstraintIndex}}()
    for (key, info) in model.affine_constraint_info
        if isempty(info.name)
            continue
        end
        S = typeof(_set(info))
        _set_name_to_constraint_index(
            model,
            info.name,
            MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S}(key.value)
        )
    end
    for (key, info) in model.variable_info
        if !isempty(info.lessthan_name)
            _set_name_to_constraint_index(
                model,
                info.lessthan_name,
                MOI.ConstraintIndex{MOI.SingleVariable,MOI.LessThan{Float64}}(key.value)
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
    model::Optimizer, name::String, index::MOI.ConstraintIndex
)
    if haskey(model.name_to_constraint_index, name)
        model.name_to_constraint_index[name] = nothing
    else
        model.name_to_constraint_index[name] = index
    end
    return
end

# TODO(odow): doesn't look like HiGHS supports these.
#
# function MOI.modify(
#     model::Optimizer,
#     c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
#     chg::MOI.ScalarCoefficientChange{Float64}
# ) where {S<:_SCALAR_SETS}
#     return
# end
#
# function MOI.set(
#     model::Optimizer,
#     ::MOI.ConstraintFunction,
#     c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:_SCALAR_SETS},
#     f::MOI.ScalarAffineFunction{Float64}
# )
#     return
# end

###
### Optimize methods.
###

function _store_solution(model::Optimizer)
    x = model.solution
    numCols = Highs_getNumCols(model)
    numRows = Highs_getNumRows(model)
    resize!(x.colvalue, numCols)
    resize!(x.coldual, numCols)
    resize!(x.rowvalue, numRows)
    resize!(x.rowdual, numRows)
    Highs_getSolution(model, x.colvalue, x.coldual, x.rowvalue, x.rowdual)
    return
end

function MOI.optimize!(model::Optimizer)
    ret = Highs_run(model)
    _check_ret(model, ret)
    model.optimize_called = true
    if MOI.get(model, MOI.ResultCount()) == 1
        _store_solution(model)
    end
    return
end

"""
    HighsModelStatus

An enum for the HiGHS simplex status codes.
"""
@enum HighsModelStatus begin
    NOTSET = 0
    LOAD_ERROR
    MODEL_ERROR
    PRESOLVE_ERROR
    SOLVE_ERROR
    POSTSOLVE_ERROR
    MODEL_EMPTY
    PRIMAL_INFEASIBLE
    PRIMAL_UNBOUNDED
    OPTIMAL
    REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND
    REACHED_TIME_LIMIT
    REACHED_ITERATION_LIMIT
    PRIMAL_DUAL_INFEASIBLE
    DUAL_INFEASIBLE
end

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if model.optimize_called == false
        return MOI.OPTIMIZE_NOT_CALLED
    end
    status = HighsModelStatus(Highs_getModelStatus(model.inner, Cint(0)))
    if status == LOAD_ERROR
        return MOI.OTHER_ERROR
    elseif status == MODEL_ERROR
        return MOI.INVALID_MODEL
    elseif status == PRESOLVE_ERROR
        return MOI.OTHER_ERROR
    elseif status == SOLVE_ERROR
        return MOI.OTHER_ERROR
    elseif status == POSTSOLVE_ERROR
        return MOI.OTHER_ERROR
    elseif status == MODEL_EMPTY
        return MOI.INVALID_MODEL
    elseif status == PRIMAL_INFEASIBLE
        return MOI.INFEASIBLE
    elseif status == PRIMAL_UNBOUNDED
        return MOI.DUAL_INFEASIBLE
    elseif status == OPTIMAL
        return MOI.OPTIMAL
    elseif status == REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND
        return MOI.OBJECTIVE_LIMIT
    elseif status == REACHED_TIME_LIMIT
        return MOI.TIME_LIMIT
    elseif status == REACHED_ITERATION_LIMIT
        return MOI.ITERATION_LIMIT
    elseif status == PRIMAL_DUAL_INFEASIBLE
        return MOI.INFEASIBLE
    else
        @assert status == DUAL_INFEASIBLE
        return MOI.DUAL_INFEASIBLE
    end
end

function MOI.get(model::Optimizer, ::MOI.ResultCount)
    status = HighsModelStatus(Highs_getModelStatus(model, Cint(0)))
    return status == OPTIMAL ? 1 : 0
end

function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    status = HighsModelStatus(Highs_getModelStatus(model, Cint(0)))
    return string(status)
end

function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    if attr.N != 1
        return MOI.NO_SOLUTION
    elseif MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
        return MOI.FEASIBLE_POINT
    end
    return MOI.NO_SOLUTION
end

function MOI.get(model::Optimizer, attr::MOI.DualStatus)
    if attr.N != 1
        return MOI.NO_SOLUTION
    elseif MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
        return MOI.FEASIBLE_POINT
    end
    return MOI.NO_SOLUTION
end

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    return Highs_getObjectiveValue(model) + model.objective_constant
end

function MOI.get(model::Optimizer, attr::MOI.DualObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    return MOI.Utilities.get_fallback(model, attr, Float64)
end

function MOI.get(model::Optimizer, ::MOI.SolveTime)
    return Highs_getHighsRunTime(model)
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
    return model.solution.colvalue[column(model, x) + 1]
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.SingleVariable,<:_SCALAR_SETS}
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
    return model.solution.rowvalue[row(model, c) + 1]
end

function _dual_multiplier(model::Optimizer)
    return MOI.get(model, MOI.ObjectiveSense()) == MOI.MAX_SENSE ? -1 : 1
end

"""
    _signed_dual(model::Optimizer, dual::Float64, ::Type{Set})

A heuristic for determining whether the dual of an interval constraint applies
to the lower or upper bound. It can be wrong by at most the solver's tolerance.
"""
function _signed_dual end

function _signed_dual(
    model::Optimizer, dual::Float64, ::Type{MOI.LessThan{Float64}}
)
    return min(_dual_multiplier(model) * dual, 0.0)
end

function _signed_dual(
    model::Optimizer, dual::Float64, ::Type{MOI.GreaterThan{Float64}}
)
    return max(_dual_multiplier(model) * dual, 0.0)
end

function _signed_dual(model::Optimizer, dual::Float64, ::Any)
    return _dual_multiplier(model) * dual
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable,S},
) where {S<:_SCALAR_SETS}
    MOI.check_result_index_bounds(model, attr)
    reduced_cost = model.solution.coldual[column(model, c) + 1]
    return _signed_dual(model, reduced_cost, S)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
) where {S<:_SCALAR_SETS}
    MOI.check_result_index_bounds(model, attr)
    dual = model.solution.rowdual[row(model, c) + 1]
    # TODO(odow): Ask HiGHS why their convention for row duals is the opposite
    # of their convention for column duals. I guess because they transform the
    # constraints into `Ax - Iy = 0`?
    return _signed_dual(model, -dual, S)
end
