import MathOptInterface

const MOI = MathOptInterface
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
        index::MOI.VariableIndex,
        column::Cint,
        bound::_BoundEnum = _BOUND_NONE,
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

"""
    _Solution

A struct to store the vector solution from HiGHS because it doesn't support
accessing them element-wise.
"""
mutable struct _Solution
    optimize_called::Bool
    colvalue::Vector{Cdouble}
    coldual::Vector{Cdouble}
    colstatus::Vector{Cint}
    rowvalue::Vector{Cdouble}
    rowdual::Vector{Cdouble}
    rowstatus::Vector{Cint}
    has_solution::Bool
    has_primal_ray::Bool
    has_dual_ray::Bool
    function _Solution()
        return new(
            false,
            Cdouble[],
            Cdouble[],
            Cint[],
            Cdouble[],
            Cdouble[],
            Cint[],
            false,
            false,
            false,
        )
    end
end

function Base.empty!(x::_Solution)
    x.optimize_called = false
    x.has_solution = false
    x.has_dual_ray = false
    x.has_primal_ray = false
    return x
end

Base.isempty(x::_Solution) = !x.optimize_called
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
    solution::_Solution

    """
        Optimizer()

    Create a new Optimizer object.
    """
    function Optimizer()
        ptr = Highs_create()
        if ptr == C_NULL
            error("Unable to create an internal model via the C API.")
        end
        model = new(
            ptr,
            "",
            false,
            true,
            0.0,
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
    if ret == Cint(1)
        return  # Cint(1) is 'true'. Nothing went wrong.
    else
        # These return codes should only ever be 0 or 1.
        @assert ret == Cint(0)
        error("Encountered an error in HiGHS. Check the log for details.")
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
    ret = Highs_clearModel(model)
    if ret == 2
        error(
            "Encountered an error in HiGHS: HighsStatus::Error. " *
            "Check the log for details",
        )
    end
    model.objective_constant = 0.0
    model.is_feasibility = true
    empty!(model.variable_info)
    empty!(model.affine_constraint_info)
    model.name_to_variable = nothing
    model.name_to_constraint_index = nothing
    empty!(model.solution)
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
           isempty(model.solution)
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
    return Highs_getHighsOptionType(model, param.name, typeP) == 0
end

function _check_option_status(ret::Cint)
    if ret != 0
        error("Encountered an error in HiGHS: Check the log for details.")
    end
    return
end

function _set_option(model::Optimizer, option::String, value::Bool)
    return Highs_setHighsBoolOptionValue(model, option, Cint(value))
end

function _set_option(model::Optimizer, option::String, value::Integer)
    return Highs_setHighsIntOptionValue(model, option, Cint(value))
end

function _set_option(model::Optimizer, option::String, value::AbstractFloat)
    return Highs_setHighsDoubleOptionValue(model, option, Cdouble(value))
end

function _set_option(model::Optimizer, option::String, value::String)
    return Highs_setHighsStringOptionValue(model, option, value)
end

function MOI.set(model::Optimizer, param::MOI.RawParameter, value)
    if !MOI.supports(model, param)
        throw(MOI.UnsupportedAttribute(param))
    end
    ret = _set_option(model, param.name, value)
    return _check_option_status(ret)
end

### MOI.get

function _get_bool_option(model::Optimizer, option::String)
    value = Ref{Cint}(0)
    ret = Highs_getHighsBoolOptionValue(model, option, value)
    _check_option_status(ret)
    return Bool(value[])
end

function _get_int_option(model::Optimizer, option::String)
    value = Ref{Cint}(0)
    ret = Highs_getHighsIntOptionValue(model, option, value)
    _check_option_status(ret)
    return value[]
end

function _get_double_option(model::Optimizer, option::String)
    value = Ref{Cdouble}()
    ret = Highs_getHighsDoubleOptionValue(model, option, value)
    _check_option_status(ret)
    return value[]
end

function _get_string_option(model::Optimizer, option::String)
    buffer = Vector{Cchar}(undef, 1024)
    bufferP = pointer(buffer)
    GC.@preserve buffer begin
        ret = Highs_getHighsStringOptionValue(model, option, bufferP)
        _check_option_status(ret)
        return unsafe_string(bufferP)
    end
end

function MOI.get(model::Optimizer, param::MOI.RawParameter)
    if !(param.name isa String)
        throw(MOI.UnsupportedAttribute(param))
    end
    typeP = Ref{Cint}()
    ret = Highs_getHighsOptionType(model, param.name, typeP)
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

MOI.get(model::Optimizer, ::MOI.Silent) = model.silent

function MOI.set(model::Optimizer, ::MOI.Silent, flag::Bool)
    if flag
        Highs_runQuiet(model)
    elseif model.silent && !flag
        @warn("Unable to restore printing. Sorry.")
    end
    model.silent = flag
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
    x = sense == MOI.MAX_SENSE ? Cint(-1) : Cint(1)
    ret = Highs_changeObjectiveSense(model, x)
    _check_ret(ret)
    if sense == MOI.FEASIBILITY_SENSE
        model.is_feasibility = true
        # TODO(odow): cache the mask.
        n = MOI.get(model, MOI.NumberOfVariables())
        ret =
            Highs_changeColsCostByMask(model, ones(Cint, n), zeros(Cdouble, n))
        _check_ret(ret)
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
    _check_ret(ret)
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
        obj[col+1] += term.coefficient
    end
    # TODO(odow): cache the mask.
    mask = ones(Cint, num_vars)
    ret = Highs_changeColsCostByMask(model, mask, obj)
    _check_ret(ret)
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
            MOI.ScalarAffineTerm{Float64}[],
            model.objective_constant,
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
        model.objective_constant,
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
    model.objective_constant = chg.new_constant
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
    bound::_BoundEnum,
    ::Type{S},
    variable::MOI.VariableIndex,
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
    bound::_BoundEnum,
    ::Type{S},
    variable::MOI.VariableIndex,
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
    model::Optimizer,
    f::MOI.SingleVariable,
    s::S,
) where {S<:_SCALAR_SETS}
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
    ::Type{<:MOI.ConstraintIndex{MOI.SingleVariable,<:_SCALAR_SETS}},
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
    chg::MOI.ScalarCoefficientChange{Float64}
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
    f::F
) where {F<:MOI.ScalarAffineFunction{Float64},S<:_SCALAR_SETS}
    if !iszero(f.constant)
        throw(MOI.ScalarFunctionConstantNotZero{Float64,F,S}(f.constant))
    end
    # This loop is a bit slow because we query the existing function first, and
    # then set each coefficient one-by-one. It's probably okay until someone
    # complainsm but it will require some upstream changes to introduce a faster
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

function _store_solution(model::Optimizer)
    x = model.solution
    x.optimize_called = true
    x.has_solution = false
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
    # Load the solution if optimal.
    if Highs_getModelStatus(model.inner, Cint(0)) == 9
        Highs_getSolution(model, x.colvalue, x.coldual, x.rowvalue, x.rowdual)
        Highs_getBasis(model, x.colstatus, x.rowstatus)
        x.has_solution = true
        return
    end
    # Check for a certificate of primal infeasibility.
    has_ray = Ref{Cint}(0)
    ret = Highs_getDualRay(model, has_ray, x.rowdual)
    x.has_dual_ray = has_ray[] == 1
    @assert ret == 0  # getDualRay only returns HighsStatus::OK
    if x.has_dual_ray
        return
    end
    # Check for a certificate of dual infeasibility.
    ret = Highs_getPrimalRay(model, has_ray, x.colvalue)
    x.has_primal_ray = has_ray[] == 1
    @assert ret == 0  # getPrimalRay only returns HighsStatus::OK
    return
end

function MOI.optimize!(model::Optimizer)
    ret = Highs_run(model)
    # `ret` is an `HighsStatus` enum: {OK = 0, Warning, Error}.
    if ret == 2
        error(
            "Encountered an error in HiGHS: HighsStatus::Error. " *
            "Check the log for details",
        )
    end
    _store_solution(model)
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
    if !model.solution.optimize_called
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
    x = model.solution
    if x.has_solution || x.has_dual_ray || x.has_primal_ray
        return 1
    else
        return 0
    end
end

function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    status = HighsModelStatus(Highs_getModelStatus(model, Cint(0)))
    return string(status)
end

function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    if attr.N != 1
        return MOI.NO_SOLUTION
    elseif model.solution.has_solution
        return MOI.FEASIBLE_POINT
    elseif model.solution.has_primal_ray
        return MOI.INFEASIBILITY_CERTIFICATE
    end
    return MOI.NO_SOLUTION
end

function MOI.get(model::Optimizer, attr::MOI.DualStatus)
    if attr.N != 1
        return MOI.NO_SOLUTION
    elseif model.solution.has_solution
        return MOI.FEASIBLE_POINT
    elseif model.solution.has_dual_ray
        return MOI.INFEASIBILITY_CERTIFICATE
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
    for i = 1:num_nz[]
        dual += -model.solution.rowdual[matrix_index[i] + 1] * matrix_value[i]
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

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable,S},
) where {S<:_SCALAR_SETS}
    MOI.check_result_index_bounds(model, attr)
    col = column(model, c)
    if model.solution.has_dual_ray[] == 1
        return _signed_dual(_farkas_variable_dual(model, col), S)
    else
        reduced_cost = _sense_corrector(model) * model.solution.coldual[col+1]
        return _signed_dual(reduced_cost, S)
    end
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
) where {S<:_SCALAR_SETS}
    MOI.check_result_index_bounds(model, attr)
    dual = model.solution.rowdual[row(model, c)+1]
    if model.solution.has_dual_ray[] == 1
        return _signed_dual(dual, S)
    else
        return _signed_dual(-_sense_corrector(model) * dual, S)
    end
end

###
### MOI.ConstraintBasisStatus
###

@enum(
    HighsBasisStatus,
    LOWER = 0,
    BASIC,
    UPPER,
    ZERO,
    NONBASIC,
    SUPER,
)

# HiGHS only reports a single basis status for each ranged constraint. Therefore
# if the status is LOWER or UPPER, we must distinguish whether or not it is
# refering to the given constraint.
function _nonbasic_status(code::HighsBasisStatus, ::Type{<:MOI.Interval})
    return code == LOWER ? MOI.NONBASIC_AT_LOWER : MOI.NONBASIC_AT_UPPER
end
function _nonbasic_status(code::HighsBasisStatus, ::Type{<:MOI.LessThan})
    return code == LOWER ? MOI.BASIC : MOI.NONBASIC
end
function _nonbasic_status(code::HighsBasisStatus, ::Type{<:MOI.GreaterThan})
    return code == LOWER ? MOI.NONBASIC : MOI.BASIC
end
_nonbasic_status(::HighsBasisStatus, ::Type{<:MOI.EqualTo}) = MOI.NONBASIC

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintBasisStatus,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
) where {S<:_SCALAR_SETS}
    MOI.check_result_index_bounds(model, attr)
    stat = HighsBasisStatus(model.solution.rowstatus[row(model, c) + 1])
    if stat == LOWER
        return _nonbasic_status(stat, S)
    elseif stat == BASIC
        return MOI.BASIC
    elseif stat == UPPER
        return _nonbasic_status(stat, S)
    elseif stat == ZERO || stat == NONBASIC
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
    stat = HighsBasisStatus(model.solution.colstatus[column(model, c) + 1])
    if stat == LOWER
        return _nonbasic_status(stat, S)
    elseif stat == BASIC
        return MOI.BASIC
    elseif stat == UPPER
        return _nonbasic_status(stat, S)
    elseif stat == ZERO || stat == NONBASIC
        return MOI.NONBASIC
    else
        @assert stat == SUPER
        return MOI.SUPER_BASIC
    end
end
