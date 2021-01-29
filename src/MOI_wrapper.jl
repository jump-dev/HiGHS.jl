
mutable struct ManagedHiGHS
    inner::Ptr{Cvoid}

    function ManagedHiGHS()
        mghs = new(
            Highs_create()
        )
        # register the memory cleanup function
        finalizer(mghs) do m
            success = free_highs(m)
            if !success
                @warn "Memory free failure, possible leak."
            end
        end
    end
end

"""
Release references and free memory and return a boolean
indicating the success of the operation.
"""
function free_highs(mhgs::ManagedHiGHS)
    # Avoid double-free (ManagedHiGHS will set the pointers to NULL).
    if mhgs.inner == C_NULL
        return false
    end
    # only mhgs.inner is GC-protected during ccall!
    GC.@preserve mhgs begin
        Highs_destroy(mhgs.inner)
    end
    mhgs.inner = C_NULL
    return true
end

"""
    reset_model!(mhgs::ManagedHiGHS)

Deletes the inner HiGHS model and recreates one.
"""
function reset_model!(mhgs::ManagedHiGHS)
    # Avoid double-free (ManagedHiGHS will set the pointers to NULL).
    if mhgs.inner == C_NULL
        return false
    end
    # only mhgs.inner is GC-protected during ccall!
    GC.@preserve mhgs begin
        Highs_destroy(mhgs.inner)
    end
    mhgs.inner = Highs_create()
    return true
end


import MathOptInterface
const MOI = MathOptInterface

mutable struct Optimizer <: MOI.AbstractOptimizer
    model::ManagedHiGHS
    objective_sense::MOI.OptimizationSense
    variable_map::Dict{MOI.VariableIndex, String}
    objective_constant::Float64
    Optimizer() = new(ManagedHiGHS(), MOI.FEASIBILITY_SENSE, Dict{MOI.VariableIndex, String}(), 0.0)
end

function MOI.empty!(o::Optimizer)
    reset_model!(o.model)
    o.objective_sense = MOI.FEASIBILITY_SENSE
    empty!(o.variable_map)
    o.objective_constant = 0.0
    return
end

function MOI.is_empty(o::Optimizer)
    return Highs_getNumRows(o.model.inner) == 0 &&
      Highs_getNumCols(o.model.inner) == 0 &&
      o.objective_sense == MOI.FEASIBILITY_SENSE &&
      o.objective_constant == 0.0
end

function MOI.optimize!(o::Optimizer)
    if o.objective_sense == MOI.FEASIBILITY_SENSE
        obj_func = MOI.get(o, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
        if any(!=(0.0), obj_func.terms)
            error("Feasibility sense with non-constant objective function, set the sense to min/max, the objective to 0 or reset the sense to feasibility to erase the objective")
        end
    end
    Highs_run(o.model.inner)
    return
end

function MOI.add_variable(o::Optimizer)
    _ = Highs_addCol(o.model.inner, 0.0, -Inf, Inf, Cint(0), Cint[], Cint[])
    col_idx = MOI.get(o, MOI.NumberOfVariables()) - 1
    return MOI.VariableIndex(col_idx)
end

function MOI.add_constrained_variable(o::Optimizer, set::S) where {S <: MOI.Interval}
    _ = Highs_addCol(o.model.inner, 0.0, Cdouble(set.lower), Cdouble(set.upper), Cint(0), Cint[], Cint[])
    col_idx = MOI.get(o, MOI.NumberOfVariables()) - 1
    return (MOI.VariableIndex(col_idx), MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}(col_idx))
end

MOI.supports_constraint(::Optimizer, ::MOI.SingleVariable, ::MOI.Interval) = true

function MOI.add_constraint(o::Optimizer, sg::MOI.SingleVariable, set::MOI.Interval)
    var_idx = Cint(sg.variable.value)
    _ = Highs_changeColBounds(o.model.inner, var_idx, Cdouble(set.lower), Cdouble(set.upper))
    return
end

function MOI.get(o::Optimizer, ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64}})
    nrows = Highs_getNumRows(o.model.inner)
    return Int(nrows)
end

MOI.supports_constraint(::Optimizer, ::MOI.ScalarAffineFunction, ::MOI.Interval) = true

function MOI.add_constraint(o::Optimizer, func::MOI.ScalarAffineFunction, set::MOI.Interval)
    number_nonzeros = length(func.terms)
    coefficients = Vector{Cdouble}(undef, number_nonzeros)
    col_indices = Vector{Cint}(undef, number_nonzeros)
    for j in Base.OneTo(number_nonzeros)
        term = func.terms[j]
        coefficients[j] = term.coefficient
        col_indices[j] = Cint(term.variable_index.value)
    end
    coefficients_ptr = pointer(coefficients)
    col_indices_ptr = pointer(col_indices)
    lower = convert(Cdouble, set.lower - func.constant)
    upper = convert(Cdouble, set.upper - func.constant)
    row_idx = Highs_addRow(o.model.inner, lower, upper, Cint(number_nonzeros), col_indices_ptr, coefficients_ptr)
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(row_idx)
end

MOI.supports(::Optimizer, ::MOI.SolverName) = true
MOI.get(::Optimizer, ::MOI.SolverName) = "HiGHS"

MOI.supports(::Optimizer, param::MOI.RawParameter) = true

# setting HiGHS options
function MOI.set(o::Optimizer, param::MOI.RawParameter, value)
    set_option(o.model.inner, param.name, value)
    return nothing
end

function MOI.get(o::Optimizer, param::MOI.RawParameter)
    return if param.name in options_string
        get_option(o.model.inner, param.name, String)
    elseif param.name in options_int
        get_option(o.model.inner, param.name, Int)
    elseif param.name in options_double
        get_option(o.model.inner, param.name, Float64)
    else
        throw(ArgumentError("Parameter $(param.name) is not supported"))
    end
end

const SUPPORTED_MODEL_ATTRIBUTES = Union{
    MOI.ObjectiveSense,
    MOI.NumberOfVariables,
    MOI.ListOfVariableIndices,
    MOI.ListOfConstraintIndices,
    MOI.NumberOfConstraints, # TODO single variables
    MOI.ListOfConstraints, # TODO single variables
    MOI.ObjectiveFunctionType,
    MOI.ObjectiveValue,
    # MOI.DualObjectiveValue, # TODO
    # MOI.SolveTime,  # TODO
    MOI.SimplexIterations,
    MOI.BarrierIterations,
    MOI.RawSolver,
    # MOI.RawStatusString,  # TODO
    MOI.ResultCount,
    # MOI.Silent, # TODO
    # MOI.TerminationStatus,
    # MOI.PrimalStatus,
    # MOI.DualStatus
}

MOI.supports(::Optimizer, ::SUPPORTED_MODEL_ATTRIBUTES) = true

MOI.supports(::Optimizer, ::MOI.VariableName, ::Type{MOI.VariableIndex}) = true

MOI.get(o::Optimizer, ::MOI.RawSolver) = o.model

function MOI.get(o::Optimizer, ::MOI.ResultCount)
    status = Highs_getModelStatus(o.model.inner)
    status == 9 ? 1 : 0
end

function MOI.get(o::Optimizer, ::MOI.ObjectiveSense)
    return o.objective_sense
end

function MOI.set(o::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    sense_code = sense == MOI.MAX_SENSE ? Cint(-1) : Cint(1)
    _ = Highs_changeObjectiveSense(o.model.inner, sense_code)
    o.objective_sense = sense
    # if feasibility sense set, erase the function
    if sense == MOI.FEASIBILITY_SENSE
        MOI.set(o, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarAffineFunction{Float64}([], 0))
    end
    return nothing
end

function MOI.get(o::Optimizer, ::MOI.NumberOfVariables)
    return Int(Highs_getNumCols(o.model.inner))
end

function MOI.get(o::Optimizer, ::MOI.ListOfVariableIndices)
    ncols = MOI.get(o, MOI.NumberOfVariables())
    return [MOI.VariableIndex(j) for j in 0:(ncols-1)]
end

MOI.get(::Optimizer, ::MOI.ObjectiveFunctionType) = MOI.ScalarAffineFunction{Float64}

MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}) = true

function MOI.set(o::Optimizer, ::MOI.ObjectiveFunction{F}, func::F) where {F <: MOI.ScalarAffineFunction{Float64}}
    total_ncols = MOI.get(o, MOI.NumberOfVariables())
    coefficients = zeros(Cdouble, total_ncols)
    for term in func.terms
        j = term.variable_index.value
        coefficients[j+1] = Cdouble(term.coefficient)
    end
    coefficients_ptr = pointer(coefficients)
    mask = pointer(ones(Cint, total_ncols))
    Highs_changeColsCostByMask(o.model.inner, mask, coefficients_ptr)
    o.objective_constant = MOI.constant(func)
    return
end

function MOI.get(o::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}})
    ncols = MOI.get(o, MOI.NumberOfVariables())
    num_cols = Ref{Cint}(0)
    costs = Vector{Cdouble}(undef, ncols)
    _ = Highs_getColsByRange(
        o.model.inner,
        Cint(0), Cint(ncols-1), # column range
        num_cols, costs,
        C_NULL, C_NULL, # lower, upper
        Ref{Cint}(0), C_NULL, C_NULL, C_NULL # coefficients
    )
    num_cols[] == ncols || error("Unexpected number of columns, inconsistent HiGHS state")
    terms = Vector{MOI.ScalarAffineTerm{Float64}}()
    for (j, cost) in enumerate(costs)
        if cost != 0.0
            var_idx = MOI.VariableIndex(j-1)
            push!(terms, MOI.ScalarAffineTerm(cost, var_idx))
        end
    end
    return MOI.ScalarAffineFunction{Float64}(terms, o.objective_constant)
end

function MOI.get(o::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(o, attr)
    value = Ref{Cdouble}()
    Highs_getHighsDoubleInfoValue(o.model.inner, "objective_function_value", value);
    return o.objective_constant + value[]
end

function MOI.get(o::Optimizer, ::MOI.SimplexIterations)
    simplex_iteration_count = Ref{Cint}(0)
    Highs_getHighsIntInfoValue(o.model.inner, "simplex_iteration_count", simplex_iteration_count)
    return Int(simplex_iteration_count[])
end

function MOI.get(o::Optimizer, ::MOI.BarrierIterations)
    return 0
end

function MOI.get(o::Optimizer, ::MOI.VariableName, v::MOI.VariableIndex)
    return get(o.variable_map, v, "")
end

function MOI.set(o::Optimizer, ::MOI.VariableName, v::MOI.VariableIndex, name::String)
    o.variable_map[v] = name
    return
end

function MOI.get(o::Optimizer, ::Type{MOI.VariableIndex}, name::String)
    for (vi, vname) in o.variable_map
        if vname == name
            return vi
        end
    end
    return nothing
end

function MOI.get(o::Optimizer, ::MOI.TimeLimitSec)
    value = Ref{Cdouble}()
    Highs_getHighsDoubleOptionValue(o.model.inner, "time_limit", value)
    return value[]
end

function MOI.set(o::Optimizer, ::MOI.TimeLimitSec, value)
    set_option(o.model.inner, "time_limit", Float64(value))
    return
end

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

# add all options methods here

const options_string = ["presolve", "solver", "parallel"]
const options_int = ["simplex_strategy", "simplex_iteration_limit", "highs_min_threads", "message_level"]
const options_double = ["time_limit"]

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
