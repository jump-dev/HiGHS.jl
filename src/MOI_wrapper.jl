import MathOptInterface
const MOI = MathOptInterface

mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::Ptr{Cvoid}
    objective_sense::MOI.OptimizationSense
    variable_map::Dict{MOI.VariableIndex, String}
    objective_constant::Float64

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
            MOI.FEASIBILITY_SENSE,
            Dict{MOI.VariableIndex, String}(),
            0.0,
        )
        finalizer(Highs_destroy, model)
        return model
    end
end

Base.cconvert(::Type{Ptr{Cvoid}}, model::Optimizer) = model
Base.unsafe_convert(::Type{Ptr{Cvoid}}, model::Optimizer) = model.inner

function _check_ret(model::Optimizer, ret::Cint)
    # TODO: errors for invalid return codes!
    return nothing
end

function MOI.empty!(o::Optimizer)
    Highs_destroy(o)
    o.inner = Highs_create()
    o.objective_sense = MOI.FEASIBILITY_SENSE
    empty!(o.variable_map)
    o.objective_constant = 0.0
    return
end

function MOI.is_empty(o::Optimizer)
    return Highs_getNumRows(o) == 0 &&
      Highs_getNumCols(o) == 0 &&
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
    Highs_run(o)
    return
end

function MOI.add_variable(o::Optimizer)
    _ = Highs_addCol(o, 0.0, -Inf, Inf, Cint(0), Cint[], Cint[])
    col_idx = MOI.get(o, MOI.NumberOfVariables()) - 1
    return MOI.VariableIndex(col_idx)
end

function MOI.add_constrained_variable(o::Optimizer, set::S) where {S <: MOI.Interval}
    _ = Highs_addCol(o, 0.0, Cdouble(set.lower), Cdouble(set.upper), Cint(0), Cint[], Cint[])
    col_idx = MOI.get(o, MOI.NumberOfVariables()) - 1
    return (MOI.VariableIndex(col_idx), MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}(col_idx))
end

MOI.supports_constraint(::Optimizer, ::MOI.SingleVariable, ::MOI.Interval) = true

function MOI.add_constraint(o::Optimizer, sg::MOI.SingleVariable, set::MOI.Interval)
    var_idx = Cint(sg.variable.value)
    _ = Highs_changeColBounds(o, var_idx, Cdouble(set.lower), Cdouble(set.upper))
    return
end

function MOI.get(o::Optimizer, ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64}})
    nrows = Highs_getNumRows(o)
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
    row_idx = Highs_addRow(o, lower, upper, Cint(number_nonzeros), col_indices_ptr, coefficients_ptr)
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(row_idx)
end

MOI.supports(::Optimizer, ::MOI.SolverName) = true
MOI.get(::Optimizer, ::MOI.SolverName) = "HiGHS"

MOI.supports(::Optimizer, param::MOI.RawParameter) = true

# setting HiGHS options
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

function MOI.get(o::Optimizer, param::MOI.RawParameter)
    param_type = get(_OPTIONS, param.name, nothing)
    if param_type === nothing
        throw(ArgumentError("Parameter $(param.name) is not supported"))
    end
    return _get_option(o, param.name, param_type)
end

const _SUPPORTED_MODEL_ATTRIBUTES = Union{
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
    MOI.TerminationStatus,
    # MOI.PrimalStatus,
    # MOI.DualStatus
}
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
MOI.supports(::Optimizer, ::_SUPPORTED_MODEL_ATTRIBUTES) = true

MOI.supports(::Optimizer, ::MOI.VariableName, ::Type{MOI.VariableIndex}) = true

MOI.get(o::Optimizer, ::MOI.RawSolver) = o

function MOI.get(o::Optimizer, ::MOI.ResultCount)
    status = HighsModelStatus(Highs_getModelStatus(o, Cint(0)))
    return status == OPTIMAL ? 1 : 0
end
function MOI.get(o::Optimizer, ::MOI.ObjectiveSense)
    return o.objective_sense
end
function MOI.get(o::Optimizer, ::MOI.TerminationStatus)
    status = HighsModelStatus(Highs_getModelStatus(o.inner, Cint(0)))
    if status == OPTIMAL
        return MOI.OPTIMAL
    elseif status == PRIMAL_INFEASIBLE
        return MOI.INFEASIBLE
    elseif status == DUAL_INFEASIBLE
        return MOI.DUAL_INFEASIBLE
    elseif status == REACHED_ITERATION_LIMIT
        return MOI.ITERATION_LIMIT
    elseif status == REACHED_TIME_LIMIT
        return MOI.TIME_LIMIT
    elseif status == MODEL_ERROR
        return MOI.INVALID_MODEL
    else
        return MOI.OTHER_ERROR
    end
end

function MOI.set(o::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    sense_code = sense == MOI.MAX_SENSE ? Cint(-1) : Cint(1)
    _ = Highs_changeObjectiveSense(o, sense_code)
    o.objective_sense = sense
    # if feasibility sense set, erase the function
    if sense == MOI.FEASIBILITY_SENSE
        MOI.set(o, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarAffineFunction{Float64}([], 0))
    end
    return nothing
end

function MOI.get(o::Optimizer, ::MOI.NumberOfVariables)
    return Int(Highs_getNumCols(o))
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
    Highs_changeColsCostByMask(o, mask, coefficients_ptr)
    o.objective_constant = MOI.constant(func)
    return
end

function MOI.get(o::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}})
    ncols = MOI.get(o, MOI.NumberOfVariables())
    num_cols = Ref{Cint}(0)
    costs = Vector{Cdouble}(undef, ncols)
    if ncols > 0
        _ = Highs_getColsByRange(
            o,
            Cint(0), Cint(ncols-1), # column range
            num_cols, costs,
            C_NULL, C_NULL, # lower, upper
            Ref{Cint}(0), C_NULL, C_NULL, C_NULL # coefficients
        )
        num_cols[] == ncols || error("Unexpected number of columns, inconsistent HiGHS state")
    end
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
    Highs_getHighsDoubleInfoValue(o, "objective_function_value", value);
    return o.objective_constant + value[]
end

function MOI.get(o::Optimizer, ::MOI.SimplexIterations)
    simplex_iteration_count = Ref{Cint}(0)
    Highs_getHighsIntInfoValue(o, "simplex_iteration_count", simplex_iteration_count)
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
    return MOI.get(o, MOI.RawParameter("time_limit"))
end

function MOI.set(o::Optimizer, ::MOI.TimeLimitSec, value::Real)
    return MOI.set(o, MOI.RawParameter("time_limit"), Cdouble(value))
end

function MOI.set(o::Optimizer, ::MOI.TimeLimitSec, ::Nothing)
    # TODO(odow): handle default time limit?
    return
end
