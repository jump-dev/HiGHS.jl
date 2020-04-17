import MathOptInterface
const MOI = MathOptInterface

mutable struct Optimizer <: MOI.AbstractOptimizer
    model::ManagedHiGHS
    objective_sense::MOI.OptimizationSense
    variable_map::Dict{MOI.VariableIndex, String}
    Optimizer() = new(ManagedHiGHS(), MOI.FEASIBILITY_SENSE, Dict{MOI.VariableIndex, String}())
end

function MOI.empty!(o::Optimizer)
    reset_model!(o.model)
    o.objective_sense = MOI.FEASIBILITY_SENSE
    empty!(o.variable_map)
    return
end

function MOI.is_empty(o::Optimizer)
    return CWrapper.Highs_getNumRows(o.model.inner) == 0 &&
      CWrapper.Highs_getNumCols(o.model.inner) == 0 &&
      o.objective_sense == MOI.FEASIBILITY_SENSE
end

function MOI.optimize!(o::Optimizer)
    CWrapper.Highs_run(o.model.inner)
    return
end

function MOI.add_variable(o::Optimizer)
    _ = CWrapper.Highs_addCol(o.model.inner, 0.0, -Inf, Inf, Cint(0), Cint[], Cint[])
    col_idx = MOI.get(o, MOI.NumberOfVariables()) - 1
    return MOI.VariableIndex(col_idx)
end

function MOI.add_constrained_variable(o::Optimizer, set::S) where {S <: MOI.Interval}
    _ = CWrapper.Highs_addCol(o.model.inner, 0.0, Cdouble(set.lower), Cdouble(set.upper), Cint(0), Cint[], Cint[])
    col_idx = MOI.get(o, MOI.NumberOfVariables()) - 1
    return (MOI.VariableIndex(col_idx), MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}(col_idx))
end

MOI.supports_constraint(::Optimizer, ::MOI.SingleVariable, ::MOI.Interval) = true

function MOI.add_constraint(o::Optimizer, sg::MOI.SingleVariable, set::MOI.Interval)
    var_idx = Cint(sg.variable.value)
    _ = CWrapper.Highs_changeColBounds(o.model.inner, var_idx, Cdouble(set.lower), Cdouble(set.upper))
    return
end

function MOI.get(o::Optimizer, ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64}})
    nrows = CWrapper.Highs_getNumRows(o.model.inner)
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
    row_idx = CWrapper.Highs_addRow(o.model.inner, lower, upper, Cint(number_nonzeros), col_indices_ptr, coefficients_ptr)
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(row_idx)
end

MOI.supports(::Optimizer, ::MOI.SolverName) = true
MOI.get(::Optimizer, ::MOI.SolverName) = "HiGHS"

MOI.supports(::Optimizer, param::MOI.RawParameter) = true

# setting HiGHS options
function MOI.set(o::Optimizer, param::MOI.RawParameter, value)
    CWrapper.set_option(o.model.inner, param.name, value)
    return nothing
end

function MOI.get(o::Optimizer, param::MOI.RawParameter)
    return if param.name in options_string
        CWrapper.get_option(o.model.inner, param.name, String)
    elseif param.name in options_int
        CWrapper.get_option(o.model.inner, param.name, Int)
    elseif param.name in options_double
        CWrapper.get_option(o.model.inner, param.name, Float64)        
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
    status = CWrapper.Highs_getModelStatus(o.model.inner)
    status == 9 ? 1 : 0
end

function MOI.get(o::Optimizer, ::MOI.ObjectiveSense)
    return o.objective_sense
end

function MOI.set(o::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    sense_code = sense == MOI.MAX_SENSE ? Cint(-1) : Cint(1)
    _ = CWrapper.Highs_changeObjectiveSense(o.model.inner, sense_code)
    o.objective_sense = sense
    return nothing
end

function MOI.get(o::Optimizer, ::MOI.NumberOfVariables)
    return Int(CWrapper.Highs_getNumCols(o.model.inner))
end

function MOI.get(o::Optimizer, ::MOI.ListOfVariableIndices)
    ncols = MOI.get(o, MOI.NumberOfVariables())
    return [MOI.VariableIndex(j) for j in 0:(ncols-1)]
end

MOI.get(::Optimizer, ::MOI.ObjectiveFunctionType) = MOI.ScalarAffineFunction{Float64}

MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}) = true

function MOI.set(o::Optimizer, ::MOI.ObjectiveFunction{F}, func::F) where {F <: MOI.ScalarAffineFunction{Float64}}
    # TODO treat constant in objective
    total_ncols = MOI.get(o, MOI.NumberOfVariables())
    coefficients = zeros(Cdouble, total_ncols)
    for term in func.terms
        j = term.variable_index.value
        coefficients[j+1] = Cdouble(term.coefficient)
    end
    coefficients_ptr = pointer(coefficients)
    mask = pointer(ones(Cint, total_ncols))
    CWrapper.Highs_changeColsCostByMask(o.model.inner, mask, coefficients_ptr)
    return
end

function MOI.get(o::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}})
    # TODO treat constant in objective
    ncols = MOI.get(o, MOI.NumberOfVariables())
    num_cols = Ref{Cint}(0)
    costs = Vector{Cdouble}(undef, ncols)
    _ = CWrapper.Highs_getColsByRange(
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
    return MOI.ScalarAffineFunction{Float64}(terms, 0.0)
end

function MOI.get(o::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(o, attr)
    value = Ref{Cdouble}()
    CWrapper.Highs_getHighsDoubleInfoValue(o.model.inner, "objective_function_value", value);
    return value[]
end

function MOI.get(o::Optimizer, ::MOI.SimplexIterations)
    simplex_iteration_count = Ref{Cint}(0)
    CWrapper.Highs_getHighsIntInfoValue(o.model.inner, "simplex_iteration_count", simplex_iteration_count)
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
    CWrapper.Highs_getHighsDoubleOptionValue(o.model.inner, "time_limit", value)
    return value[]
end

function MOI.set(o::Optimizer, ::MOI.TimeLimitSec, value)
    CWrapper.set_option(o.model.inner, "time_limit", Float64(value))
    return
end
