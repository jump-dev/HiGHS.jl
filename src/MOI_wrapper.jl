import MathOptInterface
const MOI = MathOptInterface

mutable struct Optimizer <: MOI.AbstractOptimizer
    model::ManagedHiGHS
    Optimizer() = new(ManagedHiGHS())
end

function MOI.empty!(o::Optimizer)
    CWrapper.reset_model!(o.model)
end

function MOI.add_variable(o::Optimizer)
    idx = CWrapper.Highs_addCol(o.model.inner, 0.0, -Inf, Inf, Cint(0), Cint[], Cint[])
    return MOI.VariableIndex(Int(idx))
end

function MOI.add_constrained_variable(o::Optimizer, set::S) where {S <: MOI.Interval}
    idx = CWrapper.Highs_addCol(o.model.inner, 0.0, Cdouble(set.lower), Cdouble(set.upper), Cint(0), Cint[], Cint[])
    return (MOI.VariableIndex(Int(idx)), MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}(idx))
end

function MOI.add_constraint(o::Optimizer, sg::MOI.SingleVariable, set::MOI.Interval)
    var_idx = sg.variable.value
    _ = CWrapper.Highs_changeColBounds(o.model.inner, var_idx, set.lower, set.upper)
    return nothing
end

function MOI.get(o::Optimizer, ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64}})
    nrows = CWrapper.Highs_getNumRows(o.model.inner)
    return Int(nrows)
end

function MOI.add_constraint(o::Optimizer, func::MOI.ScalarAffineFunction, set::MOI.Interval)
    number_nonzeros = length(func.terms)
    coefficients = Vector{Cdouble}(undef, number_nonzeros)
    col_indices = Vector{Cint}(undef, number_nonzeros)
    for j in Base.OneTo(number_nonzeros)
        term = func.terms[i]
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

const SUPPORTED_MODEL_ATTRIBUTES = Union{
    MOI.ObjectiveSense,
    MOI.NumberOfVariables,
    MOI.ListOfVariableIndices,
    MOI.ListOfConstraintIndices,
    MOI.NumberOfConstraints,
    MOI.ListOfConstraints,
    MOI.ObjectiveFunctionType,
    MOI.ObjectiveValue,
    MOI.DualObjectiveValue,
    # MOI.SolveTime,  # TODO
    MOI.SimplexIterations,
    MOI.BarrierIterations,
    MOI.RawSolver,
    # MOI.RawStatusString,  # TODO
    # MOI.ResultCount,
    # MOI.TerminationStatus,
    # MOI.PrimalStatus,
    # MOI.DualStatus
}

MOI.get(o::Optimizer, ::MOI.RawSolver) = o.model

# TODO get that?

function MOI.set(o::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    sense_code = sense == MOI.MAX_SENSE ? Cint(-1) : Cint(1)
    _ = CWrapper.Highs_changeObjectiveSense(o.model.inner, sense_code)
    return nothing
end

function MOI.get(o::Optimizer, ::MOI.NumberOfVariables)
    return Int(CWrapper.Highs_getNumCols(o.model.inner))
end

MOI.get(::Optimizer, ::MOI.ObjectiveFunctionType) = MOI.ScalarAffineFunction{Float64}

function MOI.get(o::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(o, attr)
    return CWrapper.Highs_getObjectiveValue(o.model.inner)
end

function MOI.get(o::Optimizer, ::MOI.SimplexIterations)
    res = Highs_getIterationCount(o.model.inner)
    return Int(res)
end

function MOI.get(o::Optimizer, ::MOI.BarrierIterations)
    return 0
end
