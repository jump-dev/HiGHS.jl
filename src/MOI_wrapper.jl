import MathOptInterface
const MOI = MathOptInterface

mutable struct Optimizer <: MOI.AbstractOptimizer
    model::ManagedHiGHS
    objective_sense::MOI.OptimizationSense
    Optimizer() = new(ManagedHiGHS(), MOI.FEASIBILITY_SENSE)
end

function MOI.empty!(o::Optimizer)
    reset_model!(o.model)
    o.objective_sense = MOI.FEASIBILITY_SENSE
    return nothing
end

function MOI.optimize!(o::Optimizer)
    CWrapper.Highs_run(o.model.inner)
    return nothing
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

MOI.supports(::Optimizer, param::MOI.RawParameter) = true

function MOI.set(o::Optimizer, param::MOI.RawParameter, value)
    # TODO check parameter
    CWrapper.Highs_setOptionValue(o, param.name, value)
    return nothing
end

const SUPPORTED_MODEL_ATTRIBUTES = Union{
    MOI.ObjectiveSense,
    MOI.NumberOfVariables,
    MOI.ListOfVariableIndices,
    MOI.ListOfConstraintIndices,
    MOI.NumberOfConstraints,
    MOI.ListOfConstraints,
    MOI.ObjectiveFunctionType,
    MOI.ObjectiveValue,
    MOI.DualObjectiveValue, # TODO
    MOI.SolveTime,  # TODO
    MOI.SimplexIterations,
    MOI.BarrierIterations,
    MOI.RawSolver,
    # MOI.RawStatusString,  # TODO
    MOI.ResultCount,
    # MOI.TerminationStatus,
    # MOI.PrimalStatus,
    # MOI.DualStatus
}

MOI.get(o::Optimizer, ::MOI.RawSolver) = o.model

function MOI.get(o::Optimizer, ::MOI.ResultCount)
    status = CWrapper.Highs_getModelStatus(o.model.inner)
    status == 11 ? 1 : 0
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

MOI.get(::Optimizer, ::MOI.ObjectiveFunctionType) = MOI.ScalarAffineFunction{Float64}

function MOI.set(o::Optimizer, ::MOI.ObjectiveFunction{F}, func::F) where {F <: MOI.ScalarAffineFunction{Float64}}
    # TODO treat constant in objective
    total_ncols = MOI.get(o, MOI.NumberOfVariables())
    coefficients  = zeros(Cdouble, total_ncols)

    for term in func.terms
        j = term.variable_index.value
        coefficients[j] = Cdouble(term.coefficient)
    end

    coefficients_ptr = pointer(coefficients)
    mask = pointer(ones(Cint, total_ncols))

    CWrapper.Highs_changeColsCostByMask(o.model.inner, mask, coefficients_ptr)
    return nothing
end

function MOI.get(o::Optimizer, ::MOI.ObjectiveFunction{F}) where {F <: MOI.ScalarAffineFunction{Float64}}
    total_ncols = MOI.get(o, MOI.NumberOfVariables())
    coefficients  = zeros(Cdouble, total_ncols)

    coefficients_ptr = pointer(coefficients)
    num_cols_obtained = Ref(Cint(0))
    num_nz = Ref(Cint(0))
    null_ptr = Ptr{Cint}()
    null_ptr_double = Ptr{Cdouble}()
    _ = CWrapper.Highs_getColsByRange(
        o.model.inner,
        Cint(0), Cint(total_ncols),
        pointer_from_objref(num_cols_obtained), coefficients_ptr,
        null_ptr_double, null_ptr_double, # lower, upper
        pointer_from_objref(num_nz), null_ptr, null_ptr, null_ptr_double
    )
    # TODO continue
    terms = Vector{MOI.ScalarAffineTerm{Float64}}()
    sizehint!(terms, total_ncols)
    for col_idx in Base.OneTo(total_ncols)
        if !≈(coefficients[col_idx], 0)
            push!(
                terms,
                MOI.ScalarAffineTerm{Float64}(coefficients[col_idx], MOI.VariableIndex(col_idx - 1))
            )
        end
    end
    return MOI.ScalarAffineFunction{Float64}(terms, 0.0)
end

function Highs_getColsByRange(highs, from_col::Cint, to_col::Cint, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getColsByRange, libhighs), Cint, (Ptr{Cvoid}, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), highs, from_col, to_col, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

function MOI.get(o::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(o, attr)
    res = CWrapper.Highs_getObjectiveValue(o.model.inner)
    # BUG linked to https://github.com/ERGO-Code/HiGHS/issues/209
    return o.objective_sense === MOI.MAX_SENSE ? -res : res
end

function MOI.get(o::Optimizer, ::MOI.SimplexIterations)
    res = Highs_getIterationCount(o.model.inner)
    return Int(res)
end

function MOI.get(o::Optimizer, ::MOI.BarrierIterations)
    return 0
end
