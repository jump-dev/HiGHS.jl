# Copyright (c) 2019 Mathieu Besançon, Oscar Dowson, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Disable JuliaFormatter for this file.
#!format:off

const HighsInt = Cint

struct HighsCallbackDataOut
    log_type::Cint
    running_time::Cdouble
    simplex_iteration_count::HighsInt
    ipm_iteration_count::HighsInt
    pdlp_iteration_count::HighsInt
    objective_function_value::Cdouble
    mip_node_count::Int64
    mip_primal_bound::Cdouble
    mip_dual_bound::Cdouble
    mip_gap::Cdouble
    mip_solution::Ptr{Cdouble}
    cutpool_num_col::HighsInt
    cutpool_num_cut::HighsInt
    cutpool_num_nz::HighsInt
    cutpool_start::Ptr{HighsInt}
    cutpool_index::Ptr{HighsInt}
    cutpool_value::Ptr{Cdouble}
    cutpool_lower::Ptr{Cdouble}
    cutpool_upper::Ptr{Cdouble}
end

struct HighsCallbackDataIn
    user_interrupt::Cint
end

# typedef void ( * HighsCCallbackType ) ( int , const char * , const HighsCallbackDataOut * , HighsCallbackDataIn * , void * )
const HighsCCallbackType = Ptr{Cvoid}

"""
    Highs_lpCall(num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, col_value, col_dual, row_value, row_dual, col_basis_status, row_basis_status, model_status)

Formulate and solve a linear program using HiGHS.

### Parameters
* `num_col`: The number of columns.
* `num_row`: The number of rows.
* `num_nz`: The number of nonzeros in the constraint matrix.
* `a_format`: The format of the constraint matrix as a `kHighsMatrixFormat` constant.
* `sense`: The optimization sense as a `kHighsObjSense` constant.
* `offset`: The objective constant.
* `col_cost`: An array of length [num\\_col] with the column costs.
* `col_lower`: An array of length [num\\_col] with the column lower bounds.
* `col_upper`: An array of length [num\\_col] with the column upper bounds.
* `row_lower`: An array of length [num\\_row] with the row lower bounds.
* `row_upper`: An array of length [num\\_row] with the row upper bounds.
* `a_start`: The constraint matrix is provided to HiGHS in compressed sparse column form (if `a_format` is `kHighsMatrixFormatColwise`, otherwise compressed sparse row form). The sparse matrix consists of three arrays, `a_start`, `a_index`, and `a_value`. `a_start` is an array of length [num\\_col] containing the starting index of each column in `a_index`. If `a_format` is `kHighsMatrixFormatRowwise` the array is of length [num\\_row] corresponding to each row.
* `a_index`: An array of length [num\\_nz] with indices of matrix entries.
* `a_value`: An array of length [num\\_nz] with values of matrix entries.
* `col_value`: An array of length [num\\_col], to be filled with the primal column solution.
* `col_dual`: An array of length [num\\_col], to be filled with the dual column solution.
* `row_value`: An array of length [num\\_row], to be filled with the primal row solution.
* `row_dual`: An array of length [num\\_row], to be filled with the dual row solution.
* `col_basis_status`: An array of length [num\\_col], to be filled with the basis status of the columns in the form of a `kHighsBasisStatus` constant.
* `row_basis_status`: An array of length [num\\_row], to be filled with the basis status of the rows in the form of a `kHighsBasisStatus` constant.
* `model_status`: The location in which to place the termination status of the model after the solve in the form of a `kHighsModelStatus` constant.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_lpCall(num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, col_value, col_dual, row_value, row_dual, col_basis_status, row_basis_status, model_status)
    ccall((:Highs_lpCall, libhighs), HighsInt, (HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}), num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, col_value, col_dual, row_value, row_dual, col_basis_status, row_basis_status, model_status)
end

"""
    Highs_mipCall(num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality, col_value, row_value, model_status)

Formulate and solve a mixed-integer linear program using HiGHS.

The signature of this method is identical to [`Highs_lpCall`](@ref), except that it has an additional `integrality` argument, and that it is missing the `col_dual`, `row_dual`, `col_basis_status` and `row_basis_status` arguments.

### Parameters
* `integrality`: An array of length [num\\_col], containing a `kHighsVarType` constant for each column.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_mipCall(num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality, col_value, row_value, model_status)
    ccall((:Highs_mipCall, libhighs), HighsInt, (HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}), num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality, col_value, row_value, model_status)
end

"""
    Highs_qpCall(num_col, num_row, num_nz, q_num_nz, a_format, q_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value, col_value, col_dual, row_value, row_dual, col_basis_status, row_basis_status, model_status)

Formulate and solve a quadratic program using HiGHS.

The signature of this method is identical to [`Highs_lpCall`](@ref), except that it has additional arguments for specifying the Hessian matrix.

### Parameters
* `q_num_nz`: The number of nonzeros in the Hessian matrix.
* `q_format`: The format of the Hessian matrix in the form of a `kHighsHessianStatus` constant. If q\\_num\\_nz > 0, this must be `kHighsHessianFormatTriangular`.
* `q_start`: The Hessian matrix is provided in the same format as the constraint matrix, using `q_start`, `q_index`, and `q_value` in the place of `a_start`, `a_index`, and `a_value`.
* `q_index`: An array of length [q\\_num\\_nz] with indices of matrix sentries.
* `q_value`: An array of length [q\\_num\\_nz] with values of matrix entries.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_qpCall(num_col, num_row, num_nz, q_num_nz, a_format, q_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value, col_value, col_dual, row_value, row_dual, col_basis_status, row_basis_status, model_status)
    ccall((:Highs_qpCall, libhighs), HighsInt, (HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}), num_col, num_row, num_nz, q_num_nz, a_format, q_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value, col_value, col_dual, row_value, row_dual, col_basis_status, row_basis_status, model_status)
end

"""
    Highs_create()

Create a Highs instance and return the reference.

Call [`Highs_destroy`](@ref) on the returned reference to clean up allocated memory.

### Returns
A pointer to the Highs instance.
"""
function Highs_create()
    ccall((:Highs_create, libhighs), Ptr{Cvoid}, ())
end

"""
    Highs_destroy(highs)

Destroy the model `highs` created by [`Highs_create`](@ref) and free all corresponding memory. Future calls using `highs` are not allowed.

To empty a model without invalidating `highs`, see [`Highs_clearModel`](@ref).

### Parameters
* `highs`: A pointer to the Highs instance.
"""
function Highs_destroy(highs)
    ccall((:Highs_destroy, libhighs), Cvoid, (Ptr{Cvoid},), highs)
end

"""
    Highs_version()

Return the HiGHS version number as a string of the form "vX.Y.Z".

### Returns
The HiGHS version as a `char*`.
"""
function Highs_version()
    ccall((:Highs_version, libhighs), Ptr{Cchar}, ())
end

"""
    Highs_versionMajor()

Return the HiGHS major version number.

### Returns
The HiGHS major version number.
"""
function Highs_versionMajor()
    ccall((:Highs_versionMajor, libhighs), HighsInt, ())
end

"""
    Highs_versionMinor()

Return the HiGHS minor version number.

### Returns
The HiGHS minor version number.
"""
function Highs_versionMinor()
    ccall((:Highs_versionMinor, libhighs), HighsInt, ())
end

"""
    Highs_versionPatch()

Return the HiGHS patch version number.

### Returns
The HiGHS patch version number.
"""
function Highs_versionPatch()
    ccall((:Highs_versionPatch, libhighs), HighsInt, ())
end

"""
    Highs_githash()

Return the HiGHS githash.

### Returns
The HiGHS githash.
"""
function Highs_githash()
    ccall((:Highs_githash, libhighs), Ptr{Cchar}, ())
end

"""
    Highs_compilationDate()

Return the HiGHS compilation date.

### Returns
Thse HiGHS compilation date.
"""
function Highs_compilationDate()
    ccall((:Highs_compilationDate, libhighs), Ptr{Cchar}, ())
end

"""
    Highs_readModel(highs, filename)

Read a model from `filename` into `highs`.

### Parameters
* `highs`: A pointer to the Highs instance.
* `filename`: The filename to read.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_readModel(highs, filename)
    ccall((:Highs_readModel, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}), highs, filename)
end

"""
    Highs_writeModel(highs, filename)

Write the model in `highs` to `filename`.

### Parameters
* `highs`: A pointer to the Highs instance.
* `filename`: The filename to write.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_writeModel(highs, filename)
    ccall((:Highs_writeModel, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}), highs, filename)
end

"""
    Highs_clear(highs)

Reset the options and then call `clearModel`.

See [`Highs_destroy`](@ref) to free all associated memory.

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_clear(highs)
    ccall((:Highs_clear, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_clearModel(highs)

Remove all variables and constraints from the model `highs`, but do not invalidate the pointer `highs`. Future calls (for example, adding new variables and constraints) are allowed.

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_clearModel(highs)
    ccall((:Highs_clearModel, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_clearSolver(highs)

Clear all solution data associated with the model.

See [`Highs_destroy`](@ref) to clear the model and free all associated memory.

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_clearSolver(highs)
    ccall((:Highs_clearSolver, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_presolve(highs)

Presolve a model.

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_presolve(highs)
    ccall((:Highs_presolve, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_run(highs)

Optimize a model. The algorithm used by HiGHS depends on the options that have been set.

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_run(highs)
    ccall((:Highs_run, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_postsolve(highs, col_value, col_dual, row_dual)

Postsolve a model using a primal (and possibly dual) solution.

### Parameters
* `highs`: A pointer to the Highs instance.
* `col_value`: An array of length [num\\_col] with the column solution values.
* `col_dual`: An array of length [num\\_col] with the column dual values, or a null pointer if not known.
* `row_dual`: An array of length [num\\_row] with the row dual values, or a null pointer if not known.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_postsolve(highs, col_value, col_dual, row_dual)
    ccall((:Highs_postsolve, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), highs, col_value, col_dual, row_dual)
end

"""
    Highs_writeSolution(highs, filename)

Write the solution information (including dual and basis status, if available) to a file.

See also: [`Highs_writeSolutionPretty`](@ref).

### Parameters
* `highs`: A pointer to the Highs instance.
* `filename`: The name of the file to write the results to.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_writeSolution(highs, filename)
    ccall((:Highs_writeSolution, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}), highs, filename)
end

"""
    Highs_writeSolutionPretty(highs, filename)

Write the solution information (including dual and basis status, if available) to a file in a human-readable format.

The method identical to [`Highs_writeSolution`](@ref), except that the printout is in a human-readable format.

### Parameters
* `highs`: A pointer to the Highs instance.
* `filename`: The name of the file to write the results to.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_writeSolutionPretty(highs, filename)
    ccall((:Highs_writeSolutionPretty, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}), highs, filename)
end

"""
    Highs_passLp(highs, num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value)

Pass a linear program (LP) to HiGHS in a single function call.

The signature of this function is identical to [`Highs_passModel`](@ref), without the arguments for passing the Hessian matrix of a quadratic program and the integrality vector.

### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_passLp(highs, num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value)
    ccall((:Highs_passLp, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value)
end

"""
    Highs_passMip(highs, num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality)

Pass a mixed-integer linear program (MILP) to HiGHS in a single function call.

The signature of function is identical to [`Highs_passModel`](@ref), without the arguments for passing the Hessian matrix of a quadratic program.

### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_passMip(highs, num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality)
    ccall((:Highs_passMip, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}), highs, num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality)
end

"""
    Highs_passModel(highs, num_col, num_row, num_nz, q_num_nz, a_format, q_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value, integrality)

Pass a model to HiGHS in a single function call. This is faster than constructing the model using [`Highs_addRow`](@ref) and [`Highs_addCol`](@ref).

### Parameters
* `highs`: A pointer to the Highs instance.
* `num_col`: The number of columns.
* `num_row`: The number of rows.
* `num_nz`: The number of elements in the constraint matrix.
* `q_num_nz`: The number of elements in the Hessian matrix.
* `a_format`: The format of the constraint matrix to use in the form of a `kHighsMatrixFormat` constant.
* `q_format`: The format of the Hessian matrix to use in the form of a `kHighsHessianFormat` constant.
* `sense`: The optimization sense in the form of a `kHighsObjSense` constant.
* `offset`: The constant term in the objective function.
* `col_cost`: An array of length [num\\_col] with the objective coefficients.
* `col_lower`: An array of length [num\\_col] with the lower column bounds.
* `col_upper`: An array of length [num\\_col] with the upper column bounds.
* `row_lower`: An array of length [num\\_row] with the upper row bounds.
* `row_upper`: An array of length [num\\_row] with the upper row bounds.
* `a_start`: The constraint matrix is provided to HiGHS in compressed sparse column form (if `a_format` is `kHighsMatrixFormatColwise`, otherwise compressed sparse row form). The sparse matrix consists of three arrays, `a_start`, `a_index`, and `a_value`. `a_start` is an array of length [num\\_col] containing the starting index of each column in `a_index`. If `a_format` is `kHighsMatrixFormatRowwise` the array is of length [num\\_row] corresponding to each row.
* `a_index`: An array of length [num\\_nz] with indices of matrix entries.
* `a_value`: An array of length [num\\_nz] with values of matrix entries.
* `q_start`: The Hessian matrix is provided in the same format as the constraint matrix, using `q_start`, `q_index`, and `q_value` in the place of `a_start`, `a_index`, and `a_value`. If the model is linear, pass NULL.
* `q_index`: An array of length [q\\_num\\_nz] with indices of matrix entries. If the model is linear, pass NULL.
* `q_value`: An array of length [q\\_num\\_nz] with values of matrix entries. If the model is linear, pass NULL.
* `integrality`: An array of length [num\\_col] containing a `kHighsVarType` constant for each column.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_passModel(highs, num_col, num_row, num_nz, q_num_nz, a_format, q_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value, integrality)
    ccall((:Highs_passModel, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}), highs, num_col, num_row, num_nz, q_num_nz, a_format, q_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value, integrality)
end

"""
    Highs_passHessian(highs, dim, num_nz, format, start, index, value)

Set the Hessian matrix for a quadratic objective.

### Parameters
* `highs`: A pointer to the Highs instance.
* `dim`: The dimension of the Hessian matrix. Should be [num\\_col].
* `num_nz`: The number of non-zero elements in the Hessian matrix.
* `format`: The format of the Hessian matrix as a `kHighsHessianFormat` constant. This must be `kHighsHessianFormatTriangular`.
* `start`: The Hessian matrix is provided to HiGHS as the upper triangular component in compressed sparse column form. The sparse matrix consists of three arrays, `start`, `index`, and `value`. `start` is an array of length [num\\_col] containing the starting index of each column in `index`.
* `index`: An array of length [num\\_nz] with indices of matrix entries.
* `value`: An array of length [num\\_nz] with values of matrix entries.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_passHessian(highs, dim, num_nz, format, start, index, value)
    ccall((:Highs_passHessian, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, dim, num_nz, format, start, index, value)
end

"""
    Highs_passRowName(highs, row, name)

Pass the name of a row.

### Parameters
* `highs`: A pointer to the Highs instance.
* `row`: The row for which the name is supplied.
* `name`: The name of the row.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_passRowName(highs, row, name)
    ccall((:Highs_passRowName, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cchar}), highs, row, name)
end

"""
    Highs_passColName(highs, col, name)

Pass the name of a column.

### Parameters
* `highs`: A pointer to the Highs instance.
* `col`: The column for which the name is supplied.
* `name`: The name of the column.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_passColName(highs, col, name)
    ccall((:Highs_passColName, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cchar}), highs, col, name)
end

"""
    Highs_readOptions(highs, filename)

Read the option values from file.

### Parameters
* `highs`: A pointer to the Highs instance.
* `filename`: The filename from which to read the option values.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_readOptions(highs, filename)
    ccall((:Highs_readOptions, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}), highs, filename)
end

"""
    Highs_setBoolOptionValue(highs, option, value)

Set a boolean-valued option.

### Parameters
* `highs`: A pointer to the Highs instance.
* `option`: The name of the option.
* `value`: The new value of the option.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_setBoolOptionValue(highs, option, value)
    ccall((:Highs_setBoolOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, HighsInt), highs, option, value)
end

"""
    Highs_setIntOptionValue(highs, option, value)

Set an int-valued option.

### Parameters
* `highs`: A pointer to the Highs instance.
* `option`: The name of the option.
* `value`: The new value of the option.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_setIntOptionValue(highs, option, value)
    ccall((:Highs_setIntOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, HighsInt), highs, option, value)
end

"""
    Highs_setDoubleOptionValue(highs, option, value)

Set a double-valued option.

### Parameters
* `highs`: A pointer to the Highs instance.
* `option`: The name of the option.
* `value`: The new value of the option.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_setDoubleOptionValue(highs, option, value)
    ccall((:Highs_setDoubleOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Cdouble), highs, option, value)
end

"""
    Highs_setStringOptionValue(highs, option, value)

Set a string-valued option.

### Parameters
* `highs`: A pointer to the Highs instance.
* `option`: The name of the option.
* `value`: The new value of the option.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_setStringOptionValue(highs, option, value)
    ccall((:Highs_setStringOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cchar}), highs, option, value)
end

"""
    Highs_getBoolOptionValue(highs, option, value)

Get a boolean-valued option.

### Parameters
* `highs`: A pointer to the Highs instance.
* `option`: The name of the option.
* `value`: The location in which the current value of the option should be placed.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getBoolOptionValue(highs, option, value)
    ccall((:Highs_getBoolOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{HighsInt}), highs, option, value)
end

"""
    Highs_getIntOptionValue(highs, option, value)

Get an int-valued option.

### Parameters
* `highs`: A pointer to the Highs instance.
* `option`: The name of the option.
* `value`: The location in which the current value of the option should be placed.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getIntOptionValue(highs, option, value)
    ccall((:Highs_getIntOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{HighsInt}), highs, option, value)
end

"""
    Highs_getDoubleOptionValue(highs, option, value)

Get a double-valued option.

### Parameters
* `highs`: A pointer to the Highs instance.
* `option`: The name of the option.
* `value`: The location in which the current value of the option should be placed.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getDoubleOptionValue(highs, option, value)
    ccall((:Highs_getDoubleOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cdouble}), highs, option, value)
end

"""
    Highs_getStringOptionValue(highs, option, value)

Get a string-valued option.

### Parameters
* `highs`: A pointer to the Highs instance.
* `option`: The name of the option.
* `value`: A pointer to allocated memory (of at least `kMaximumStringLength`) to store the current value of the option.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getStringOptionValue(highs, option, value)
    ccall((:Highs_getStringOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cchar}), highs, option, value)
end

"""
    Highs_getOptionType(highs, option, type)

Get the type expected by an option.

### Parameters
* `highs`: A pointer to the Highs instance.
* `option`: The name of the option.
* `type`: An int in which the corresponding `kHighsOptionType` constant should be placed.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getOptionType(highs, option, type)
    ccall((:Highs_getOptionType, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{HighsInt}), highs, option, type)
end

"""
    Highs_resetOptions(highs)

Reset all options to their default value.

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_resetOptions(highs)
    ccall((:Highs_resetOptions, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_writeOptions(highs, filename)

Write the current options to file.

### Parameters
* `highs`: A pointer to the Highs instance.
* `filename`: The filename to write the options to.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_writeOptions(highs, filename)
    ccall((:Highs_writeOptions, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}), highs, filename)
end

"""
    Highs_writeOptionsDeviations(highs, filename)

Write the value of non-default options to file.

This is similar to [`Highs_writeOptions`](@ref), except only options with non-default value are written to `filename`.

### Parameters
* `highs`: A pointer to the Highs instance.
* `filename`: The filename to write the options to.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_writeOptionsDeviations(highs, filename)
    ccall((:Highs_writeOptionsDeviations, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}), highs, filename)
end

"""
    Highs_getNumOptions(highs)

Return the number of options

### Parameters
* `highs`: A pointer to the Highs instance.
"""
function Highs_getNumOptions(highs)
    ccall((:Highs_getNumOptions, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_getOptionName(highs, index, name)

Get the name of an option identified by index

### Parameters
* `highs`: A pointer to the Highs instance.
* `index`: The index of the option.
* `name`: The name of the option.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getOptionName(highs, index, name)
    ccall((:Highs_getOptionName, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Ptr{Cchar}}), highs, index, name)
end

"""
    Highs_getBoolOptionValues(highs, option, current_value, default_value)

Get the current and default values of a bool option

### Parameters
* `highs`: A pointer to the Highs instance.
* `current_value`: A pointer to the current value of the option.
* `default_value`: A pointer to the default value of the option.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getBoolOptionValues(highs, option, current_value, default_value)
    ccall((:Highs_getBoolOptionValues, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{HighsInt}, Ptr{HighsInt}), highs, option, current_value, default_value)
end

"""
    Highs_getIntOptionValues(highs, option, current_value, min_value, max_value, default_value)

Get the current and default values of an int option

### Parameters
* `highs`: A pointer to the Highs instance.
* `current_value`: A pointer to the current value of the option.
* `min_value`: A pointer to the minimum value of the option.
* `max_value`: A pointer to the maximum value of the option.
* `default_value`: A pointer to the default value of the option.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getIntOptionValues(highs, option, current_value, min_value, max_value, default_value)
    ccall((:Highs_getIntOptionValues, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}), highs, option, current_value, min_value, max_value, default_value)
end

"""
    Highs_getDoubleOptionValues(highs, option, current_value, min_value, max_value, default_value)

Get the current and default values of a double option

### Parameters
* `highs`: A pointer to the Highs instance.
* `current_value`: A pointer to the current value of the option.
* `min_value`: A pointer to the minimum value of the option.
* `max_value`: A pointer to the maximum value of the option.
* `default_value`: A pointer to the default value of the option.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getDoubleOptionValues(highs, option, current_value, min_value, max_value, default_value)
    ccall((:Highs_getDoubleOptionValues, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), highs, option, current_value, min_value, max_value, default_value)
end

"""
    Highs_getStringOptionValues(highs, option, current_value, default_value)

Get the current and default values of a string option

### Parameters
* `highs`: A pointer to the Highs instance.
* `current_value`: A pointer to the current value of the option.
* `default_value`: A pointer to the default value of the option.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getStringOptionValues(highs, option, current_value, default_value)
    ccall((:Highs_getStringOptionValues, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}), highs, option, current_value, default_value)
end

"""
    Highs_getIntInfoValue(highs, info, value)

Get an int-valued info value.

### Parameters
* `highs`: A pointer to the Highs instance.
* `info`: The name of the info item.
* `value`: A reference to an integer that the result will be stored in.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getIntInfoValue(highs, info, value)
    ccall((:Highs_getIntInfoValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{HighsInt}), highs, info, value)
end

"""
    Highs_getDoubleInfoValue(highs, info, value)

Get a double-valued info value.

### Parameters
* `highs`: A pointer to the Highs instance.
* `info`: The name of the info item.
* `value`: A reference to a double that the result will be stored in.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getDoubleInfoValue(highs, info, value)
    ccall((:Highs_getDoubleInfoValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cdouble}), highs, info, value)
end

"""
    Highs_getInt64InfoValue(highs, info, value)

Get an int64-valued info value.

### Parameters
* `highs`: A pointer to the Highs instance.
* `info`: The name of the info item.
* `value`: A reference to an int64 that the result will be stored in.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getInt64InfoValue(highs, info, value)
    ccall((:Highs_getInt64InfoValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Int64}), highs, info, value)
end

"""
    Highs_getInfoType(highs, info, type)

Get the type expected by an info item.

### Parameters
* `highs`: A pointer to the Highs instance.
* `info`: The name of the info item.
* `type`: An int in which the corresponding `kHighsOptionType` constant is stored.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getInfoType(highs, info, type)
    ccall((:Highs_getInfoType, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{HighsInt}), highs, info, type)
end

"""
    Highs_getSolution(highs, col_value, col_dual, row_value, row_dual)

Get the primal and dual solution from an optimized model.

### Parameters
* `highs`: A pointer to the Highs instance.
* `col_value`: An array of length [num\\_col], to be filled with primal column values.
* `col_dual`: An array of length [num\\_col], to be filled with dual column values.
* `row_value`: An array of length [num\\_row], to be filled with primal row values.
* `row_dual`: An array of length [num\\_row], to be filled with dual row values.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getSolution(highs, col_value, col_dual, row_value, row_dual)
    ccall((:Highs_getSolution, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), highs, col_value, col_dual, row_value, row_dual)
end

"""
    Highs_getBasis(highs, col_status, row_status)

Given a linear program with a basic feasible solution, get the column and row basis statuses.

### Parameters
* `highs`: A pointer to the Highs instance.
* `col_status`: An array of length [num\\_col], to be filled with the column basis statuses in the form of a `kHighsBasisStatus` constant.
* `row_status`: An array of length [num\\_row], to be filled with the row basis statuses in the form of a `kHighsBasisStatus` constant.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getBasis(highs, col_status, row_status)
    ccall((:Highs_getBasis, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{HighsInt}), highs, col_status, row_status)
end

"""
    Highs_getModelStatus(highs)

Return the optimization status of the model in the form of a `kHighsModelStatus` constant.

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
An integer corresponding to the `kHighsModelStatus` constant
"""
function Highs_getModelStatus(highs)
    ccall((:Highs_getModelStatus, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_getDualRay(highs, has_dual_ray, dual_ray_value)

Get an unbounded dual ray that is a certificate of primal infeasibility.

### Parameters
* `highs`: A pointer to the Highs instance.
* `has_dual_ray`: A pointer to an int to store 1 if the dual ray exists.
* `dual_ray_value`: An array of length [num\\_row] filled with the unbounded ray.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getDualRay(highs, has_dual_ray, dual_ray_value)
    ccall((:Highs_getDualRay, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{Cdouble}), highs, has_dual_ray, dual_ray_value)
end

"""
    Highs_getPrimalRay(highs, has_primal_ray, primal_ray_value)

Get an unbounded primal ray that is a certificate of dual infeasibility.

### Parameters
* `highs`: A pointer to the Highs instance.
* `has_primal_ray`: A pointer to an int to store 1 if the primal ray exists.
* `primal_ray_value`: An array of length [num\\_col] filled with the unbounded ray.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getPrimalRay(highs, has_primal_ray, primal_ray_value)
    ccall((:Highs_getPrimalRay, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{Cdouble}), highs, has_primal_ray, primal_ray_value)
end

"""
    Highs_getObjectiveValue(highs)

Get the primal objective function value.

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
The primal objective function value
"""
function Highs_getObjectiveValue(highs)
    ccall((:Highs_getObjectiveValue, libhighs), Cdouble, (Ptr{Cvoid},), highs)
end

"""
    Highs_getBasicVariables(highs, basic_variables)

Get the indices of the rows and columns that make up the basis matrix ``B`` of a basic feasible solution.

Non-negative entries are indices of columns, and negative entries are `-row\\_index - 1`. For example, `{1, -1}` would be the second column and first row.

The order of these rows and columns is important for calls to the functions:

- [`Highs_getBasisInverseRow`](@ref) - [`Highs_getBasisInverseCol`](@ref) - [`Highs_getBasisSolve`](@ref) - [`Highs_getBasisTransposeSolve`](@ref) - [`Highs_getReducedRow`](@ref) - [`Highs_getReducedColumn`](@ref)

### Parameters
* `highs`: A pointer to the Highs instance.
* `basic_variables`: An array of size [num\\_rows], filled with the indices of the basic variables.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getBasicVariables(highs, basic_variables)
    ccall((:Highs_getBasicVariables, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}), highs, basic_variables)
end

"""
    Highs_getBasisInverseRow(highs, row, row_vector, row_num_nz, row_index)

Get a row of the inverse basis matrix ``B^{-1}``.

See [`Highs_getBasicVariables`](@ref) for a description of the ``B`` matrix.

The arrays `row_vector` and `row_index` must have an allocated length of [num\\_row]. However, check `row_num_nz` to see how many non-zero elements are actually stored.

### Parameters
* `highs`: A pointer to the Highs instance.
* `row`: The index of the row to compute.
* `row_vector`: An array of length [num\\_row] in which to store the values of the non-zero elements.
* `row_num_nz`: The number of non-zeros in the row.
* `row_index`: An array of length [num\\_row] in which to store the indices of the non-zero elements.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getBasisInverseRow(highs, row, row_vector, row_num_nz, row_index)
    ccall((:Highs_getBasisInverseRow, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, row, row_vector, row_num_nz, row_index)
end

"""
    Highs_getBasisInverseCol(highs, col, col_vector, col_num_nz, col_index)

Get a column of the inverse basis matrix ``B^{-1}``.

See [`Highs_getBasicVariables`](@ref) for a description of the ``B`` matrix.

The arrays `col_vector` and `col_index` must have an allocated length of [num\\_row]. However, check `col_num_nz` to see how many non-zero elements are actually stored.

### Parameters
* `highs`: A pointer to the Highs instance.
* `col`: The index of the column to compute.
* `col_vector`: An array of length [num\\_row] in which to store the values of the non-zero elements.
* `col_num_nz`: The number of non-zeros in the column.
* `col_index`: An array of length [num\\_row] in which to store the indices of the non-zero elements.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getBasisInverseCol(highs, col, col_vector, col_num_nz, col_index)
    ccall((:Highs_getBasisInverseCol, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, col, col_vector, col_num_nz, col_index)
end

"""
    Highs_getBasisSolve(highs, rhs, solution_vector, solution_num_nz, solution_index)

Compute ``{x}=B^{-1}{b}`` for a given vector ``{b}``.

See [`Highs_getBasicVariables`](@ref) for a description of the ``B`` matrix.

The arrays `solution_vector` and `solution_index` must have an allocated length of [num\\_row]. However, check `solution_num_nz` to see how many non-zero elements are actually stored.

### Parameters
* `highs`: A pointer to the Highs instance.
* `rhs`: The right-hand side vector ``b``.
* `solution_vector`: An array of length [num\\_row] in which to store the values of the non-zero elements.
* `solution_num_nz`: The number of non-zeros in the solution.
* `solution_index`: An array of length [num\\_row] in which to store the indices of the non-zero elements.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getBasisSolve(highs, rhs, solution_vector, solution_num_nz, solution_index)
    ccall((:Highs_getBasisSolve, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, rhs, solution_vector, solution_num_nz, solution_index)
end

"""
    Highs_getBasisTransposeSolve(highs, rhs, solution_vector, solution_nz, solution_index)

Compute ``{x}=B^{-T}{b}`` for a given vector ``{b}``.

See [`Highs_getBasicVariables`](@ref) for a description of the ``B`` matrix.

The arrays `solution_vector` and `solution_index` must have an allocated length of [num\\_row]. However, check `solution_num_nz` to see how many non-zero elements are actually stored.

### Parameters
* `highs`: A pointer to the Highs instance.
* `rhs`: The right-hand side vector ``b``
* `solution_vector`: An array of length [num\\_row] in which to store the values of the non-zero elements.
* `solution_num_nz`: The number of non-zeros in the solution.
* `solution_index`: An array of length [num\\_row] in which to store the indices of the non-zero elements.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getBasisTransposeSolve(highs, rhs, solution_vector, solution_nz, solution_index)
    ccall((:Highs_getBasisTransposeSolve, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, rhs, solution_vector, solution_nz, solution_index)
end

"""
    Highs_getReducedRow(highs, row, row_vector, row_num_nz, row_index)

Compute a row of ``B^{-1}A``.

See [`Highs_getBasicVariables`](@ref) for a description of the ``B`` matrix.

The arrays `row_vector` and `row_index` must have an allocated length of [num\\_row]. However, check `row_num_nz` to see how many non-zero elements are actually stored.

### Parameters
* `highs`: A pointer to the Highs instance.
* `row`: The index of the row to compute.
* `row_vector`: An array of length [num\\_row] in which to store the values of the non-zero elements.
* `row_num_nz`: The number of non-zeros in the row.
* `row_index`: An array of length [num\\_row] in which to store the indices of the non-zero elements.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getReducedRow(highs, row, row_vector, row_num_nz, row_index)
    ccall((:Highs_getReducedRow, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, row, row_vector, row_num_nz, row_index)
end

"""
    Highs_getReducedColumn(highs, col, col_vector, col_num_nz, col_index)

Compute a column of ``B^{-1}A``.

See [`Highs_getBasicVariables`](@ref) for a description of the ``B`` matrix.

The arrays `col_vector` and `col_index` must have an allocated length of [num\\_row]. However, check `col_num_nz` to see how many non-zero elements are actually stored.

### Parameters
* `highs`: A pointer to the Highs instance.
* `col`: The index of the column to compute.
* `col_vector`: An array of length [num\\_row] in which to store the values of the non-zero elements.
* `col_num_nz`: The number of non-zeros in the column.
* `col_index`: An array of length [num\\_row] in which to store the indices of the non-zero elements.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getReducedColumn(highs, col, col_vector, col_num_nz, col_index)
    ccall((:Highs_getReducedColumn, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, col, col_vector, col_num_nz, col_index)
end

"""
    Highs_setBasis(highs, col_status, row_status)

Set a basic feasible solution by passing the column and row basis statuses to the model.

### Parameters
* `highs`: A pointer to the Highs instance.
* `col_status`: an array of length [num\\_col] with the column basis status in the form of `kHighsBasisStatus` constants
* `row_status`: an array of length [num\\_row] with the row basis status in the form of `kHighsBasisStatus` constants
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_setBasis(highs, col_status, row_status)
    ccall((:Highs_setBasis, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{HighsInt}), highs, col_status, row_status)
end

"""
    Highs_setLogicalBasis(highs)

Set a logical basis in the model.

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_setLogicalBasis(highs)
    ccall((:Highs_setLogicalBasis, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_setSolution(highs, col_value, row_value, col_dual, row_dual)

Set a solution by passing the column and row primal and dual solution values.

For any values that are unavailable, pass NULL.

### Parameters
* `highs`: A pointer to the Highs instance.
* `col_value`: An array of length [num\\_col] with the column solution values.
* `row_value`: An array of length [num\\_row] with the row solution values.
* `col_dual`: An array of length [num\\_col] with the column dual values.
* `row_dual`: An array of length [num\\_row] with the row dual values.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_setSolution(highs, col_value, row_value, col_dual, row_dual)
    ccall((:Highs_setSolution, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), highs, col_value, row_value, col_dual, row_dual)
end

"""
    Highs_setCallback(highs, user_callback, user_callback_data)

Set the callback method to use for HiGHS

### Parameters
* `highs`: A pointer to the Highs instance.
* `user_callback`: A pointer to the user callback
* `user_callback_data`: A pointer to the user callback data
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_setCallback(highs, user_callback, user_callback_data)
    ccall((:Highs_setCallback, libhighs), HighsInt, (Ptr{Cvoid}, HighsCCallbackType, Ptr{Cvoid}), highs, user_callback, user_callback_data)
end

"""
    Highs_startCallback(highs, callback_type)

Start callback of given type

### Parameters
* `highs`: A pointer to the Highs instance.
* `callback_type`: The type of callback to be started
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_startCallback(highs, callback_type)
    ccall((:Highs_startCallback, libhighs), HighsInt, (Ptr{Cvoid}, Cint), highs, callback_type)
end

"""
    Highs_stopCallback(highs, callback_type)

Stop callback of given type

### Parameters
* `highs`: A pointer to the Highs instance.
* `callback_type`: The type of callback to be stopped
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_stopCallback(highs, callback_type)
    ccall((:Highs_stopCallback, libhighs), HighsInt, (Ptr{Cvoid}, Cint), highs, callback_type)
end

"""
    Highs_getRunTime(highs)

Return the cumulative wall-clock time spent in [`Highs_run`](@ref).

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
The cumulative wall-clock time spent in [`Highs_run`](@ref)
"""
function Highs_getRunTime(highs)
    ccall((:Highs_getRunTime, libhighs), Cdouble, (Ptr{Cvoid},), highs)
end

"""
    Highs_zeroAllClocks(highs)

Reset the clocks in a `highs` model.

Each `highs` model contains a single instance of clock that records how much time is spent in various parts of the algorithm. This clock is not reset on entry to [`Highs_run`](@ref), so repeated calls to [`Highs_run`](@ref) report the cumulative time spent in the algorithm. A side-effect is that this will trigger a time limit termination once the cumulative run time exceeds the time limit, rather than the run time of each individual call to [`Highs_run`](@ref).

As a work-around, call [`Highs_zeroAllClocks`](@ref) before each call to [`Highs_run`](@ref).

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_zeroAllClocks(highs)
    ccall((:Highs_zeroAllClocks, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_addCol(highs, cost, lower, upper, num_new_nz, index, value)

Add a new column (variable) to the model.

### Parameters
* `highs`: A pointer to the Highs instance.
* `cost`: The objective coefficient of the column.
* `lower`: The lower bound of the column.
* `upper`: The upper bound of the column.
* `num_new_nz`: The number of non-zeros in the column.
* `index`: An array of size [num\\_new\\_nz] with the row indices.
* `value`: An array of size [num\\_new\\_nz] with row values.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_addCol(highs, cost, lower, upper, num_new_nz, index, value)
    ccall((:Highs_addCol, libhighs), HighsInt, (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}), highs, cost, lower, upper, num_new_nz, index, value)
end

"""
    Highs_addCols(highs, num_new_col, costs, lower, upper, num_new_nz, starts, index, value)

Add multiple columns (variables) to the model.

### Parameters
* `highs`: A pointer to the Highs instance.
* `num_new_col`: The number of new columns to add.
* `costs`: An array of size [num\\_new\\_col] with objective coefficients.
* `lower`: An array of size [num\\_new\\_col] with lower bounds.
* `upper`: An array of size [num\\_new\\_col] with upper bounds.
* `num_new_nz`: The number of new nonzeros in the constraint matrix.
* `starts`: The constraint coefficients are given as a matrix in compressed sparse column form by the arrays `starts`, `index`, and `value`. `starts` is an array of size [num\\_new\\_cols] with the start index of each row in indices and values.
* `index`: An array of size [num\\_new\\_nz] with row indices.
* `value`: An array of size [num\\_new\\_nz] with row values.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_addCols(highs, num_new_col, costs, lower, upper, num_new_nz, starts, index, value)
    ccall((:Highs_addCols, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, num_new_col, costs, lower, upper, num_new_nz, starts, index, value)
end

"""
    Highs_addVar(highs, lower, upper)

Add a new variable to the model.

### Parameters
* `highs`: A pointer to the Highs instance.
* `lower`: The lower bound of the column.
* `upper`: The upper bound of the column.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_addVar(highs, lower, upper)
    ccall((:Highs_addVar, libhighs), HighsInt, (Ptr{Cvoid}, Cdouble, Cdouble), highs, lower, upper)
end

"""
    Highs_addVars(highs, num_new_var, lower, upper)

Add multiple variables to the model.

### Parameters
* `highs`: A pointer to the Highs instance.
* `num_new_var`: The number of new variables to add.
* `lower`: An array of size [num\\_new\\_var] with lower bounds.
* `upper`: An array of size [num\\_new\\_var] with upper bounds.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_addVars(highs, num_new_var, lower, upper)
    ccall((:Highs_addVars, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{Cdouble}), highs, num_new_var, lower, upper)
end

"""
    Highs_addRow(highs, lower, upper, num_new_nz, index, value)

Add a new row (a linear constraint) to the model.

### Parameters
* `highs`: A pointer to the Highs instance.
* `lower`: The lower bound of the row.
* `upper`: The upper bound of the row.
* `num_new_nz`: The number of non-zeros in the row
* `index`: An array of size [num\\_new\\_nz] with column indices.
* `value`: An array of size [num\\_new\\_nz] with column values.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_addRow(highs, lower, upper, num_new_nz, index, value)
    ccall((:Highs_addRow, libhighs), HighsInt, (Ptr{Cvoid}, Cdouble, Cdouble, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}), highs, lower, upper, num_new_nz, index, value)
end

"""
    Highs_addRows(highs, num_new_row, lower, upper, num_new_nz, starts, index, value)

Add multiple rows (linear constraints) to the model.

### Parameters
* `highs`: A pointer to the Highs instance.
* `num_new_row`: The number of new rows to add
* `lower`: An array of size [num\\_new\\_row] with the lower bounds of the rows.
* `upper`: An array of size [num\\_new\\_row] with the upper bounds of the rows.
* `num_new_nz`: The number of non-zeros in the rows.
* `starts`: The constraint coefficients are given as a matrix in compressed sparse row form by the arrays `starts`, `index`, and `value`. `starts` is an array of size [num\\_new\\_rows] with the start index of each row in indices and values.
* `index`: An array of size [num\\_new\\_nz] with column indices.
* `value`: An array of size [num\\_new\\_nz] with column values.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_addRows(highs, num_new_row, lower, upper, num_new_nz, starts, index, value)
    ccall((:Highs_addRows, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{Cdouble}, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, num_new_row, lower, upper, num_new_nz, starts, index, value)
end

"""
    Highs_changeObjectiveSense(highs, sense)

Change the objective sense of the model.

### Parameters
* `highs`: A pointer to the Highs instance.
* `sense`: The new optimization sense in the form of a `kHighsObjSense` constant.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_changeObjectiveSense(highs, sense)
    ccall((:Highs_changeObjectiveSense, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt), highs, sense)
end

"""
    Highs_changeObjectiveOffset(highs, offset)

Change the objective offset of the model.

### Parameters
* `highs`: A pointer to the Highs instance.
* `offset`: The new objective offset.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_changeObjectiveOffset(highs, offset)
    ccall((:Highs_changeObjectiveOffset, libhighs), HighsInt, (Ptr{Cvoid}, Cdouble), highs, offset)
end

"""
    Highs_changeColIntegrality(highs, col, integrality)

Change the integrality of a column.

### Parameters
* `highs`: A pointer to the Highs instance.
* `col`: The column index to change.
* `integrality`: The new integrality of the column in the form of a `kHighsVarType` constant.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_changeColIntegrality(highs, col, integrality)
    ccall((:Highs_changeColIntegrality, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt), highs, col, integrality)
end

"""
    Highs_changeColsIntegralityByRange(highs, from_col, to_col, integrality)

Change the integrality of multiple adjacent columns.

### Parameters
* `highs`: A pointer to the Highs instance.
* `from_col`: The index of the first column whose integrality changes.
* `to_col`: The index of the last column whose integrality changes.
* `integrality`: An array of length [to\\_col - from\\_col + 1] with the new integralities of the columns in the form of `kHighsVarType` constants.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_changeColsIntegralityByRange(highs, from_col, to_col, integrality)
    ccall((:Highs_changeColsIntegralityByRange, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Ptr{HighsInt}), highs, from_col, to_col, integrality)
end

"""
    Highs_changeColsIntegralityBySet(highs, num_set_entries, set, integrality)

Change the integrality of multiple columns given by an array of indices.

### Parameters
* `highs`: A pointer to the Highs instance.
* `num_set_entries`: The number of columns to change.
* `set`: An array of size [num\\_set\\_entries] with the indices of the columns to change.
* `integrality`: An array of length [num\\_set\\_entries] with the new integralities of the columns in the form of `kHighsVarType` constants.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_changeColsIntegralityBySet(highs, num_set_entries, set, integrality)
    ccall((:Highs_changeColsIntegralityBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}), highs, num_set_entries, set, integrality)
end

"""
    Highs_changeColsIntegralityByMask(highs, mask, integrality)

Change the integrality of multiple columns given by a mask.

### Parameters
* `highs`: A pointer to the Highs instance.
* `mask`: An array of length [num\\_col] with 1 if the column integrality should be changed and 0 otherwise.
* `integrality`: An array of length [num\\_col] with the new integralities of the columns in the form of `kHighsVarType` constants.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_changeColsIntegralityByMask(highs, mask, integrality)
    ccall((:Highs_changeColsIntegralityByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{HighsInt}), highs, mask, integrality)
end

"""
    Highs_clearIntegrality(highs)

Clear the integrality of all columns

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_clearIntegrality(highs)
    ccall((:Highs_clearIntegrality, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_changeColCost(highs, col, cost)

Change the objective coefficient of a column.

### Parameters
* `highs`: A pointer to the Highs instance.
* `col`: The index of the column fo change.
* `cost`: The new objective coefficient.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_changeColCost(highs, col, cost)
    ccall((:Highs_changeColCost, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Cdouble), highs, col, cost)
end

"""
    Highs_changeColsCostByRange(highs, from_col, to_col, cost)

Change the cost coefficients of multiple adjacent columns.

### Parameters
* `highs`: A pointer to the Highs instance.
* `from_col`: The index of the first column whose cost changes.
* `to_col`: The index of the last column whose cost changes.
* `cost`: An array of length [to\\_col - from\\_col + 1] with the new objective coefficients.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_changeColsCostByRange(highs, from_col, to_col, cost)
    ccall((:Highs_changeColsCostByRange, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Ptr{Cdouble}), highs, from_col, to_col, cost)
end

"""
    Highs_changeColsCostBySet(highs, num_set_entries, set, cost)

Change the cost of multiple columns given by an array of indices.

### Parameters
* `highs`: A pointer to the Highs instance.
* `num_set_entries`: The number of columns to change.
* `set`: An array of size [num\\_set\\_entries] with the indices of the columns to change.
* `cost`: An array of length [num\\_set\\_entries] with the new costs of the columns.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_changeColsCostBySet(highs, num_set_entries, set, cost)
    ccall((:Highs_changeColsCostBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}), highs, num_set_entries, set, cost)
end

"""
    Highs_changeColsCostByMask(highs, mask, cost)

Change the cost of multiple columns given by a mask.

### Parameters
* `highs`: A pointer to the Highs instance.
* `mask`: An array of length [num\\_col] with 1 if the column cost should be changed and 0 otherwise.
* `cost`: An array of length [num\\_col] with the new costs.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_changeColsCostByMask(highs, mask, cost)
    ccall((:Highs_changeColsCostByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{Cdouble}), highs, mask, cost)
end

"""
    Highs_changeColBounds(highs, col, lower, upper)

Change the variable bounds of a column.

### Parameters
* `highs`: A pointer to the Highs instance.
* `col`: The index of the column whose bounds are to change.
* `lower`: The new lower bound.
* `upper`: The new upper bound.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_changeColBounds(highs, col, lower, upper)
    ccall((:Highs_changeColBounds, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Cdouble, Cdouble), highs, col, lower, upper)
end

"""
    Highs_changeColsBoundsByRange(highs, from_col, to_col, lower, upper)

Change the variable bounds of multiple adjacent columns.

### Parameters
* `highs`: A pointer to the Highs instance.
* `from_col`: The index of the first column whose bound changes.
* `to_col`: The index of the last column whose bound changes.
* `lower`: An array of length [to\\_col - from\\_col + 1] with the new lower bounds.
* `upper`: An array of length [to\\_col - from\\_col + 1] with the new upper bounds.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_changeColsBoundsByRange(highs, from_col, to_col, lower, upper)
    ccall((:Highs_changeColsBoundsByRange, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Ptr{Cdouble}, Ptr{Cdouble}), highs, from_col, to_col, lower, upper)
end

"""
    Highs_changeColsBoundsBySet(highs, num_set_entries, set, lower, upper)

Change the bounds of multiple columns given by an array of indices.

### Parameters
* `highs`: A pointer to the Highs instance.
* `num_set_entries`: The number of columns to change.
* `set`: An array of size [num\\_set\\_entries] with the indices of the columns to change.
* `lower`: An array of length [num\\_set\\_entries] with the new lower bounds.
* `upper`: An array of length [num\\_set\\_entries] with the new upper bounds.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_changeColsBoundsBySet(highs, num_set_entries, set, lower, upper)
    ccall((:Highs_changeColsBoundsBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}), highs, num_set_entries, set, lower, upper)
end

"""
    Highs_changeColsBoundsByMask(highs, mask, lower, upper)

Change the variable bounds of multiple columns given by a mask.

### Parameters
* `highs`: A pointer to the Highs instance.
* `mask`: An array of length [num\\_col] with 1 if the column bounds should be changed and 0 otherwise.
* `lower`: An array of length [num\\_col] with the new lower bounds.
* `upper`: An array of length [num\\_col] with the new upper bounds.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_changeColsBoundsByMask(highs, mask, lower, upper)
    ccall((:Highs_changeColsBoundsByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}), highs, mask, lower, upper)
end

"""
    Highs_changeRowBounds(highs, row, lower, upper)

Change the bounds of a row.

### Parameters
* `highs`: A pointer to the Highs instance.
* `row`: The index of the row whose bounds are to change.
* `lower`: The new lower bound.
* `upper`: The new upper bound.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_changeRowBounds(highs, row, lower, upper)
    ccall((:Highs_changeRowBounds, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Cdouble, Cdouble), highs, row, lower, upper)
end

"""
    Highs_changeRowsBoundsBySet(highs, num_set_entries, set, lower, upper)

Change the bounds of multiple rows given by an array of indices.

### Parameters
* `highs`: A pointer to the Highs instance.
* `num_set_entries`: The number of rows to change.
* `set`: An array of size [num\\_set\\_entries] with the indices of the rows to change.
* `lower`: An array of length [num\\_set\\_entries] with the new lower bounds.
* `upper`: An array of length [num\\_set\\_entries] with the new upper bounds.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_changeRowsBoundsBySet(highs, num_set_entries, set, lower, upper)
    ccall((:Highs_changeRowsBoundsBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}), highs, num_set_entries, set, lower, upper)
end

"""
    Highs_changeRowsBoundsByMask(highs, mask, lower, upper)

Change the bounds of multiple rows given by a mask.

### Parameters
* `highs`: A pointer to the Highs instance.
* `mask`: An array of length [num\\_row] with 1 if the row bounds should be changed and 0 otherwise.
* `lower`: An array of length [num\\_row] with the new lower bounds.
* `upper`: An array of length [num\\_row] with the new upper bounds.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_changeRowsBoundsByMask(highs, mask, lower, upper)
    ccall((:Highs_changeRowsBoundsByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}), highs, mask, lower, upper)
end

"""
    Highs_changeCoeff(highs, row, col, value)

Change a coefficient in the constraint matrix.

### Parameters
* `highs`: A pointer to the Highs instance.
* `row`: The index of the row to change.
* `col`: The index of the column to change.
* `value`: The new constraint coefficient.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_changeCoeff(highs, row, col, value)
    ccall((:Highs_changeCoeff, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Cdouble), highs, row, col, value)
end

"""
    Highs_getObjectiveSense(highs, sense)

Get the objective sense.

### Parameters
* `highs`: A pointer to the Highs instance.
* `sense`: The location in which the current objective sense should be placed. The sense is a `kHighsObjSense` constant.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getObjectiveSense(highs, sense)
    ccall((:Highs_getObjectiveSense, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}), highs, sense)
end

"""
    Highs_getObjectiveOffset(highs, offset)

Get the objective offset.

### Parameters
* `highs`: A pointer to the Highs instance.
* `offset`: The location in which the current objective offset should be placed.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getObjectiveOffset(highs, offset)
    ccall((:Highs_getObjectiveOffset, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cdouble}), highs, offset)
end

"""
    Highs_getColsByRange(highs, from_col, to_col, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)

Get data associated with multiple adjacent columns from the model.

To query the constraint coefficients, this function should be called twice.

First, call this function with `matrix_start`, `matrix_index`, and `matrix_value` as `NULL`. This call will populate `num_nz` with the number of nonzero elements in the corresponding section of the constraint matrix.

Second, allocate new `matrix_index` and `matrix_value` arrays of length `num_nz` and call this function again to populate the new arrays with their contents.

### Parameters
* `highs`: A pointer to the Highs instance.
* `from_col`: The first column for which to query data for.
* `to_col`: The last column (inclusive) for which to query data for.
* `num_col`: An integer populated with the number of columns got from the model (this should equal `to\\_col - from\\_col + 1`).
* `costs`: An array of size [to\\_col - from\\_col + 1] for the column cost coefficients.
* `lower`: An array of size [to\\_col - from\\_col + 1] for the column lower bounds.
* `upper`: An array of size [to\\_col - from\\_col + 1] for the column upper bounds.
* `num_nz`: An integer to be populated with the number of non-zero elements in the constraint matrix.
* `matrix_start`: An array of size [to\\_col - from\\_col + 1] with the start indices of each column in `matrix_index` and `matrix_value`.
* `matrix_index`: An array of size [num\\_nz] with the row indices of each element in the constraint matrix.
* `matrix_value`: An array of size [num\\_nz] with the non-zero elements of the constraint matrix.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getColsByRange(highs, from_col, to_col, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getColsByRange, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, from_col, to_col, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

"""
    Highs_getColsBySet(highs, num_set_entries, set, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)

Get data associated with multiple columns given by an array.

This function is identical to [`Highs_getColsByRange`](@ref), except for how the columns are specified.

### Parameters
* `num_set_indices`: The number of indices in `set`.
* `set`: An array of size [num\\_set\\_entries] with the column indices to get.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getColsBySet(highs, num_set_entries, set, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getColsBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, num_set_entries, set, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

"""
    Highs_getColsByMask(highs, mask, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)

Get data associated with multiple columns given by a mask.

This function is identical to [`Highs_getColsByRange`](@ref), except for how the columns are specified.

### Parameters
* `mask`: An array of length [num\\_col] containing a `1` to get the column and `0` otherwise.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getColsByMask(highs, mask, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getColsByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, mask, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

"""
    Highs_getRowsByRange(highs, from_row, to_row, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)

Get data associated with multiple adjacent rows from the model.

To query the constraint coefficients, this function should be called twice.

First, call this function with `matrix_start`, `matrix_index`, and `matrix_value` as `NULL`. This call will populate `num_nz` with the number of nonzero elements in the corresponding section of the constraint matrix.

Second, allocate new `matrix_index` and `matrix_value` arrays of length `num_nz` and call this function again to populate the new arrays with their contents.

### Parameters
* `highs`: A pointer to the Highs instance.
* `from_row`: The first row for which to query data for.
* `to_row`: The last row (inclusive) for which to query data for.
* `num_row`: An integer to be populated with the number of rows got from the smodel.
* `lower`: An array of size [to\\_row - from\\_row + 1] for the row lower bounds.
* `upper`: An array of size [to\\_row - from\\_row + 1] for the row upper bounds.
* `num_nz`: An integer to be populated with the number of non-zero elements in the constraint matrix.
* `matrix_start`: An array of size [to\\_row - from\\_row + 1] with the start indices of each row in `matrix_index` and `matrix_value`.
* `matrix_index`: An array of size [num\\_nz] with the column indices of each element in the constraint matrix.
* `matrix_value`: An array of size [num\\_nz] with the non-zero elements of the constraint matrix.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getRowsByRange(highs, from_row, to_row, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getRowsByRange, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, from_row, to_row, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

"""
    Highs_getRowsBySet(highs, num_set_entries, set, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)

Get data associated with multiple rows given by an array.

This function is identical to [`Highs_getRowsByRange`](@ref), except for how the rows are specified.

### Parameters
* `num_set_indices`: The number of indices in `set`.
* `set`: An array of size [num\\_set\\_entries] containing the row indices to get.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getRowsBySet(highs, num_set_entries, set, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getRowsBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, num_set_entries, set, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

"""
    Highs_getRowsByMask(highs, mask, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)

Get data associated with multiple rows given by a mask.

This function is identical to [`Highs_getRowsByRange`](@ref), except for how the rows are specified.

### Parameters
* `mask`: An array of length [num\\_row] containing a `1` to get the row and `0` otherwise.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getRowsByMask(highs, mask, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getRowsByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, mask, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

"""
    Highs_getRowName(highs, row, name)

Get the name of a row.

### Parameters
* `row`: The index of the row to query.
* `name`: A pointer in which to store the name of the row. This must have length `kHighsMaximumStringLength`.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getRowName(highs, row, name)
    ccall((:Highs_getRowName, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cchar}), highs, row, name)
end

"""
    Highs_getRowByName(highs, name, row)

Get the index of a row from its name.

If multiple rows have the same name, or if no row exists with `name`, this function returns `kHighsStatusError`.

### Parameters
* `name`: A pointer of the name of the row to query.
* `row`: A pointer in which to store the index of the row
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getRowByName(highs, name, row)
    ccall((:Highs_getRowByName, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{HighsInt}), highs, name, row)
end

"""
    Highs_getColName(highs, col, name)

Get the name of a column.

### Parameters
* `col`: The index of the column to query.
* `name`: A pointer in which to store the name of the column. This must have length `kHighsMaximumStringLength`.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getColName(highs, col, name)
    ccall((:Highs_getColName, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cchar}), highs, col, name)
end

"""
    Highs_getColByName(highs, name, col)

Get the index of a column from its name.

If multiple columns have the same name, or if no column exists with `name`, this function returns `kHighsStatusError`.

### Parameters
* `name`: A pointer of the name of the column to query.
* `col`: A pointer in which to store the index of the column
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getColByName(highs, name, col)
    ccall((:Highs_getColByName, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{HighsInt}), highs, name, col)
end

"""
    Highs_getColIntegrality(highs, col, integrality)

Get the integrality of a column.

### Parameters
* `col`: The index of the column to query.
* `integrality`: An integer in which the integrality of the column should be placed. The integer is one of the `kHighsVarTypeXXX` constants.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getColIntegrality(highs, col, integrality)
    ccall((:Highs_getColIntegrality, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}), highs, col, integrality)
end

"""
    Highs_deleteColsByRange(highs, from_col, to_col)

Delete multiple adjacent columns.

### Parameters
* `highs`: A pointer to the Highs instance.
* `from_col`: The index of the first column to delete.
* `to_col`: The index of the last column to delete.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_deleteColsByRange(highs, from_col, to_col)
    ccall((:Highs_deleteColsByRange, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt), highs, from_col, to_col)
end

"""
    Highs_deleteColsBySet(highs, num_set_entries, set)

Delete multiple columns given by an array of indices.

### Parameters
* `highs`: A pointer to the Highs instance.
* `num_set_entries`: The number of columns to delete.
* `set`: An array of size [num\\_set\\_entries] with the indices of the columns to delete.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_deleteColsBySet(highs, num_set_entries, set)
    ccall((:Highs_deleteColsBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}), highs, num_set_entries, set)
end

"""
    Highs_deleteColsByMask(highs, mask)

Delete multiple columns given by a mask.

### Parameters
* `highs`: A pointer to the Highs instance.
* `mask`: An array of length [num\\_col] with 1 if the column should be deleted and 0 otherwise.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_deleteColsByMask(highs, mask)
    ccall((:Highs_deleteColsByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}), highs, mask)
end

"""
    Highs_deleteRowsByRange(highs, from_row, to_row)

Delete multiple adjacent rows.

### Parameters
* `highs`: A pointer to the Highs instance.
* `from_row`: The index of the first row to delete.
* `to_row`: The index of the last row to delete.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_deleteRowsByRange(highs, from_row, to_row)
    ccall((:Highs_deleteRowsByRange, libhighs), HighsInt, (Ptr{Cvoid}, Cint, HighsInt), highs, from_row, to_row)
end

"""
    Highs_deleteRowsBySet(highs, num_set_entries, set)

Delete multiple rows given by an array of indices.

### Parameters
* `highs`: A pointer to the Highs instance.
* `num_set_entries`: The number of rows to delete.
* `set`: An array of size [num\\_set\\_entries] with the indices of the rows to delete.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_deleteRowsBySet(highs, num_set_entries, set)
    ccall((:Highs_deleteRowsBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}), highs, num_set_entries, set)
end

"""
    Highs_deleteRowsByMask(highs, mask)

Delete multiple rows given by a mask.

### Parameters
* `highs`: A pointer to the Highs instance.
* `mask`: An array of length [num\\_row] with `1` if the row should be deleted and `0` otherwise. The new index of any column not deleted is stored in place of the value `0`.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_deleteRowsByMask(highs, mask)
    ccall((:Highs_deleteRowsByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}), highs, mask)
end

"""
    Highs_scaleCol(highs, col, scaleval)

Scale a column by a constant.

Scaling a column modifies the elements in the constraint matrix, the variable bounds, and the objective coefficient.

### Parameters
* `highs`: A pointer to the Highs instance.
* `col`: The index of the column to scale.
* `scaleval`: The value by which to scale the column. If `scaleval < 0`, the variable bounds flipped.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_scaleCol(highs, col, scaleval)
    ccall((:Highs_scaleCol, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Cdouble), highs, col, scaleval)
end

"""
    Highs_scaleRow(highs, row, scaleval)

Scale a row by a constant.

### Parameters
* `highs`: A pointer to the Highs instance.
* `row`: The index of the row to scale.
* `scaleval`: The value by which to scale the row. If `scaleval < 0`, the row bounds are flipped.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_scaleRow(highs, row, scaleval)
    ccall((:Highs_scaleRow, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Cdouble), highs, row, scaleval)
end

"""
    Highs_getInfinity(highs)

Return the value of infinity used by HiGHS.

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
The value of infinity used by HiGHS.
"""
function Highs_getInfinity(highs)
    ccall((:Highs_getInfinity, libhighs), Cdouble, (Ptr{Cvoid},), highs)
end

"""
    Highs_getSizeofHighsInt(highs)

Return the size of integers used by HiGHS.

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
The size of integers used by HiGHS.
"""
function Highs_getSizeofHighsInt(highs)
    ccall((:Highs_getSizeofHighsInt, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_getNumCol(highs)

Return the number of columns in the model.

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
The number of columns in the model.
"""
function Highs_getNumCol(highs)
    ccall((:Highs_getNumCol, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_getNumRow(highs)

Return the number of rows in the model.

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
The number of rows in the model.
"""
function Highs_getNumRow(highs)
    ccall((:Highs_getNumRow, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_getNumNz(highs)

Return the number of nonzeros in the constraint matrix of the model.

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
The number of nonzeros in the constraint matrix of the model.
"""
function Highs_getNumNz(highs)
    ccall((:Highs_getNumNz, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_getHessianNumNz(highs)

Return the number of nonzeroes in the Hessian matrix of the model.

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
The number of nonzeroes in the Hessian matrix of the model.
"""
function Highs_getHessianNumNz(highs)
    ccall((:Highs_getHessianNumNz, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_getPresolvedNumCol(highs)

Return the number of columns in the presolved model.

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
The number of columns in the presolved model.
"""
function Highs_getPresolvedNumCol(highs)
    ccall((:Highs_getPresolvedNumCol, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_getPresolvedNumRow(highs)

Return the number of rows in the presolved model.

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
The number of rows in the presolved model.
"""
function Highs_getPresolvedNumRow(highs)
    ccall((:Highs_getPresolvedNumRow, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_getPresolvedNumNz(highs)

Return the number of nonzeros in the constraint matrix of the presolved model.

### Parameters
* `highs`: A pointer to the Highs instance.
### Returns
The number of nonzeros in the constraint matrix of the presolved model.
"""
function Highs_getPresolvedNumNz(highs)
    ccall((:Highs_getPresolvedNumNz, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_getModel(highs, a_format, q_format, num_col, num_row, num_nz, hessian_num_nz, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value, integrality)

Get the data from a HiGHS model.

The input arguments have the same meaning (in a different order) to those used in [`Highs_passModel`](@ref).

Note that all arrays must be pre-allocated to the correct size before calling [`Highs_getModel`](@ref). Use the following query methods to check the appropriate size: - [`Highs_getNumCol`](@ref) - [`Highs_getNumRow`](@ref) - [`Highs_getNumNz`](@ref) - [`Highs_getHessianNumNz`](@ref)

### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getModel(highs, a_format, q_format, num_col, num_row, num_nz, hessian_num_nz, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value, integrality)
    ccall((:Highs_getModel, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}), highs, a_format, q_format, num_col, num_row, num_nz, hessian_num_nz, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value, integrality)
end

"""
    Highs_getLp(highs, a_format, num_col, num_row, num_nz, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality)

Get the data from a HiGHS LP.

The input arguments have the same meaning (in a different order) to those used in [`Highs_passModel`](@ref).

Note that all arrays must be pre-allocated to the correct size before calling [`Highs_getModel`](@ref). Use the following query methods to check the appropriate size: - [`Highs_getNumCol`](@ref) - [`Highs_getNumRow`](@ref) - [`Highs_getNumNz`](@ref)

### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getLp(highs, a_format, num_col, num_row, num_nz, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality)
    ccall((:Highs_getLp, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}), highs, a_format, num_col, num_row, num_nz, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality)
end

"""
    Highs_getPresolvedLp(highs, a_format, num_col, num_row, num_nz, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality)

Get the data from a HiGHS presolved LP.

The input arguments have the same meaning (in a different order) to those used in [`Highs_passModel`](@ref).

Note that all arrays must be pre-allocated to the correct size before calling [`Highs_getModel`](@ref). Use the following query methods to check the appropriate size: - [`Highs_getPresolvedNumCol`](@ref) - [`Highs_getPresolvedNumRow`](@ref) - [`Highs_getPresolvedNumNz`](@ref)

### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getPresolvedLp(highs, a_format, num_col, num_row, num_nz, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality)
    ccall((:Highs_getPresolvedLp, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}), highs, a_format, num_col, num_row, num_nz, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality)
end

"""
    Highs_crossover(highs, num_col, num_row, col_value, col_dual, row_dual)

Set a primal (and possibly dual) solution as a starting point, then run crossover to compute a basic feasible solution.

### Parameters
* `highs`: A pointer to the Highs instance.
* `num_col`: The number of variables.
* `num_row`: The number of rows.
* `col_value`: An array of length [num\\_col] with optimal primal solution for each column.
* `col_dual`: An array of length [num\\_col] with optimal dual solution for each column. May be `NULL`, in which case no dual solution is passed.
* `row_dual`: An array of length [num\\_row] with optimal dual solution for each row. . May be `NULL`, in which case no dual solution is passed.
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_crossover(highs, num_col, num_row, col_value, col_dual, row_dual)
    ccall((:Highs_crossover, libhighs), HighsInt, (Ptr{Cvoid}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), highs, num_col, num_row, col_value, col_dual, row_dual)
end

"""
    Highs_getRanging(highs, col_cost_up_value, col_cost_up_objective, col_cost_up_in_var, col_cost_up_ou_var, col_cost_dn_value, col_cost_dn_objective, col_cost_dn_in_var, col_cost_dn_ou_var, col_bound_up_value, col_bound_up_objective, col_bound_up_in_var, col_bound_up_ou_var, col_bound_dn_value, col_bound_dn_objective, col_bound_dn_in_var, col_bound_dn_ou_var, row_bound_up_value, row_bound_up_objective, row_bound_up_in_var, row_bound_up_ou_var, row_bound_dn_value, row_bound_dn_objective, row_bound_dn_in_var, row_bound_dn_ou_var)

Compute the ranging information for all costs and bounds. For nonbasic variables the ranging information is relative to the active bound. For basic variables the ranging information relates to...

For any values that are not required, pass NULL.

### Parameters
* `highs`: A pointer to the Highs instance.
* `col_cost_up_value`: The upper range of the cost value
* `col_cost_up_objective`: The objective at the upper cost range
* `col_cost_up_in_var`: The variable entering the basis at the upper cost range
* `col_cost_up_ou_var`: The variable leaving the basis at the upper cost range
* `col_cost_dn_value`: The lower range of the cost value
* `col_cost_dn_objective`: The objective at the lower cost range
* `col_cost_dn_in_var`: The variable entering the basis at the lower cost range
* `col_cost_dn_ou_var`: The variable leaving the basis at the lower cost range
* `col_bound_up_value`: The upper range of the column bound value
* `col_bound_up_objective`: The objective at the upper column bound range
* `col_bound_up_in_var`: The variable entering the basis at the upper column bound range
* `col_bound_up_ou_var`: The variable leaving the basis at the upper column bound range
* `col_bound_dn_value`: The lower range of the column bound value
* `col_bound_dn_objective`: The objective at the lower column bound range
* `col_bound_dn_in_var`: The variable entering the basis at the lower column bound range
* `col_bound_dn_ou_var`: The variable leaving the basis at the lower column bound range
* `row_bound_up_value`: The upper range of the row bound value
* `row_bound_up_objective`: The objective at the upper row bound range
* `row_bound_up_in_var`: The variable entering the basis at the upper row bound range
* `row_bound_up_ou_var`: The variable leaving the basis at the upper row bound range
* `row_bound_dn_value`: The lower range of the row bound value
* `row_bound_dn_objective`: The objective at the lower row bound range
* `row_bound_dn_in_var`: The variable entering the basis at the lower row bound range
* `row_bound_dn_ou_var`: The variable leaving the basis at the lower row bound range
### Returns
A `kHighsStatus` constant indicating whether the call succeeded.
"""
function Highs_getRanging(highs, col_cost_up_value, col_cost_up_objective, col_cost_up_in_var, col_cost_up_ou_var, col_cost_dn_value, col_cost_dn_objective, col_cost_dn_in_var, col_cost_dn_ou_var, col_bound_up_value, col_bound_up_objective, col_bound_up_in_var, col_bound_up_ou_var, col_bound_dn_value, col_bound_dn_objective, col_bound_dn_in_var, col_bound_dn_ou_var, row_bound_up_value, row_bound_up_objective, row_bound_up_in_var, row_bound_up_ou_var, row_bound_dn_value, row_bound_dn_objective, row_bound_dn_in_var, row_bound_dn_ou_var)
    ccall((:Highs_getRanging, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, col_cost_up_value, col_cost_up_objective, col_cost_up_in_var, col_cost_up_ou_var, col_cost_dn_value, col_cost_dn_objective, col_cost_dn_in_var, col_cost_dn_ou_var, col_bound_up_value, col_bound_up_objective, col_bound_up_in_var, col_bound_up_ou_var, col_bound_dn_value, col_bound_dn_objective, col_bound_dn_in_var, col_bound_dn_ou_var, row_bound_up_value, row_bound_up_objective, row_bound_up_in_var, row_bound_up_ou_var, row_bound_dn_value, row_bound_dn_objective, row_bound_dn_in_var, row_bound_dn_ou_var)
end

"""
    Highs_resetGlobalScheduler(blocking)

Releases all resources held by the global scheduler instance.

It is not thread-safe to call this function while calling [`Highs_run`](@ref) or one of the `Highs_XXXcall` methods on any other Highs instance in any thread.

After this function has terminated, it is guaranteed that eventually all previously created scheduler threads will terminate and allocated memory will be released.

After this function has returned, the option value for the number of threads may be altered to a new value before the next call to [`Highs_run`](@ref) or one of the `Highs_XXXcall` methods.

### Parameters
* `blocking`: If the `blocking` parameter has a nonzero value, then this function will not return until all memory is freed, which might be desirable when debugging heap memory, but it requires the calling thread to wait for all scheduler threads to wake-up which is usually not necessary.
### Returns
No status is returned since the function call cannot fail. Calling this function while any Highs instance is in use on any thread is undefined behavior and may cause crashes, but cannot be detected and hence is fully in the callers responsibility.
"""
function Highs_resetGlobalScheduler(blocking)
    ccall((:Highs_resetGlobalScheduler, libhighs), Cvoid, (HighsInt,), blocking)
end

function Highs_call(num_col, num_row, num_nz, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, col_value, col_dual, row_value, row_dual, col_basis_status, row_basis_status, model_status)
    ccall((:Highs_call, libhighs), HighsInt, (HighsInt, HighsInt, HighsInt, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}), num_col, num_row, num_nz, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, col_value, col_dual, row_value, row_dual, col_basis_status, row_basis_status, model_status)
end

function Highs_runQuiet(highs)
    ccall((:Highs_runQuiet, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

function Highs_setHighsLogfile(highs, logfile)
    ccall((:Highs_setHighsLogfile, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cvoid}), highs, logfile)
end

function Highs_setHighsOutput(highs, outputfile)
    ccall((:Highs_setHighsOutput, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cvoid}), highs, outputfile)
end

function Highs_getIterationCount(highs)
    ccall((:Highs_getIterationCount, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

function Highs_getSimplexIterationCount(highs)
    ccall((:Highs_getSimplexIterationCount, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

function Highs_setHighsBoolOptionValue(highs, option, value)
    ccall((:Highs_setHighsBoolOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, HighsInt), highs, option, value)
end

function Highs_setHighsIntOptionValue(highs, option, value)
    ccall((:Highs_setHighsIntOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, HighsInt), highs, option, value)
end

function Highs_setHighsDoubleOptionValue(highs, option, value)
    ccall((:Highs_setHighsDoubleOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Cdouble), highs, option, value)
end

function Highs_setHighsStringOptionValue(highs, option, value)
    ccall((:Highs_setHighsStringOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cchar}), highs, option, value)
end

function Highs_setHighsOptionValue(highs, option, value)
    ccall((:Highs_setHighsOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cchar}), highs, option, value)
end

function Highs_getHighsBoolOptionValue(highs, option, value)
    ccall((:Highs_getHighsBoolOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{HighsInt}), highs, option, value)
end

function Highs_getHighsIntOptionValue(highs, option, value)
    ccall((:Highs_getHighsIntOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{HighsInt}), highs, option, value)
end

function Highs_getHighsDoubleOptionValue(highs, option, value)
    ccall((:Highs_getHighsDoubleOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cdouble}), highs, option, value)
end

function Highs_getHighsStringOptionValue(highs, option, value)
    ccall((:Highs_getHighsStringOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cchar}), highs, option, value)
end

function Highs_getHighsOptionType(highs, option, type)
    ccall((:Highs_getHighsOptionType, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{HighsInt}), highs, option, type)
end

function Highs_resetHighsOptions(highs)
    ccall((:Highs_resetHighsOptions, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

function Highs_getHighsIntInfoValue(highs, info, value)
    ccall((:Highs_getHighsIntInfoValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{HighsInt}), highs, info, value)
end

function Highs_getHighsDoubleInfoValue(highs, info, value)
    ccall((:Highs_getHighsDoubleInfoValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cdouble}), highs, info, value)
end

function Highs_getNumCols(highs)
    ccall((:Highs_getNumCols, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

function Highs_getNumRows(highs)
    ccall((:Highs_getNumRows, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

function Highs_getHighsInfinity(highs)
    ccall((:Highs_getHighsInfinity, libhighs), Cdouble, (Ptr{Cvoid},), highs)
end

function Highs_getHighsRunTime(highs)
    ccall((:Highs_getHighsRunTime, libhighs), Cdouble, (Ptr{Cvoid},), highs)
end

function Highs_setOptionValue(highs, option, value)
    ccall((:Highs_setOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cchar}), highs, option, value)
end

function Highs_getScaledModelStatus(highs)
    ccall((:Highs_getScaledModelStatus, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

const HighsUInt = Cuint

const CMAKE_BUILD_TYPE = "Release"

const CMAKE_INSTALL_PREFIX = "/workspace/destdir"

const HIGHS_GITHASH = "50670fd4c"

const HIGHS_COMPILATION_DATE = "1970-01-01"

const HIGHS_VERSION_MAJOR = 1

const HIGHS_VERSION_MINOR = 7

const HIGHS_VERSION_PATCH = 0

const HIGHS_DIR = "/workspace/srcdir/HiGHS"

const HIGHSINT_FORMAT = "d"

const kHighsMaximumStringLength = HighsInt(512)
const kHighsStatusError = HighsInt(-1)
const kHighsStatusOk = HighsInt(0)
const kHighsStatusWarning = HighsInt(1)
const kHighsVarTypeContinuous = HighsInt(0)
const kHighsVarTypeInteger = HighsInt(1)
const kHighsVarTypeSemiContinuous = HighsInt(2)
const kHighsVarTypeSemiInteger = HighsInt(3)
const kHighsVarTypeImplicitInteger = HighsInt(4)
const kHighsOptionTypeBool = HighsInt(0)
const kHighsOptionTypeInt = HighsInt(1)
const kHighsOptionTypeDouble = HighsInt(2)
const kHighsOptionTypeString = HighsInt(3)
const kHighsInfoTypeInt = HighsInt(1)
const kHighsInfoTypeDouble = HighsInt(2)
const kHighsObjSenseMinimize = HighsInt(1)
const kHighsObjSenseMaximize = HighsInt(-1)
const kHighsMatrixFormatColwise = HighsInt(1)
const kHighsMatrixFormatRowwise = HighsInt(2)
const kHighsHessianFormatTriangular = HighsInt(1)
const kHighsHessianFormatSquare = HighsInt(2)
const kHighsSolutionStatusNone = HighsInt(0)
const kHighsSolutionStatusInfeasible = HighsInt(1)
const kHighsSolutionStatusFeasible = HighsInt(2)
const kHighsBasisValidityInvalid = HighsInt(0)
const kHighsBasisValidityValid = HighsInt(1)
const kHighsPresolveStatusNotPresolved = HighsInt(-1)
const kHighsPresolveStatusNotReduced = HighsInt(0)
const kHighsPresolveStatusInfeasible = HighsInt(1)
const kHighsPresolveStatusUnboundedOrInfeasible = HighsInt(2)
const kHighsPresolveStatusReduced = HighsInt(3)
const kHighsPresolveStatusReducedToEmpty = HighsInt(4)
const kHighsPresolveStatusTimeout = HighsInt(5)
const kHighsPresolveStatusNullError = HighsInt(6)
const kHighsPresolveStatusOptionsError = HighsInt(7)
const kHighsModelStatusNotset = HighsInt(0)
const kHighsModelStatusLoadError = HighsInt(1)
const kHighsModelStatusModelError = HighsInt(2)
const kHighsModelStatusPresolveError = HighsInt(3)
const kHighsModelStatusSolveError = HighsInt(4)
const kHighsModelStatusPostsolveError = HighsInt(5)
const kHighsModelStatusModelEmpty = HighsInt(6)
const kHighsModelStatusOptimal = HighsInt(7)
const kHighsModelStatusInfeasible = HighsInt(8)
const kHighsModelStatusUnboundedOrInfeasible = HighsInt(9)
const kHighsModelStatusUnbounded = HighsInt(10)
const kHighsModelStatusObjectiveBound = HighsInt(11)
const kHighsModelStatusObjectiveTarget = HighsInt(12)
const kHighsModelStatusTimeLimit = HighsInt(13)
const kHighsModelStatusIterationLimit = HighsInt(14)
const kHighsModelStatusUnknown = HighsInt(15)
const kHighsModelStatusSolutionLimit = HighsInt(16)
const kHighsModelStatusInterrupt = HighsInt(17)
const kHighsBasisStatusLower = HighsInt(0)
const kHighsBasisStatusBasic = HighsInt(1)
const kHighsBasisStatusUpper = HighsInt(2)
const kHighsBasisStatusZero = HighsInt(3)
const kHighsBasisStatusNonbasic = HighsInt(4)
const kHighsCallbackLogging = HighsInt(0)
const kHighsCallbackSimplexInterrupt = HighsInt(1)
const kHighsCallbackIpmInterrupt = HighsInt(2)
const kHighsCallbackMipSolution = HighsInt(3)
const kHighsCallbackMipImprovingSolution = HighsInt(4)
const kHighsCallbackMipLogging = HighsInt(5)
const kHighsCallbackMipInterrupt = HighsInt(6)
