#!format:off

const HighsInt = Cint

"""
    Highs_lpCall(num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, col_value, col_dual, row_value, row_dual, col_basis_status, row_basis_status, model_status)

Formulate and solve a linear program using HiGHS.

### Parameters
* `num_col`: the number of columns

* `num_row`: the number of rows

* `num_nz`: the number of nonzeros in the constraint matrix

* `a_format`: the format of the constraint matrix as a `kHighsMatrixFormat` constant

* `sense`: the optimization sense as a `kHighsObjSense` constant

* `offset`: the objective constant

* `col_cost`: array of length [num\\_col] with the column costs

* `col_lower`: array of length [num\\_col] with the column lower bounds

* `col_upper`: array of length [num\\_col] with the column upper bounds

* `row_lower`: array of length [num\\_row] with the row lower bounds

* `row_upper`: array of length [num\\_row] with the row upper bounds

* `a_start`: the constraint matrix is provided to HiGHS in compressed sparse column form (if `a_format` is `kHighsMatrixFormatColwise`, otherwise compressed sparse row form). The sparse matrix consists of three arrays, `a_start`, `a_index`, and `a_value`. `a_start` is an array of length [num\\_col] containing the starting index of each column in `a_index`. If `a_format` is `kHighsMatrixFormatRowwise` the array is of length [num\\_row] corresponding to each row.

* `a_index`: array of length [num\\_nz] with indices of matrix entries

* `a_value`: array of length [num\\_nz] with values of matrix entries

* `col_value`: array of length [num\\_col], filled with the primal column solution

* `col_dual`: array of length [num\\_col], filled with the dual column solution

* `row_value`: array of length [num\\_row], filled with the primal row solution

* `row_dual`: array of length [num\\_row], filled with the dual row solution

* `col_basis_status`: array of length [num\\_col], filled with the basis status of the columns in the form of a `kHighsBasisStatus` constant

* `row_basis_status`: array of length [num\\_row], filled with the basis status of the rows in the form of a `kHighsBasisStatus` constant

* `model_status`: termination status of the model after the solve in the form of a `kHighsModelStatus` constant

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_lpCall(num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, col_value, col_dual, row_value, row_dual, col_basis_status, row_basis_status, model_status)
    ccall((:Highs_lpCall, libhighs), HighsInt, (HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}), num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, col_value, col_dual, row_value, row_dual, col_basis_status, row_basis_status, model_status)
end

"""
    Highs_mipCall(num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality, col_value, row_value, model_status)

Formulate and solve a mixed-integer linear program using HiGHS.

The signature of this method is identical to [`Highs_lpCall`](@ref), except that it has an additional `integrality` argument, and that it is missing the `col_dual`, `row_dual`, `col_basis_status` and `row_basis_status` arguments.

### Parameters
* `integrality`: array of length [num\\_col] containing a `kHighsVarType` constant for each column

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_mipCall(num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality, col_value, row_value, model_status)
    ccall((:Highs_mipCall, libhighs), HighsInt, (HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}), num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality, col_value, row_value, model_status)
end

"""
    Highs_qpCall(num_col, num_row, num_nz, q_num_nz, a_format, q_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value, col_value, col_dual, row_value, row_dual, col_basis_status, row_basis_status, model_status)

Formulate and solve a quadratic program using HiGHS.

The signature of this method is identical to [`Highs_lpCall`](@ref), except that it has additional arguments for specifying the Hessian matrix.

### Parameters
* `q_num_nz`: the number of nonzeros in the Hessian matrix

* `q_format`: the format of the Hessian matrix in the form of a `kHighsHessianStatus` constant. If q\\_num\\_nz > 0, this must be `kHighsHessianFormatTriangular`

* `q_start`: the Hessian matrix is provided in the same format as the constraint matrix, using `q_start`, `q_index`, and `q_value` in the place of `a_start`, `a_index`, and `a_value`

* `q_index`: array of length [q\\_num\\_nz] with indices of matrix entries

* `q_value`: array of length [q\\_num\\_nz] with values of matrix entries

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_qpCall(num_col, num_row, num_nz, q_num_nz, a_format, q_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value, col_value, col_dual, row_value, row_dual, col_basis_status, row_basis_status, model_status)
    ccall((:Highs_qpCall, libhighs), HighsInt, (HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}), num_col, num_row, num_nz, q_num_nz, a_format, q_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value, col_value, col_dual, row_value, row_dual, col_basis_status, row_basis_status, model_status)
end

"""
    Highs_create()

Create a Highs instance and return the reference.

Call [`Highs_destroy`](@ref) on the returned reference to clean up allocated memory.

### Returns
A pointer to the Highs instance
"""
function Highs_create()
    ccall((:Highs_create, libhighs), Ptr{Cvoid}, ())
end

"""
    Highs_destroy(highs)

Destroy the model `highs` created by [`Highs_create`](@ref) and free all corresponding memory. Future calls using `highs` are not allowed.

To empty a model without invalidating `highs`, see [`Highs_clearModel`](@ref).

### Parameters
* `highs`: a pointer to the Highs instance
"""
function Highs_destroy(highs)
    ccall((:Highs_destroy, libhighs), Cvoid, (Ptr{Cvoid},), highs)
end

"""
    Highs_readModel(highs, filename)

Read a model from `filename` into `highs`.

### Parameters
* `highs`: a pointer to the Highs instance

* `filename`: the filename to read

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_readModel(highs, filename)
    ccall((:Highs_readModel, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}), highs, filename)
end

"""
    Highs_writeModel(highs, filename)

Write the model in `highs` to `filename`.

### Parameters
* `highs`: a pointer to the Highs instance

* `filename`: the filename to write.

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_writeModel(highs, filename)
    ccall((:Highs_writeModel, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}), highs, filename)
end

"""
    Highs_clearModel(highs)

Remove all variables and constraints from the model `highs`, but do not invalidate the pointer `highs`. Future calls (for example, adding new variables and constraints) are allowed.

See [`Highs_destroy`](@ref) to clear the model and free all associated memory.

### Parameters
* `highs`: a pointer to the Highs instance

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_clearModel(highs)
    ccall((:Highs_clearModel, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_run(highs)

Optimize a model. The algorithm used by HiGHS depends on the options that have been set.

### Parameters
* `highs`: a pointer to the Highs instance

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_run(highs)
    ccall((:Highs_run, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_writeSolution(highs, filename)

Write the solution information (including dual and basis status, if available) to a file.

See also: [`Highs_writeSolutionPretty`](@ref).

### Parameters
* `highs`: a pointer to the Highs instance

* `filename`: the name of the file to write the results to

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_writeSolution(highs, filename)
    ccall((:Highs_writeSolution, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}), highs, filename)
end

"""
    Highs_writeSolutionPretty(highs, filename)

Write the solution information (including dual and basis status, if available) to a file in a human-readable format.

The method identical to [`Highs_writeSolution`](@ref), except that the printout is in a human-readiable format.

### Parameters
* `highs`: a pointer to the Highs instance

* `filename`: the name of the file to write the results to

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_writeSolutionPretty(highs, filename)
    ccall((:Highs_writeSolutionPretty, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}), highs, filename)
end

"""
    Highs_passLp(highs, num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value)

Pass a linear program (LP) to HiGHS in a single function call.

The signature of this function is identical to [`Highs_passModel`](@ref), without the arguments for passing the Hessian matrix of a quadratic program and the integrality vector.

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_passLp(highs, num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value)
    ccall((:Highs_passLp, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value)
end

"""
    Highs_passMip(highs, num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality)

Pass a mixed-integer linear program (MILP) to HiGHS in a single function call.

The signature of function is identical to [`Highs_passModel`](@ref), without the arguments for passing the Hessian matrix of a quadratic program.

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_passMip(highs, num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality)
    ccall((:Highs_passMip, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}), highs, num_col, num_row, num_nz, a_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, integrality)
end

"""
    Highs_passModel(highs, num_col, num_row, num_nz, q_num_nz, a_format, q_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value, integrality)

Pass a model to HiGHS in a single function call. This is faster than constructing the model using [`Highs_addRow`](@ref) and [`Highs_addCol`](@ref).

### Parameters
* `highs`: a pointer to the Highs instance

* `num_col`: the number of columns

* `num_row`: the number of rows

* `num_nz`: the number of elements in the constraint matrix

* `q_num_nz`: the number of elements in the Hessian matrix

* `a_format`: the format of the constraint matrix to use in th form of a `kHighsMatrixFormat` constant

* `q_format`: the format of the Hessian matrix to use in the form of a `kHighsHessianFormat` constant

* `sense`: the optimization sense in the form of a `kHighsObjSense` constant

* `offset`: the constant term in the objective function

* `col_cost`: array of length [num\\_col] with the objective coefficients

* `col_lower`: array of length [num\\_col] with the lower column bounds

* `col_upper`: array of length [num\\_col] with the upper column bounds

* `row_lower`: array of length [num\\_row] with the upper row bounds

* `row_upper`: array of length [num\\_row] with the upper row bounds

* `a_start`: the constraint matrix is provided to HiGHS in compressed sparse column form (if `a_format` is `kHighsMatrixFormatColwise`, otherwise compressed sparse row form). The sparse matrix consists of three arrays, `a_start`, `a_index`, and `a_value`. `a_start` is an array of length [num\\_col] containing the starting index of each column in `a_index`. If `a_format` is `kHighsMatrixFormatRowwise` the array is of length [num\\_row] corresponding to each row.

* `a_index`: array of length [num\\_nz] with indices of matrix entries

* `a_value`: array of length [num\\_nz] with values of matrix entries

* `q_start`: the Hessian matrix is provided in the same format as the constraint matrix, using `q_start`, `q_index`, and `q_value` in the place of `a_start`, `a_index`, and `a_value`. If the model is linear, pass NULL.

* `q_index`: array of length [q\\_num\\_nz] with indices of matrix entries. If the model is linear, pass NULL.

* `q_value`: array of length [q\\_num\\_nz] with values of matrix entries. If the model is linear, pass NULL.

* `integrality`: an array of length [num\\_col] containing a `kHighsVarType` consatnt for each column

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_passModel(highs, num_col, num_row, num_nz, q_num_nz, a_format, q_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value, integrality)
    ccall((:Highs_passModel, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}), highs, num_col, num_row, num_nz, q_num_nz, a_format, q_format, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value, integrality)
end

"""
    Highs_passHessian(highs, dim, num_nz, format, start, index, value)

Set the Hessian matrix for a quadratic objective.

### Parameters
* `highs`: a pointer to the Highs instance

* `dim`: the dimension of the Hessian matrix. Should be [num\\_col].

* `num_nz`: the number of non-zero elements in the Hessian matrix

* `format`: the format of the Hessian matrix as a `kHighsHessianFormat` constant. This must be `kHighsHessianFormatTriangular`.

* `start`: the Hessian matrix is provided to HiGHS as the upper triangular component in compressed sparse column form. The sparse matrix consists of three arrays, `start`, `index`, and `value`. `start` is an array of length [num\\_col] containing the starting index of each column in `index`.

* `index`: array of length [num\\_nz] with indices of matrix entries

* `value`: array of length [num\\_nz] with values of matrix entries

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_passHessian(highs, dim, num_nz, format, start, index, value)
    ccall((:Highs_passHessian, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, dim, num_nz, format, start, index, value)
end

"""
    Highs_setBoolOptionValue(highs, option, value)

Set a boolean-valued option.

### Parameters
* `highs`: a pointer to the Highs instance

* `option`: the name of the option

* `value`: the value of the option

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_setBoolOptionValue(highs, option, value)
    ccall((:Highs_setBoolOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, HighsInt), highs, option, value)
end

"""
    Highs_setIntOptionValue(highs, option, value)

Set an int-valued option.

### Parameters
* `highs`: a pointer to the Highs instance

* `option`: the name of the option

* `value`: the value of the option

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_setIntOptionValue(highs, option, value)
    ccall((:Highs_setIntOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, HighsInt), highs, option, value)
end

"""
    Highs_setDoubleOptionValue(highs, option, value)

Set a double-valued option.

### Parameters
* `highs`: a pointer to the Highs instance

* `option`: the name of the option

* `value`: the value of the option

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_setDoubleOptionValue(highs, option, value)
    ccall((:Highs_setDoubleOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Cdouble), highs, option, value)
end

"""
    Highs_setStringOptionValue(highs, option, value)

Set a string-valued option.

### Parameters
* `highs`: a pointer to the Highs instance

* `option`: the name of the option

* `value`: the value of the option

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_setStringOptionValue(highs, option, value)
    ccall((:Highs_setStringOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cchar}), highs, option, value)
end

"""
    Highs_getBoolOptionValue(highs, option, value)

Get a boolean-valued option.

### Parameters
* `highs`: a pointer to the Highs instance

* `option`: the name of the option

* `value`: storage for the value of the option

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getBoolOptionValue(highs, option, value)
    ccall((:Highs_getBoolOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{HighsInt}), highs, option, value)
end

"""
    Highs_getIntOptionValue(highs, option, value)

Get an int-valued option.

### Parameters
* `highs`: a pointer to the Highs instance

* `option`: the name of the option

* `value`: storage for the value of the option

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getIntOptionValue(highs, option, value)
    ccall((:Highs_getIntOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{HighsInt}), highs, option, value)
end

"""
    Highs_getDoubleOptionValue(highs, option, value)

Get a double-valued option.

### Parameters
* `highs`: a pointer to the Highs instance

* `option`: the name of the option

* `value`: storage for the value of the option

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getDoubleOptionValue(highs, option, value)
    ccall((:Highs_getDoubleOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cdouble}), highs, option, value)
end

"""
    Highs_getStringOptionValue(highs, option, value)

Get a string-valued option.

### Parameters
* `highs`: a pointer to the Highs instance

* `option`: the name of the option

* `value`: pointer to allocated memory to store the value of the option

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getStringOptionValue(highs, option, value)
    ccall((:Highs_getStringOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cchar}), highs, option, value)
end

"""
    Highs_getOptionType(highs, option, type)

Get the type expected by an option.

### Parameters
* `highs`: a pointer to the Highs instance

* `option`: the name of the option

* `type`: int in which the corresponding `kHighsOptionType` constant is stored

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getOptionType(highs, option, type)
    ccall((:Highs_getOptionType, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{HighsInt}), highs, option, type)
end

"""
    Highs_resetOptions(highs)

Reset all options to their default value.

### Parameters
* `highs`: a pointer to the Highs instance

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_resetOptions(highs)
    ccall((:Highs_resetOptions, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_writeOptions(highs, filename)

Write the current options to file.

### Parameters
* `highs`: a pointer to the Highs instance

* `filename`: the filename to write the options to

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_writeOptions(highs, filename)
    ccall((:Highs_writeOptions, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}), highs, filename)
end

"""
    Highs_writeOptionsDeviations(highs, filename)

Write the value of non-default options to file.

This is similar to [`Highs_writeOptions`](@ref), except only options with non-default value are written to `filename`.

### Parameters
* `highs`: a pointer to the Highs instance

* `filename`: the filename to write the options to

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_writeOptionsDeviations(highs, filename)
    ccall((:Highs_writeOptionsDeviations, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}), highs, filename)
end

"""
    Highs_getIntInfoValue(highs, info, value)

Get an int-valued info value.

### Parameters
* `highs`: a pointer to the Highs instance

* `info`: the name of the info item

* `value`: a reference to an integer that the result will be stored in

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getIntInfoValue(highs, info, value)
    ccall((:Highs_getIntInfoValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{HighsInt}), highs, info, value)
end

"""
    Highs_getDoubleInfoValue(highs, info, value)

Get a double-valued info value.

### Parameters
* `highs`: a pointer to the Highs instance

* `info`: the name of the info item

* `value`: a reference to an double that the result will be stored in

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getDoubleInfoValue(highs, info, value)
    ccall((:Highs_getDoubleInfoValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Cdouble}), highs, info, value)
end

"""
    Highs_getInt64InfoValue(highs, info, value)

Get an int64-valued info value.

### Parameters
* `highs`: a pointer to the Highs instance

* `info`: the name of the info item

* `value`: a reference to a int64 that the result will be stored in

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getInt64InfoValue(highs, info, value)
    ccall((:Highs_getInt64InfoValue, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cchar}, Ptr{Int64}), highs, info, value)
end

"""
    Highs_getSolution(highs, col_value, col_dual, row_value, row_dual)

Get the primal and dual solution from an optimized model.

### Parameters
* `highs`: a pointer to the Highs instance

* `col_value`: array of length [num\\_col], filled with primal column values

* `col_dual`: array of length [num\\_col], filled with dual column values

* `row_value`: array of length [num\\_row], filled with primal row values

* `row_dual`: array of length [num\\_row], filled with dual row values

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getSolution(highs, col_value, col_dual, row_value, row_dual)
    ccall((:Highs_getSolution, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), highs, col_value, col_dual, row_value, row_dual)
end

"""
    Highs_getBasis(highs, col_status, row_status)

Given a linear program with a basic feasible solution, get the column and row basis statuses.

### Parameters
* `highs`: a pointer to the Highs instance

* `col_status`: array of length [num\\_col], to be filled with the column basis statuses in the form of a `kHighsBasisStatus` constant

* `row_status`: array of length [num\\_row], to be filled with the row basis statuses in the form of a `kHighsBasisStatus` constant

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getBasis(highs, col_status, row_status)
    ccall((:Highs_getBasis, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{HighsInt}), highs, col_status, row_status)
end

"""
    Highs_getModelStatus(highs)

Return the optimization status of the model in the form of a `kHighsModelStatus` constant.

### Parameters
* `highs`: a pointer to the Highs instance

### Returns
an integer corresponding to the `kHighsModelStatus` constant
"""
function Highs_getModelStatus(highs)
    ccall((:Highs_getModelStatus, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_getDualRay(highs, has_dual_ray, dual_ray_value)

Get an unbounded dual ray that is a certificate of primal infeasibility.

### Parameters
* `highs`: a pointer to the Highs instance

* `has_dual_ray`: a pointer to an int to store 1 if the dual ray exists

* `dual_ray_value`: an array of length [num\\_row] filled with the unbounded ray

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getDualRay(highs, has_dual_ray, dual_ray_value)
    ccall((:Highs_getDualRay, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{Cdouble}), highs, has_dual_ray, dual_ray_value)
end

"""
    Highs_getPrimalRay(highs, has_primal_ray, primal_ray_value)

Get an unbounded primal ray that is a certificate of dual infeasibility.

### Parameters
* `highs`: a pointer to the Highs instance

* `has_primal_ray`: a pointer to an int to store 1 if the primal ray exists

* `primal_ray_value`: an array of length [num\\_col] filled with the unbounded ray

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getPrimalRay(highs, has_primal_ray, primal_ray_value)
    ccall((:Highs_getPrimalRay, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{Cdouble}), highs, has_primal_ray, primal_ray_value)
end

"""
    Highs_getObjectiveValue(highs)

Get the primal objective function value.

### Parameters
* `highs`: a pointer to the Highs instance

### Returns
the primal objective function value
"""
function Highs_getObjectiveValue(highs)
    ccall((:Highs_getObjectiveValue, libhighs), Cdouble, (Ptr{Cvoid},), highs)
end

"""
    Highs_getBasicVariables(highs, basic_variables)

Get the indices of the rows and columns that make up the basis matrix of a basic feasible solution.

Non-negative entries are indices of columns, and negative entries are `-row\\_index - 1`. For example, `{1, -1}` would be the second column and first row.

The order of these rows and columns is important for calls to the functions: - [`Highs_getBasisInverseRow`](@ref) - [`Highs_getBasisInverseCol`](@ref) - [`Highs_getBasisSolve`](@ref) - [`Highs_getBasisTransposeSolve`](@ref) - [`Highs_getReducedRow`](@ref) - [`Highs_getReducedColumn`](@ref)

### Parameters
* `highs`: a pointer to the Highs instance

* `basic_variables`: array of size [num\\_rows], filled with the indices of the basic variables

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getBasicVariables(highs, basic_variables)
    ccall((:Highs_getBasicVariables, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}), highs, basic_variables)
end

"""
    Highs_getBasisInverseRow(highs, row, row_vector, row_num_nz, row_index)

Get a row of the inverse basis matrix

```c++
B^{-1}
```

.

See [`Highs_getBasicVariables`](@ref) for a description of the `B` matrix.

The arrays `row_vector` and `row_index` must have an allocated length of [num\\_row]. However, check `row_num_nz` to see how many non-zero elements are actually stored.

### Parameters
* `highs`: a pointer to the Highs instance

* `row`: index of the row to compute

* `row_vector`: values of the non-zero elements

* `row_num_nz`: the number of non-zeros in the row

* `row_index`: indices of the non-zero elements

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getBasisInverseRow(highs, row, row_vector, row_num_nz, row_index)
    ccall((:Highs_getBasisInverseRow, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, row, row_vector, row_num_nz, row_index)
end

"""
    Highs_getBasisInverseCol(highs, col, col_vector, col_num_nz, col_index)

Get a column of the inverse basis matrix

```c++
B^{-1}
```

.

See [`Highs_getBasicVariables`](@ref) for a description of the `B` matrix.

The arrays `col_vector` and `col_index` must have an allocated length of [num\\_row]. However, check `col_num_nz` to see how many non-zero elements are actually stored.

### Parameters
* `highs`: a pointer to the Highs instance

* `col`: index of the column to compute

* `col_vector`: values of the non-zero elements

* `col_num_nz`: the number of non-zeros in the column

* `col_index`: indices of the non-zero elements

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getBasisInverseCol(highs, col, col_vector, col_num_nz, col_index)
    ccall((:Highs_getBasisInverseCol, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, col, col_vector, col_num_nz, col_index)
end

"""
    Highs_getBasisSolve(highs, rhs, solution_vector, solution_num_nz, solution_index)

Compute

```c++
\\mathbf{x}=B^{-1}\\mathbf{b}
```

for a given vector

```c++
\\mathbf{b}
```

.

See [`Highs_getBasicVariables`](@ref) for a description of the `B` matrix.

The arrays `solution_vector` and `solution_index` must have an allocated length of [num\\_row]. However, check `solution_num_nz` to see how many non-zero elements are actually stored.

### Parameters
* `highs`: a pointer to the Highs instance

* `rhs`: the right-hand side vector `b`

* `solution_vector`: values of the non-zero elements

* `solution_num_nz`: the number of non-zeros in the solution

* `solution_index`: indices of the non-zero elements

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getBasisSolve(highs, rhs, solution_vector, solution_num_nz, solution_index)
    ccall((:Highs_getBasisSolve, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, rhs, solution_vector, solution_num_nz, solution_index)
end

"""
    Highs_getBasisTransposeSolve(highs, rhs, solution_vector, solution_nz, solution_index)

Compute

```c++
\\mathbf{x}=B^{-T}\\mathbf{b}
```

for a given vector

```c++
\\mathbf{b}
```

.

See [`Highs_getBasicVariables`](@ref) for a description of the `B` matrix.

The arrays `solution_vector` and `solution_index` must have an allocated length of [num\\_row]. However, check `solution_num_nz` to see how many non-zero elements are actually stored.

### Parameters
* `highs`: a pointer to the Highs instance

* `rhs`: the right-hand side vector `b`

* `solution_vector`: values of the non-zero elements

* `solution_num_nz`: the number of non-zeros in the solution

* `solution_index`: indices of the non-zero elements

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getBasisTransposeSolve(highs, rhs, solution_vector, solution_nz, solution_index)
    ccall((:Highs_getBasisTransposeSolve, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, rhs, solution_vector, solution_nz, solution_index)
end

"""
    Highs_getReducedRow(highs, row, row_vector, row_num_nz, row_index)

Compute a row of

```c++
B^{-1}A
```

.

See [`Highs_getBasicVariables`](@ref) for a description of the `B` matrix.

The arrays `row_vector` and `row_index` must have an allocated length of [num\\_row]. However, check `row_num_nz` to see how many non-zero elements are actually stored.

### Parameters
* `highs`: a pointer to the Highs instance

* `row`: index of the row to compute

* `row_vector`: values of the non-zero elements

* `row_num_nz`: the number of non-zeros in the row

* `row_index`: indices of the non-zero elements

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getReducedRow(highs, row, row_vector, row_num_nz, row_index)
    ccall((:Highs_getReducedRow, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, row, row_vector, row_num_nz, row_index)
end

"""
    Highs_getReducedColumn(highs, col, col_vector, col_num_nz, col_index)

Compute a column of

```c++
B^{-1}A
```

.

See [`Highs_getBasicVariables`](@ref) for a description of the `B` matrix.

The arrays `col_vector` and `col_index` must have an allocated length of [num\\_row]. However, check `col_num_nz` to see how many non-zero elements are actually stored.

### Parameters
* `highs`: a pointer to the Highs instance

* `col`: index of the column to compute

* `col_vector`: values of the non-zero elements

* `col_num_nz`: the number of non-zeros in the column

* `col_index`: indices of the non-zero elements

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getReducedColumn(highs, col, col_vector, col_num_nz, col_index)
    ccall((:Highs_getReducedColumn, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, col, col_vector, col_num_nz, col_index)
end

"""
    Highs_setBasis(highs, col_status, row_status)

Set a basic feasible solution by passing the column and row basis statuses to the model.

### Parameters
* `highs`: a pointer to the Highs instance

* `col_status`: an array of length [num\\_col] with the column basis status in the form of `kHighsBasisStatus` constants

* `row_status`: an array of length [num\\_row] with the row basis status in the form of `kHighsBasisStatus` constants

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_setBasis(highs, col_status, row_status)
    ccall((:Highs_setBasis, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{HighsInt}), highs, col_status, row_status)
end

"""
    Highs_setLogicalBasis(highs)

Set a logical basis in the model.

### Parameters
* `highs`: a pointer to the Highs instance

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_setLogicalBasis(highs)
    ccall((:Highs_setLogicalBasis, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_getRunTime(highs)

Return the cumulative wall-clock time spent in [`Highs_run`](@ref).

### Parameters
* `highs`: a pointer to the Highs instance

### Returns
the cumulative wall-clock time spent in [`Highs_run`](@ref)
"""
function Highs_getRunTime(highs)
    ccall((:Highs_getRunTime, libhighs), Cdouble, (Ptr{Cvoid},), highs)
end

"""
    Highs_addRow(highs, lower, upper, num_new_nz, index, value)

Add a new row (a linear constraint) to the model.

### Parameters
* `highs`: a pointer to the Highs instance

* `lower`: lower bound of the row

* `upper`: upper bound of the row

* `num_new_nz`: number of non-zeros in the row

* `index`: array of size [num\\_new\\_nz] with column indices

* `value`: array of size [num\\_new\\_nz] with column values

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_addRow(highs, lower, upper, num_new_nz, index, value)
    ccall((:Highs_addRow, libhighs), HighsInt, (Ptr{Cvoid}, Cdouble, Cdouble, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}), highs, lower, upper, num_new_nz, index, value)
end

"""
    Highs_addRows(highs, num_new_row, lower, upper, num_new_nz, starts, index, value)

Add multiple rows (linear constraints) to the model.

### Parameters
* `highs`: a pointer to the Highs instance

* `num_new_row`: the number of new rows to add

* `lower`: array of size [num\\_new\\_row] with the lower bounds of the rows

* `upper`: array of size [num\\_new\\_row] with the upper bounds of the rows

* `num_new_nz`: number of non-zeros in the rows

* `starts`: the constraint coefficients are given as a matrix in compressed sparse row form by the arrays `starts`, `index`, and `value`. `starts` is an array of size [num\\_new\\_rows] with the start index of each row in indices and values.

* `index`: array of size [num\\_new\\_nz] with column indices

* `value`: array of size [num\\_new\\_nz] with column values

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_addRows(highs, num_new_row, lower, upper, num_new_nz, starts, index, value)
    ccall((:Highs_addRows, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{Cdouble}, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, num_new_row, lower, upper, num_new_nz, starts, index, value)
end

"""
    Highs_addCol(highs, cost, lower, upper, num_new_nz, index, value)

Add a new column (variable) to the model.

### Parameters
* `highs`: a pointer to the Highs instance

* `cost`: objective coefficient of the column

* `lower`: lower bound of the column

* `upper`: upper bound of the column

* `num_new_nz`: number of non-zeros in the column

* `index`: array of size [num\\_new\\_nz] with the row indices

* `value`: array of size [num\\_new\\_nz] with row values

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_addCol(highs, cost, lower, upper, num_new_nz, index, value)
    ccall((:Highs_addCol, libhighs), HighsInt, (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}), highs, cost, lower, upper, num_new_nz, index, value)
end

"""
    Highs_addCols(highs, num_new_col, costs, lower, upper, num_new_nz, starts, index, value)

Add multiple columns (linear constraints) to the model.

### Parameters
* `highs`: a pointer to the Highs instance

* `num_new_col`: number of new columns to add

* `costs`: array of size [num\\_new\\_col] with objective coefficients

* `lower`: array of size [num\\_new\\_col] with lower bounds

* `upper`: array of size [num\\_new\\_col] with upper bounds

* `num_new_nz`: number of new nonzeros in the constraint matrix

* `starts`: the constraint coefficients are given as a matrix in compressed sparse column form by the arrays `starts`, `index`, and `value`. `starts` is an array of size [num\\_new\\_cols] with the start index of each row in indices and values.

* `index`: array of size [num\\_new\\_nz] with row indices

* `value`: array of size [num\\_new\\_nz] with row values

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_addCols(highs, num_new_col, costs, lower, upper, num_new_nz, starts, index, value)
    ccall((:Highs_addCols, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, num_new_col, costs, lower, upper, num_new_nz, starts, index, value)
end

"""
    Highs_changeObjectiveSense(highs, sense)

Change the objective sense of the model.

### Parameters
* `highs`: a pointer to the Highs instance

* `sense`: the new optimization sense in the form of a `kHighsObjSense` constant

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_changeObjectiveSense(highs, sense)
    ccall((:Highs_changeObjectiveSense, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt), highs, sense)
end

"""
    Highs_changeObjectiveOffset(highs, offset)

Change the objective offset of the model.

### Parameters
* `highs`: a pointer to the Highs instance

* `offset`: the new objective offset

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_changeObjectiveOffset(highs, offset)
    ccall((:Highs_changeObjectiveOffset, libhighs), HighsInt, (Ptr{Cvoid}, Cdouble), highs, offset)
end

"""
    Highs_changeColIntegrality(highs, col, integrality)

Change the integrality of a column.

### Parameters
* `highs`: a pointer to the Highs instance

* `col`: the column index to change

* `integrality`: the new integrality of the column in the form of a `kHighsVarType` constant

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_changeColIntegrality(highs, col, integrality)
    ccall((:Highs_changeColIntegrality, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt), highs, col, integrality)
end

"""
    Highs_changeColsIntegralityByRange(highs, from_col, to_col, integrality)

Change the integrality of multiple adjacent columns.

### Parameters
* `highs`: a pointer to the Highs instance

* `from_col`: the index of the first column whose integrality changes

* `to_col`: the index of the last column whose integrality changes

* `integrality`: an array of length [to\\_col - from\\_col + 1] with the new integralities of the columns in the form of `kHighsVarType` constants

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_changeColsIntegralityByRange(highs, from_col, to_col, integrality)
    ccall((:Highs_changeColsIntegralityByRange, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Ptr{HighsInt}), highs, from_col, to_col, integrality)
end

"""
    Highs_changeColsIntegralityBySet(highs, num_set_entries, set, integrality)

Change the integrality of multiple columns given by an array of indices.

### Parameters
* `highs`: a pointer to the Highs instance

* `num_set_entries`: the number of columns to change

* `set`: an array of size [num\\_set\\_entries] with the indices of the columns to change

* `integrality`: an array of length [num\\_set\\_entries] with the new integralities of the columns in the form of `kHighsVarType` constants

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_changeColsIntegralityBySet(highs, num_set_entries, set, integrality)
    ccall((:Highs_changeColsIntegralityBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}), highs, num_set_entries, set, integrality)
end

"""
    Highs_changeColsIntegralityByMask(highs, mask, integrality)

Change the integrality of multiple columns given by a mask.

### Parameters
* `highs`: a pointer to the Highs instance

* `mask`: an array of length [num\\_col] with 1 if the column integrality should be changed and 0 otherwise

* `integrality`: an array of length [num\\_col] with the new integralities of the columns in the form of `kHighsVarType` constants

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_changeColsIntegralityByMask(highs, mask, integrality)
    ccall((:Highs_changeColsIntegralityByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{HighsInt}), highs, mask, integrality)
end

"""
    Highs_changeColCost(highs, col, cost)

Change the objective coefficient of a column.

### Parameters
* `highs`: a pointer to the Highs instance

* `col`: the index of the column fo change

* `cost`: the new objective coefficient

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_changeColCost(highs, col, cost)
    ccall((:Highs_changeColCost, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Cdouble), highs, col, cost)
end

"""
    Highs_changeColsCostByRange(highs, from_col, to_col, cost)

Change the cost coefficients of multiple adjacent columns.

### Parameters
* `highs`: a pointer to the Highs instance

* `from_col`: the index of the first column whose cost changes

* `to_col`: the index of the last column whose cost changes

* `cost`: an array of length [to\\_col - from\\_col + 1] with the new objective coefficients

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_changeColsCostByRange(highs, from_col, to_col, cost)
    ccall((:Highs_changeColsCostByRange, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Ptr{Cdouble}), highs, from_col, to_col, cost)
end

"""
    Highs_changeColsCostBySet(highs, num_set_entries, set, cost)

Change the cost of multiple columns given by an array of indices.

### Parameters
* `highs`: a pointer to the Highs instance

* `num_set_entries`: the number of columns to change

* `set`: an array of size [num\\_set\\_entries] with the indices of the columns to change

* `cost`: an array of length [num\\_set\\_entries] with the new costs of the columns.

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_changeColsCostBySet(highs, num_set_entries, set, cost)
    ccall((:Highs_changeColsCostBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}), highs, num_set_entries, set, cost)
end

"""
    Highs_changeColsCostByMask(highs, mask, cost)

Change the cost of multiple columns given by a mask.

### Parameters
* `highs`: a pointer to the Highs instance

* `mask`: an array of length [num\\_col] with 1 if the column cost should be changed and 0 otherwise

* `cost`: an array of length [num\\_col] with the new costs

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_changeColsCostByMask(highs, mask, cost)
    ccall((:Highs_changeColsCostByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{Cdouble}), highs, mask, cost)
end

"""
    Highs_changeColBounds(highs, col, lower, upper)

Change the variable bounds of a column.

### Parameters
* `highs`: a pointer to the Highs instance

* `col`: the index of the column whose bounds are to change

* `lower`: the new lower bound

* `upper`: the new upper bound

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_changeColBounds(highs, col, lower, upper)
    ccall((:Highs_changeColBounds, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Cdouble, Cdouble), highs, col, lower, upper)
end

"""
    Highs_changeColsBoundsByRange(highs, from_col, to_col, lower, upper)

Change the variable bounds of multiple adjacent columns.

### Parameters
* `highs`: a pointer to the Highs instance

* `from_col`: the index of the first column whose bound changes

* `to_col`: the index of the last column whose bound changes

* `lower`: an array of length [to\\_col - from\\_col + 1] with the new lower bounds

* `upper`: an array of length [to\\_col - from\\_col + 1] with the new upper bounds

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_changeColsBoundsByRange(highs, from_col, to_col, lower, upper)
    ccall((:Highs_changeColsBoundsByRange, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Ptr{Cdouble}, Ptr{Cdouble}), highs, from_col, to_col, lower, upper)
end

"""
    Highs_changeColsBoundsBySet(highs, num_set_entries, set, lower, upper)

Change the bounds of multiple columns given by an array of indices.

### Parameters
* `highs`: a pointer to the Highs instance

* `num_set_entries`: the number of columns to change

* `set`: an array of size [num\\_set\\_entries] with the indices of the columns to change

* `lower`: an array of length [num\\_set\\_entries] with the new lower bounds

* `upper`: an array of length [num\\_set\\_entries] with the new upper bounds

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_changeColsBoundsBySet(highs, num_set_entries, set, lower, upper)
    ccall((:Highs_changeColsBoundsBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}), highs, num_set_entries, set, lower, upper)
end

"""
    Highs_changeColsBoundsByMask(highs, mask, lower, upper)

Change the variable bounds of multiple columns given by a mask.

### Parameters
* `highs`: a pointer to the Highs instance

* `mask`: an array of length [num\\_col] with 1 if the column bounds should be changed and 0 otherwise

* `lower`: an array of length [num\\_col] with the new lower bounds

* `upper`: an array of length [num\\_col] with the new upper bounds

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_changeColsBoundsByMask(highs, mask, lower, upper)
    ccall((:Highs_changeColsBoundsByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}), highs, mask, lower, upper)
end

"""
    Highs_changeRowBounds(highs, row, lower, upper)

Change the bounds of a row.

### Parameters
* `highs`: a pointer to the Highs instance

* `row`: the index of the row whose bounds are to change

* `lower`: the new lower bound

* `upper`: the new upper bound

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_changeRowBounds(highs, row, lower, upper)
    ccall((:Highs_changeRowBounds, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Cdouble, Cdouble), highs, row, lower, upper)
end

"""
    Highs_changeRowsBoundsBySet(highs, num_set_entries, set, lower, upper)

Change the bounds of multiple rows given by an array of indices.

### Parameters
* `highs`: a pointer to the Highs instance

* `num_set_entries`: the number of rows to change

* `set`: an array of size [num\\_set\\_entries] with the indices of the rows to change

* `lower`: an array of length [num\\_set\\_entries] with the new lower bounds

* `upper`: an array of length [num\\_set\\_entries] with the new upper bounds

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_changeRowsBoundsBySet(highs, num_set_entries, set, lower, upper)
    ccall((:Highs_changeRowsBoundsBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}), highs, num_set_entries, set, lower, upper)
end

"""
    Highs_changeRowsBoundsByMask(highs, mask, lower, upper)

Change the bounds of multiple rows given by a mask.

### Parameters
* `highs`: a pointer to the Highs instance

* `mask`: an array of length [num\\_row] with 1 if the row bounds should be changed and 0 otherwise

* `lower`: an array of length [num\\_row] with the new lower bounds

* `upper`: an array of length [num\\_row] with the new upper bounds

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_changeRowsBoundsByMask(highs, mask, lower, upper)
    ccall((:Highs_changeRowsBoundsByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}), highs, mask, lower, upper)
end

"""
    Highs_changeCoeff(highs, row, col, value)

Change a coefficient in the constraint matrix.

### Parameters
* `highs`: a pointer to the Highs instance

* `row`: the index of the row to change

* `col`: the index of the col to change

* `value`: the new constraint coefficient

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_changeCoeff(highs, row, col, value)
    ccall((:Highs_changeCoeff, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Cdouble), highs, row, col, value)
end

"""
    Highs_getObjectiveSense(highs, sense)

Get the objective sense.

### Parameters
* `highs`: a pointer to the Highs instance

* `sense`: stores the current objective sense as a `kHighsObjSense` constant

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getObjectiveSense(highs, sense)
    ccall((:Highs_getObjectiveSense, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}), highs, sense)
end

"""
    Highs_getObjectiveOffset(highs, offset)

Get the objective offset.

### Parameters
* `highs`: a pointer to the Highs instance

* `offset`: stores the current objective offset

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getObjectiveOffset(highs, offset)
    ccall((:Highs_getObjectiveOffset, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cdouble}), highs, offset)
end

"""
    Highs_getColsByRange(highs, from_col, to_col, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)

Get data associated with multiple adjacent columns from the model.

To query the constraint coefficients, this function should be called twice: - First, call this function with `matrix_start`, `matrix_index`, and `matrix_value` as `NULL`. This call will populate `num_nz` with the number of nonzero elements in the corresponding section of the constraint matrix. - Second, allocate new `matrix_index` and `matrix_value` arrays of length `num_nz` and call this function again to populate the new arrays with their contents.

### Parameters
* `highs`: a pointer to the Highs instance

* `from_col`: the first column for which to query data for

* `to_col`: the last column (inclusive) for which to query data for

* `num_col`: an integer populated with the number of columns got from the model (this should equal `to\\_col - from\\_col + 1`)

* `costs`: array of size [to\\_col - from\\_col + 1] for the column cost coefficients

* `lower`: array of size [to\\_col - from\\_col + 1] for the column lower bounds

* `upper`: array of size [to\\_col - from\\_col + 1] for the column upper bounds

* `num_nz`: an integer populated with the number of non-zero elements in the constraint matrix

* `matrix_start`: array of size [to\\_col - from\\_col + 1] with the start indices of each column in `matrix_index` and `matrix_value`

* `matrix_index`: array of size [num\\_nz] with the row indices of each element in the constraint matrix

* `matrix_value`: array of size [num\\_nz] with the non-zero elements of the constraint matrix.

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getColsByRange(highs, from_col, to_col, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getColsByRange, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, from_col, to_col, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

"""
    Highs_getColsBySet(highs, num_set_entries, set, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)

Get data associated with multiple columns given by an array.

This function is identical to [`Highs_getColsByRange`](@ref), except for how the columns are specified.

### Parameters
* `num_set_indices`: the number of indices in the set

* `set`: array of size [num\\_set\\_entries] with the column indices to get

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getColsBySet(highs, num_set_entries, set, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getColsBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, num_set_entries, set, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

"""
    Highs_getColsByMask(highs, mask, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)

Get data associated with multiple columns given by a mask.

This function is identical to [`Highs_getColsByRange`](@ref), except for how the columns are specified.

### Parameters
* `mask`: array of length [num\\_col] containing a 1 to get the column and 0 otherwise

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getColsByMask(highs, mask, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getColsByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, mask, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

"""
    Highs_getRowsByRange(highs, from_row, to_row, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)

Get data associated with multiple adjacent rows from the model.

To query the constraint coefficients, this function should be called twice: - First, call this function with `matrix_start`, `matrix_index`, and `matrix_value` as `NULL`. This call will populate `num_nz` with the number of nonzero elements in the corresponding section of the constraint matrix. - Second, allocate new `matrix_index` and `matrix_value` arrays of length `num_nz` and call this function again to populate the new arrays with their contents.

### Parameters
* `highs`: a pointer to the Highs instance

* `from_row`: the first row for which to query data for

* `to_row`: the last row (inclusive) for which to query data for

* `num_row`: an integer populated with the number of row got from the model

* `lower`: array of size [to\\_row - from\\_row + 1] for the row lower bounds

* `upper`: array of size [to\\_row - from\\_row + 1] for the row upper bounds

* `num_nz`: an integer populated with the number of non-zero elements in the constraint matrix

* `matrix_start`: array of size [to\\_row - from\\_row + 1] with the start indices of each row in `matrix_index` and `matrix_value`

* `matrix_index`: array of size [num\\_nz] with the column indices of each element in the constraint matrix

* `matrix_value`: array of size [num\\_nz] with the non-zero elements of the constraint matrix.

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getRowsByRange(highs, from_row, to_row, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getRowsByRange, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, from_row, to_row, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

"""
    Highs_getRowsBySet(highs, num_set_entries, set, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)

Get data associated with multiple rows given by an array.

This function is identical to [`Highs_getRowsByRange`](@ref), except for how the rows are specified.

### Parameters
* `num_set_indices`: the number of indices in the set

* `set`: array of size [num\\_set\\_entries] with the row indices to get

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getRowsBySet(highs, num_set_entries, set, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getRowsBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, num_set_entries, set, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

"""
    Highs_getRowsByMask(highs, mask, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)

Get data associated with multiple rows given by a mask.

This function is identical to [`Highs_getRowsByRange`](@ref), except for how the rows are specified.

### Parameters
* `mask`: array of length [num\\_row] containing a 1 to get the row and 0 otherwise

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getRowsByMask(highs, mask, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getRowsByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, mask, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

"""
    Highs_deleteColsByRange(highs, from_col, to_col)

Delete multiple adjacent columns.

### Parameters
* `highs`: a pointer to the Highs instance

* `from_col`: the index of the first column to delete

* `to_col`: the index of the last column to delete

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_deleteColsByRange(highs, from_col, to_col)
    ccall((:Highs_deleteColsByRange, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt), highs, from_col, to_col)
end

"""
    Highs_deleteColsBySet(highs, num_set_entries, set)

Delete multiple columns given by an array of indices.

### Parameters
* `highs`: a pointer to the Highs instance

* `num_set_entries`: the number of columns to delete

* `set`: an array of size [num\\_set\\_entries] with the indices of the columns to delete

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_deleteColsBySet(highs, num_set_entries, set)
    ccall((:Highs_deleteColsBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}), highs, num_set_entries, set)
end

"""
    Highs_deleteColsByMask(highs, mask)

Delete multiple columns given by a mask.

### Parameters
* `highs`: a pointer to the Highs instance

* `mask`: an array of length [num\\_col] with 1 if the column should be deleted and 0 otherwise

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_deleteColsByMask(highs, mask)
    ccall((:Highs_deleteColsByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}), highs, mask)
end

"""
    Highs_deleteRowsByRange(highs, from_row, to_row)

Delete multiple adjacent rows.

### Parameters
* `highs`: a pointer to the Highs instance

* `from_row`: the index of the first row to delete

* `to_row`: the index of the last row to delete

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_deleteRowsByRange(highs, from_row, to_row)
    ccall((:Highs_deleteRowsByRange, libhighs), HighsInt, (Ptr{Cvoid}, Cint, HighsInt), highs, from_row, to_row)
end

"""
    Highs_deleteRowsBySet(highs, num_set_entries, set)

Delete multiple rows given by an array of indices.

### Parameters
* `highs`: a pointer to the Highs instance

* `num_set_entries`: the number of rows to delete

* `set`: an array of size [num\\_set\\_entries] with the indices of the rows to delete

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_deleteRowsBySet(highs, num_set_entries, set)
    ccall((:Highs_deleteRowsBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}), highs, num_set_entries, set)
end

"""
    Highs_deleteRowsByMask(highs, mask)

Delete multiple rows given by a mask.

### Parameters
* `highs`: a pointer to the Highs instance

* `mask`: an array of length [num\\_row] with 1 if the row should be deleted and 0 otherwise. New index of any column not deleted is returned in place of the value 0.

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_deleteRowsByMask(highs, mask)
    ccall((:Highs_deleteRowsByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}), highs, mask)
end

"""
    Highs_scaleCol(highs, col, scaleval)

Scale a column by a constant.

Scaling a column modifies the elements in the constraint matrix, the variable bounds, and the objective coefficient.

If scaleval < 0, the variable bounds flipped.

### Parameters
* `highs`: a pointer to the Highs instance

* `col`: the index of the column to scale

* `scaleval`: the value by which to scale the column

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_scaleCol(highs, col, scaleval)
    ccall((:Highs_scaleCol, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Cdouble), highs, col, scaleval)
end

"""
    Highs_scaleRow(highs, row, scaleval)

Scale a row by a constant.

If scaleval < 0, the row bounds are flipped.

### Parameters
* `highs`: a pointer to the Highs instance

* `row`: the index of the row to scale

* `scaleval`: the value by which to scale the row

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_scaleRow(highs, row, scaleval)
    ccall((:Highs_scaleRow, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Cdouble), highs, row, scaleval)
end

"""
    Highs_getInfinity(highs)

Return the value of infinity used by HiGHS.

### Parameters
* `highs`: a pointer to the Highs instance

### Returns
the value of infinity used by HiGHS
"""
function Highs_getInfinity(highs)
    ccall((:Highs_getInfinity, libhighs), Cdouble, (Ptr{Cvoid},), highs)
end

"""
    Highs_getNumCol(highs)

Return the number of columns in the model.

### Parameters
* `highs`: a pointer to the Highs instance

### Returns
the number of columns in the model
"""
function Highs_getNumCol(highs)
    ccall((:Highs_getNumCol, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_getNumRow(highs)

Return the number of rows in the model.

### Parameters
* `highs`: a pointer to the Highs instance

### Returns
the number of rows in the model.
"""
function Highs_getNumRow(highs)
    ccall((:Highs_getNumRow, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_getNumNz(highs)

Return the number of nonzeros in the constraint matrix of the model.

### Parameters
* `highs`: a pointer to the Highs instance

### Returns
the number of nonzeros in the constraint matrix of the model.
"""
function Highs_getNumNz(highs)
    ccall((:Highs_getNumNz, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_getHessianNumNz(highs)

Return the number of nonzeroes in the Hessian matrix of the model.

### Parameters
* `highs`: a pointer to the Highs instance

### Returns
the number of nonzeroes in the Hessian matrix of the model.
"""
function Highs_getHessianNumNz(highs)
    ccall((:Highs_getHessianNumNz, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_getModel(highs, a_format, q_format, num_col, num_row, num_nz, hessian_num_nz, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value, integrality)

Get the data from a HiGHS model.

The input arguments have the same meaning (in a different order) to those used in [`Highs_passModel`](@ref).

Note that all arrays must be pre-allocated to the correct size before calling [`Highs_getModel`](@ref). Use the following query methods to check the appropriate size: - [`Highs_getNumCol`](@ref) - [`Highs_getNumRow`](@ref) - [`Highs_getNumNz`](@ref) - [`Highs_getHessianNumNz`](@ref)

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_getModel(highs, a_format, q_format, num_col, num_row, num_nz, hessian_num_nz, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value, integrality)
    ccall((:Highs_getModel, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}), highs, a_format, q_format, num_col, num_row, num_nz, hessian_num_nz, sense, offset, col_cost, col_lower, col_upper, row_lower, row_upper, a_start, a_index, a_value, q_start, q_index, q_value, integrality)
end

"""
    Highs_crossover(highs)

Given a model solved with an interior point method, run crossover to compute a basic feasible solution.

### Parameters
* `highs`: a pointer to the Highs instance

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_crossover(highs)
    ccall((:Highs_crossover, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

"""
    Highs_crossover_set(highs, num_col, num_row, col_value, col_dual, row_dual)

Set a primal (and possibly dual) solution as a starting point, then run crossover to compute a basic feasible solution. If there is no dual solution, pass col\\_dual and row\\_dual as nullptr.

### Parameters
* `highs`: a pointer to the Highs instance

* `num_col`: the number of variables

* `num_row`: the number of rows

* `col_value`: array of length [num\\_col] with optimal primal solution for each column

* `col_dual`: array of length [num\\_col] with optimal dual solution for each column

* `row_dual`: array of length [num\\_row] with optimal dual solution for each row

### Returns
a `kHighsStatus` constant indicating whether the call succeeded
"""
function Highs_crossover_set(highs, num_col, num_row, col_value, col_dual, row_dual)
    ccall((:Highs_crossover_set, libhighs), HighsInt, (Ptr{Cvoid}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), highs, num_col, num_row, col_value, col_dual, row_dual)
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

const HIGHS_GITHASH = ""

const HIGHS_COMPILATION_DATE = ""

const HIGHS_VERSION_MAJOR = 1

const HIGHS_VERSION_MINOR = 2

const HIGHS_VERSION_PATCH = 1

const HIGHS_DIR = "/workspace/srcdir/HiGHS"

const HIGHSINT_FORMAT = "d"

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
const kHighsBasisStatusLower = HighsInt(0)
const kHighsBasisStatusBasic = HighsInt(1)
const kHighsBasisStatusUpper = HighsInt(2)
const kHighsBasisStatusZero = HighsInt(3)
const kHighsBasisStatusNonbasic = HighsInt(4)
