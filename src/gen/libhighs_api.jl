#! format: off

# Julia wrapper for header: highs_c_api.h
# Automatically generated using Clang.jl


function Highs_lpCall(numcol, numrow, numnz, a_format, sense, offset, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue, colvalue, coldual, rowvalue, rowdual, colbasisstatus, rowbasisstatus, modelstatus)
    ccall((:Highs_lpCall, libhighs), HighsInt, (HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}), numcol, numrow, numnz, a_format, sense, offset, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue, colvalue, coldual, rowvalue, rowdual, colbasisstatus, rowbasisstatus, modelstatus)
end

function Highs_mipCall(numcol, numrow, numnz, a_format, sense, offset, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue, integrality, colvalue, rowvalue, modelstatus)
    ccall((:Highs_mipCall, libhighs), HighsInt, (HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}), numcol, numrow, numnz, a_format, sense, offset, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue, integrality, colvalue, rowvalue, modelstatus)
end

function Highs_qpCall(numcol, numrow, numnz, q_numnz, a_format, q_format, sense, offset, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue, qstart, qindex, qvalue, colvalue, coldual, rowvalue, rowdual, colbasisstatus, rowbasisstatus, modelstatus)
    ccall((:Highs_qpCall, libhighs), HighsInt, (HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}), numcol, numrow, numnz, q_numnz, a_format, q_format, sense, offset, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue, qstart, qindex, qvalue, colvalue, coldual, rowvalue, rowdual, colbasisstatus, rowbasisstatus, modelstatus)
end

function Highs_lpDimMpsRead(numcol, numrow, numNz)
    ccall((:Highs_lpDimMpsRead, libhighs), HighsInt, (Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}), numcol, numrow, numNz)
end

function Highs_lpDataMpsRead(numcol, numrow, sense, offset, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue)
    ccall((:Highs_lpDataMpsRead, libhighs), HighsInt, (HighsInt, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), numcol, numrow, sense, offset, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue)
end

function Highs_create()
    ccall((:Highs_create, libhighs), Ptr{Cvoid}, ())
end

function Highs_destroy(highs)
    ccall((:Highs_destroy, libhighs), Cvoid, (Ptr{Cvoid},), highs)
end

function Highs_readModel(highs, filename)
    ccall((:Highs_readModel, libhighs), HighsInt, (Ptr{Cvoid}, Cstring), highs, filename)
end

function Highs_writeModel(highs, filename)
    ccall((:Highs_writeModel, libhighs), HighsInt, (Ptr{Cvoid}, Cstring), highs, filename)
end

function Highs_clearModel(highs)
    ccall((:Highs_clearModel, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

function Highs_run(highs)
    ccall((:Highs_run, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

function Highs_writeSolution(highs, filename)
    ccall((:Highs_writeSolution, libhighs), HighsInt, (Ptr{Cvoid}, Cstring), highs, filename)
end

function Highs_writeSolutionPretty(highs, filename)
    ccall((:Highs_writeSolutionPretty, libhighs), HighsInt, (Ptr{Cvoid}, Cstring), highs, filename)
end

function Highs_passLp(highs, numcol, numrow, numnz, a_format, sense, offset, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue)
    ccall((:Highs_passLp, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, numcol, numrow, numnz, a_format, sense, offset, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue)
end

function Highs_passMip(highs, numcol, numrow, numnz, a_format, sense, offset, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue, integrality)
    ccall((:Highs_passMip, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}), highs, numcol, numrow, numnz, a_format, sense, offset, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue, integrality)
end

function Highs_passModel(highs, numcol, numrow, numnz, q_num_nz, a_format, q_format, sense, offset, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue, qstart, qindex, qvalue, integrality)
    ccall((:Highs_passModel, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, HighsInt, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}), highs, numcol, numrow, numnz, q_num_nz, a_format, q_format, sense, offset, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue, qstart, qindex, qvalue, integrality)
end

function Highs_passHessian(highs, dim, num_nz, format, start, index, value)
    ccall((:Highs_passHessian, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, dim, num_nz, format, start, index, value)
end

function Highs_setBoolOptionValue(highs, option, value)
    ccall((:Highs_setBoolOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, HighsInt), highs, option, value)
end

function Highs_setIntOptionValue(highs, option, value)
    ccall((:Highs_setIntOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, HighsInt), highs, option, value)
end

function Highs_setDoubleOptionValue(highs, option, value)
    ccall((:Highs_setDoubleOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Cdouble), highs, option, value)
end

function Highs_setStringOptionValue(highs, option, value)
    ccall((:Highs_setStringOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Cstring), highs, option, value)
end

function Highs_setOptionValue(highs, option, value)
    ccall((:Highs_setOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Cstring), highs, option, value)
end

function Highs_getBoolOptionValue(highs, option, value)
    ccall((:Highs_getBoolOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Ptr{HighsInt}), highs, option, value)
end

function Highs_getIntOptionValue(highs, option, value)
    ccall((:Highs_getIntOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Ptr{HighsInt}), highs, option, value)
end

function Highs_getDoubleOptionValue(highs, option, value)
    ccall((:Highs_getDoubleOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Ptr{Cdouble}), highs, option, value)
end

function Highs_getStringOptionValue(highs, option, value)
    ccall((:Highs_getStringOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Cstring), highs, option, value)
end

function Highs_getOptionType(highs, option, type)
    ccall((:Highs_getOptionType, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Ptr{HighsInt}), highs, option, type)
end

function Highs_resetOptions(highs)
    ccall((:Highs_resetOptions, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

function Highs_writeOptions(highs, filename)
    ccall((:Highs_writeOptions, libhighs), HighsInt, (Ptr{Cvoid}, Cstring), highs, filename)
end

function Highs_writeOptionsDeviations(highs, filename)
    ccall((:Highs_writeOptionsDeviations, libhighs), HighsInt, (Ptr{Cvoid}, Cstring), highs, filename)
end

function Highs_getIntInfoValue(highs, info, value)
    ccall((:Highs_getIntInfoValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Ptr{HighsInt}), highs, info, value)
end

function Highs_getDoubleInfoValue(highs, info, value)
    ccall((:Highs_getDoubleInfoValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Ptr{Cdouble}), highs, info, value)
end

function Highs_getSolution(highs, colvalue, coldual, rowvalue, rowdual)
    ccall((:Highs_getSolution, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), highs, colvalue, coldual, rowvalue, rowdual)
end

function Highs_getBasis(highs, colstatus, rowstatus)
    ccall((:Highs_getBasis, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{HighsInt}), highs, colstatus, rowstatus)
end

function Highs_getModelStatus(highs)
    ccall((:Highs_getModelStatus, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

function Highs_getScaledModelStatus(highs)
    ccall((:Highs_getScaledModelStatus, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

function Highs_getDualRay(highs, has_dual_ray, dual_ray_value)
    ccall((:Highs_getDualRay, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{Cdouble}), highs, has_dual_ray, dual_ray_value)
end

function Highs_getPrimalRay(highs, has_primal_ray, primal_ray_value)
    ccall((:Highs_getPrimalRay, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{Cdouble}), highs, has_primal_ray, primal_ray_value)
end

function Highs_getObjectiveValue(highs)
    ccall((:Highs_getObjectiveValue, libhighs), Cdouble, (Ptr{Cvoid},), highs)
end

function Highs_getBasicVariables(highs, basic_variables)
    ccall((:Highs_getBasicVariables, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}), highs, basic_variables)
end

function Highs_getBasisInverseRow(highs, row, row_vector, row_num_nz, row_indices)
    ccall((:Highs_getBasisInverseRow, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, row, row_vector, row_num_nz, row_indices)
end

function Highs_getBasisInverseCol(highs, col, col_vector, col_num_nz, col_indices)
    ccall((:Highs_getBasisInverseCol, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, col, col_vector, col_num_nz, col_indices)
end

function Highs_getBasisSolve(highs, rhs, solution_vector, solution_num_nz, solution_indices)
    ccall((:Highs_getBasisSolve, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, rhs, solution_vector, solution_num_nz, solution_indices)
end

function Highs_getBasisTransposeSolve(highs, rhs, solution_vector, solution_nz, solution_indices)
    ccall((:Highs_getBasisTransposeSolve, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, rhs, solution_vector, solution_nz, solution_indices)
end

function Highs_getReducedRow(highs, row, row_vector, row_num_nz, row_indices)
    ccall((:Highs_getReducedRow, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, row, row_vector, row_num_nz, row_indices)
end

function Highs_getReducedColumn(highs, col, col_vector, col_num_nz, col_indices)
    ccall((:Highs_getReducedColumn, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}), highs, col, col_vector, col_num_nz, col_indices)
end

function Highs_setBasis(highs, colstatus, rowstatus)
    ccall((:Highs_setBasis, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{HighsInt}), highs, colstatus, rowstatus)
end

function Highs_setLogicalBasis(highs)
    ccall((:Highs_setLogicalBasis, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

function Highs_getRunTime(highs)
    ccall((:Highs_getRunTime, libhighs), Cdouble, (Ptr{Cvoid},), highs)
end

function Highs_addRow(highs, lower, upper, num_new_nz, indices, values)
    ccall((:Highs_addRow, libhighs), HighsInt, (Ptr{Cvoid}, Cdouble, Cdouble, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}), highs, lower, upper, num_new_nz, indices, values)
end

function Highs_addRows(highs, num_new_row, lower, upper, num_new_nz, starts, indices, values)
    ccall((:Highs_addRows, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{Cdouble}, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, num_new_row, lower, upper, num_new_nz, starts, indices, values)
end

function Highs_addCol(highs, cost, lower, upper, num_new_nz, indices, values)
    ccall((:Highs_addCol, libhighs), HighsInt, (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}), highs, cost, lower, upper, num_new_nz, indices, values)
end

function Highs_addCols(highs, num_new_col, costs, lower, upper, num_new_nz, starts, indices, values)
    ccall((:Highs_addCols, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, num_new_col, costs, lower, upper, num_new_nz, starts, indices, values)
end

function Highs_changeObjectiveSense(highs, sense)
    ccall((:Highs_changeObjectiveSense, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt), highs, sense)
end

function Highs_changeObjectiveOffset(highs, offset)
    ccall((:Highs_changeObjectiveOffset, libhighs), HighsInt, (Ptr{Cvoid}, Cdouble), highs, offset)
end

function Highs_changeColIntegrality(highs, col, integrality)
    ccall((:Highs_changeColIntegrality, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt), highs, col, integrality)
end

function Highs_changeColsIntegralityByRange(highs, from_col, to_col, integrality)
    ccall((:Highs_changeColsIntegralityByRange, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Ptr{HighsInt}), highs, from_col, to_col, integrality)
end

function Highs_changeColsIntegralityBySet(highs, num_set_entries, set, integrality)
    ccall((:Highs_changeColsIntegralityBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}), highs, num_set_entries, set, integrality)
end

function Highs_changeColsIntegralityByMask(highs, mask, integrality)
    ccall((:Highs_changeColsIntegralityByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{HighsInt}), highs, mask, integrality)
end

function Highs_changeColCost(highs, col, cost)
    ccall((:Highs_changeColCost, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Cdouble), highs, col, cost)
end

function Highs_changeColsCostByRange(highs, from_col, to_col, cost)
    ccall((:Highs_changeColsCostByRange, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Ptr{Cdouble}), highs, from_col, to_col, cost)
end

function Highs_changeColsCostBySet(highs, num_set_entries, set, cost)
    ccall((:Highs_changeColsCostBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}), highs, num_set_entries, set, cost)
end

function Highs_changeColsCostByMask(highs, mask, cost)
    ccall((:Highs_changeColsCostByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{Cdouble}), highs, mask, cost)
end

function Highs_changeColBounds(highs, col, lower, upper)
    ccall((:Highs_changeColBounds, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Cdouble, Cdouble), highs, col, lower, upper)
end

function Highs_changeColsBoundsByRange(highs, from_col, to_col, lower, upper)
    ccall((:Highs_changeColsBoundsByRange, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Ptr{Cdouble}, Ptr{Cdouble}), highs, from_col, to_col, lower, upper)
end

function Highs_changeColsBoundsBySet(highs, num_set_entries, set, lower, upper)
    ccall((:Highs_changeColsBoundsBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}), highs, num_set_entries, set, lower, upper)
end

function Highs_changeColsBoundsByMask(highs, mask, lower, upper)
    ccall((:Highs_changeColsBoundsByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}), highs, mask, lower, upper)
end

function Highs_changeRowBounds(highs, row, lower, upper)
    ccall((:Highs_changeRowBounds, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Cdouble, Cdouble), highs, row, lower, upper)
end

function Highs_changeRowsBoundsBySet(highs, num_set_entries, set, lower, upper)
    ccall((:Highs_changeRowsBoundsBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}), highs, num_set_entries, set, lower, upper)
end

function Highs_changeRowsBoundsByMask(highs, mask, lower, upper)
    ccall((:Highs_changeRowsBoundsByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}), highs, mask, lower, upper)
end

function Highs_changeCoeff(highs, row, col, value)
    ccall((:Highs_changeCoeff, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Cdouble), highs, row, col, value)
end

function Highs_getObjectiveSense(highs, sense)
    ccall((:Highs_getObjectiveSense, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}), highs, sense)
end

function Highs_getObjectiveOffset(highs, offset)
    ccall((:Highs_getObjectiveOffset, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{Cdouble}), highs, offset)
end

function Highs_getColsByRange(highs, from_col, to_col, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getColsByRange, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, from_col, to_col, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

function Highs_getColsBySet(highs, num_set_entries, set, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getColsBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, num_set_entries, set, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

function Highs_getColsByMask(highs, mask, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getColsByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, mask, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

function Highs_getRowsByRange(highs, from_row, to_row, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getRowsByRange, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, from_row, to_row, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

function Highs_getRowsBySet(highs, num_set_entries, set, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getRowsBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, num_set_entries, set, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

function Highs_getRowsByMask(highs, mask, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getRowsByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}), highs, mask, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

function Highs_deleteColsByRange(highs, from_col, to_col)
    ccall((:Highs_deleteColsByRange, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt), highs, from_col, to_col)
end

function Highs_deleteColsBySet(highs, num_set_entries, set)
    ccall((:Highs_deleteColsBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}), highs, num_set_entries, set)
end

function Highs_deleteColsByMask(highs, mask)
    ccall((:Highs_deleteColsByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}), highs, mask)
end

function Highs_deleteRowsByRange(highs, from_row, to_row)
    ccall((:Highs_deleteRowsByRange, libhighs), HighsInt, (Ptr{Cvoid}, Cint, HighsInt), highs, from_row, to_row)
end

function Highs_deleteRowsBySet(highs, num_set_entries, set)
    ccall((:Highs_deleteRowsBySet, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Ptr{HighsInt}), highs, num_set_entries, set)
end

function Highs_deleteRowsByMask(highs, mask)
    ccall((:Highs_deleteRowsByMask, libhighs), HighsInt, (Ptr{Cvoid}, Ptr{HighsInt}), highs, mask)
end

function Highs_scaleCol(highs, col, scaleval)
    ccall((:Highs_scaleCol, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Cdouble), highs, col, scaleval)
end

function Highs_scaleRow(highs, row, scaleval)
    ccall((:Highs_scaleRow, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, Cdouble), highs, row, scaleval)
end

function Highs_getInfinity(highs)
    ccall((:Highs_getInfinity, libhighs), Cdouble, (Ptr{Cvoid},), highs)
end

function Highs_getNumCol(highs)
    ccall((:Highs_getNumCol, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

function Highs_getNumRow(highs)
    ccall((:Highs_getNumRow, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

function Highs_getNumNz(highs)
    ccall((:Highs_getNumNz, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

function Highs_getHessianNumNz(highs)
    ccall((:Highs_getHessianNumNz, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

function Highs_getModel(highs, a_format, q_format, numcol, numrow, numnz, hessian_num_nz, sense, offset, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue, qstart, qindex, qvalue, integrality)
    ccall((:Highs_getModel, libhighs), HighsInt, (Ptr{Cvoid}, HighsInt, HighsInt, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{HighsInt}), highs, a_format, q_format, numcol, numrow, numnz, hessian_num_nz, sense, offset, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue, qstart, qindex, qvalue, integrality)
end

function Highs_crossover(highs)
    ccall((:Highs_crossover, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

function Highs_crossover_set(highs, n, m, col_value, col_dual, row_dual)
    ccall((:Highs_crossover_set, libhighs), HighsInt, (Ptr{Cvoid}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), highs, n, m, col_value, col_dual, row_dual)
end

function Highs_call(numcol, numrow, numnz, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue, colvalue, coldual, rowvalue, rowdual, colbasisstatus, rowbasisstatus, modelstatus)
    ccall((:Highs_call, libhighs), HighsInt, (HighsInt, HighsInt, HighsInt, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{HighsInt}, Ptr{HighsInt}, Ptr{HighsInt}), numcol, numrow, numnz, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue, colvalue, coldual, rowvalue, rowdual, colbasisstatus, rowbasisstatus, modelstatus)
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
    ccall((:Highs_setHighsBoolOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, HighsInt), highs, option, value)
end

function Highs_setHighsIntOptionValue(highs, option, value)
    ccall((:Highs_setHighsIntOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, HighsInt), highs, option, value)
end

function Highs_setHighsDoubleOptionValue(highs, option, value)
    ccall((:Highs_setHighsDoubleOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Cdouble), highs, option, value)
end

function Highs_setHighsStringOptionValue(highs, option, value)
    ccall((:Highs_setHighsStringOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Cstring), highs, option, value)
end

function Highs_setHighsOptionValue(highs, option, value)
    ccall((:Highs_setHighsOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Cstring), highs, option, value)
end

function Highs_getHighsBoolOptionValue(highs, option, value)
    ccall((:Highs_getHighsBoolOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Ptr{HighsInt}), highs, option, value)
end

function Highs_getHighsIntOptionValue(highs, option, value)
    ccall((:Highs_getHighsIntOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Ptr{HighsInt}), highs, option, value)
end

function Highs_getHighsDoubleOptionValue(highs, option, value)
    ccall((:Highs_getHighsDoubleOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Ptr{Cdouble}), highs, option, value)
end

function Highs_getHighsStringOptionValue(highs, option, value)
    ccall((:Highs_getHighsStringOptionValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Cstring), highs, option, value)
end

function Highs_getHighsOptionType(highs, option, type)
    ccall((:Highs_getHighsOptionType, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Ptr{HighsInt}), highs, option, type)
end

function Highs_resetHighsOptions(highs)
    ccall((:Highs_resetHighsOptions, libhighs), HighsInt, (Ptr{Cvoid},), highs)
end

function Highs_getHighsIntInfoValue(highs, info, value)
    ccall((:Highs_getHighsIntInfoValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Ptr{HighsInt}), highs, info, value)
end

function Highs_getHighsDoubleInfoValue(highs, info, value)
    ccall((:Highs_getHighsDoubleInfoValue, libhighs), HighsInt, (Ptr{Cvoid}, Cstring, Ptr{Cdouble}), highs, info, value)
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
# Julia wrapper for header: HighsInt.h
# Automatically generated using Clang.jl
