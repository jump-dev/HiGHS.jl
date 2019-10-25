# Julia wrapper for header: highs_c_api.h
# Automatically generated using Clang.jl


function Highs_call(numcol::Cint, numrow::Cint, numnz::Cint, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue, colvalue, coldual, rowvalue, rowdual, colbasisstatus, rowbasisstatus, modelstatus)
    ccall((:Highs_call, libhighs), Cint, (Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), numcol, numrow, numnz, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue, colvalue, coldual, rowvalue, rowdual, colbasisstatus, rowbasisstatus, modelstatus)
end

function Highs_create()
    ccall((:Highs_create, libhighs), Ptr{Cvoid}, ())
end

function Highs_destroy(highs)
    ccall((:Highs_destroy, libhighs), Cvoid, (Ptr{Cvoid},), highs)
end

function Highs_run(highs)
    ccall((:Highs_run, libhighs), Cint, (Ptr{Cvoid},), highs)
end

function Highs_readFromFile(highs, filename)
    ccall((:Highs_readFromFile, libhighs), Cint, (Ptr{Cvoid}, Cstring), highs, filename)
end

function Highs_writeToFile(highs, filename)
    ccall((:Highs_writeToFile, libhighs), Cint, (Ptr{Cvoid}, Cstring), highs, filename)
end

function Highs_loadModel(highs, numcol::Cint, numrow::Cint, numnz::Cint, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue)
    ccall((:Highs_loadModel, libhighs), Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), highs, numcol, numrow, numnz, colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue)
end

function Highs_setOptionValue(highs, option, value)
    ccall((:Highs_setOptionValue, libhighs), Cint, (Ptr{Cvoid}, Cstring, Cstring), highs, option, value)
end

function Highs_getSolution(highs, colvalue, coldual, rowvalue, rowdual)
    ccall((:Highs_getSolution, libhighs), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), highs, colvalue, coldual, rowvalue, rowdual)
end

function Highs_getBasis(highs, colstatus, rowstatus)
    ccall((:Highs_getBasis, libhighs), Cvoid, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}), highs, colstatus, rowstatus)
end

function Highs_getObjectiveValue(highs)
    ccall((:Highs_getObjectiveValue, libhighs), Cdouble, (Ptr{Cvoid},), highs)
end

function Highs_getIterationCount(highs)
    ccall((:Highs_getIterationCount, libhighs), Cint, (Ptr{Cvoid},), highs)
end

function Highs_addRow(highs, lower::Cdouble, upper::Cdouble, num_new_nz::Cint, indices, values)
    ccall((:Highs_addRow, libhighs), Cint, (Ptr{Cvoid}, Cdouble, Cdouble, Cint, Ptr{Cint}, Ptr{Cdouble}), highs, lower, upper, num_new_nz, indices, values)
end

function Highs_addRows(highs, num_new_row::Cint, lower, upper, num_new_nz::Cint, starts, indices, values)
    ccall((:Highs_addRows, libhighs), Cint, (Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), highs, num_new_row, lower, upper, num_new_nz, starts, indices, values)
end

function Highs_addCol(highs, cost::Cdouble, lower::Cdouble, upper::Cdouble, num_new_nz::Cint, indices, values)
    ccall((:Highs_addCol, libhighs), Cint, (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Cint, Ptr{Cint}, Ptr{Cdouble}), highs, cost, lower, upper, num_new_nz, indices, values)
end

function Highs_addCols(highs, num_new_col::Cint, costs, lower, upper, num_new_nz::Cint, starts, indices, values)
    ccall((:Highs_addCols, libhighs), Cint, (Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), highs, num_new_col, costs, lower, upper, num_new_nz, starts, indices, values)
end

function Highs_changeObjectiveSense(highs, sense::Cint)
    ccall((:Highs_changeObjectiveSense, libhighs), Cint, (Ptr{Cvoid}, Cint), highs, sense)
end

function Highs_changeColCost(highs, col::Cint, cost::Cdouble)
    ccall((:Highs_changeColCost, libhighs), Cint, (Ptr{Cvoid}, Cint, Cdouble), highs, col, cost)
end

function Highs_changeColsCostBySet(highs, num_set_entries::Cint, set, cost)
    ccall((:Highs_changeColsCostBySet, libhighs), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}), highs, num_set_entries, set, cost)
end

function Highs_changeColsCostByMask(highs, mask, cost)
    ccall((:Highs_changeColsCostByMask, libhighs), Cint, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cdouble}), highs, mask, cost)
end

function Highs_changeColBounds(highs, col::Cint, lower::Cdouble, upper::Cdouble)
    ccall((:Highs_changeColBounds, libhighs), Cint, (Ptr{Cvoid}, Cint, Cdouble, Cdouble), highs, col, lower, upper)
end

function Highs_changeColsBoundsByRange(highs, from_col::Cint, to_col::Cint, lower, upper)
    ccall((:Highs_changeColsBoundsByRange, libhighs), Cint, (Ptr{Cvoid}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}), highs, from_col, to_col, lower, upper)
end

function Highs_changeColsBoundsBySet(highs, num_set_entries::Cint, set, lower, upper)
    ccall((:Highs_changeColsBoundsBySet, libhighs), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}), highs, num_set_entries, set, lower, upper)
end

function Highs_changeColsBoundsByMask(highs, mask, lower, upper)
    ccall((:Highs_changeColsBoundsByMask, libhighs), Cint, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}), highs, mask, lower, upper)
end

function Highs_changeRowBounds(highs, row::Cint, lower::Cdouble, upper::Cdouble)
    ccall((:Highs_changeRowBounds, libhighs), Cint, (Ptr{Cvoid}, Cint, Cdouble, Cdouble), highs, row, lower, upper)
end

function Highs_changeRowsBoundsBySet(highs, num_set_entries::Cint, set, lower, upper)
    ccall((:Highs_changeRowsBoundsBySet, libhighs), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}), highs, num_set_entries, set, lower, upper)
end

function Highs_changeRowsBoundsByMask(highs, mask, lower, upper)
    ccall((:Highs_changeRowsBoundsByMask, libhighs), Cint, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}), highs, mask, lower, upper)
end

function Highs_getColsByRange(highs, from_col::Cint, to_col::Cint, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getColsByRange, libhighs), Cint, (Ptr{Cvoid}, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), highs, from_col, to_col, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

function Highs_getColsBySet(highs, num_set_entries::Cint, set, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getColsBySet, libhighs), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), highs, num_set_entries, set, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

function Highs_getColsByMask(highs, mask, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getColsByMask, libhighs), Cint, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), highs, mask, num_col, costs, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

function Highs_getRowsByRange(highs, from_row::Cint, to_row::Cint, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getRowsByRange, libhighs), Cint, (Ptr{Cvoid}, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), highs, from_row, to_row, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

function Highs_getRowsBySet(highs, num_set_entries::Cint, set, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getRowsBySet, libhighs), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), highs, num_set_entries, set, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

function Highs_getRowsByMask(highs, mask, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
    ccall((:Highs_getRowsByMask, libhighs), Cint, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), highs, mask, num_row, lower, upper, num_nz, matrix_start, matrix_index, matrix_value)
end

function Highs_deleteColsByRange(highs, from_col::Cint, to_col::Cint)
    ccall((:Highs_deleteColsByRange, libhighs), Cint, (Ptr{Cvoid}, Cint, Cint), highs, from_col, to_col)
end

function Highs_deleteColsBySet(highs, num_set_entries::Cint, set)
    ccall((:Highs_deleteColsBySet, libhighs), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}), highs, num_set_entries, set)
end

function Highs_deleteColsByMask(highs, mask)
    ccall((:Highs_deleteColsByMask, libhighs), Cint, (Ptr{Cvoid}, Ptr{Cint}), highs, mask)
end

function Highs_deleteRowsByRange(highs, from_row::Cint, to_row::Cint)
    ccall((:Highs_deleteRowsByRange, libhighs), Cint, (Ptr{Cvoid}, Cint, Cint), highs, from_row, to_row)
end

function Highs_deleteRowsBySet(highs, num_set_entries::Cint, set)
    ccall((:Highs_deleteRowsBySet, libhighs), Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}), highs, num_set_entries, set)
end

function Highs_deleteRowsByMask(highs, mask)
    ccall((:Highs_deleteRowsByMask, libhighs), Cint, (Ptr{Cvoid}, Ptr{Cint}), highs, mask)
end

function Highs_getNumCols(highs)
    ccall((:Highs_getNumCols, libhighs), Cint, (Ptr{Cvoid},), highs)
end

function Highs_getNumRows(highs)
    ccall((:Highs_getNumRows, libhighs), Cint, (Ptr{Cvoid},), highs)
end

function Highs_getNumNz(highs)
    ccall((:Highs_getNumNz, libhighs), Cint, (Ptr{Cvoid},), highs)
end

function Highs_getModelStatus(highs)
    ccall((:Highs_getModelStatus, libhighs), Cint, (Ptr{Cvoid},), highs)
end

function Highs_getBasicVariables(highs, basic_variables)
    ccall((:Highs_getBasicVariables, libhighs), Cint, (Ptr{Cvoid}, Ptr{Cint}), highs, basic_variables)
end

function Highs_getBasisInverseRow(highs, row::Cint, row_vector, row_num_nz, row_indices)
    ccall((:Highs_getBasisInverseRow, libhighs), Cint, (Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}), highs, row, row_vector, row_num_nz, row_indices)
end

function Highs_getBasisInverseCol(highs, col::Cint, col_vector, col_num_nz, col_indices)
    ccall((:Highs_getBasisInverseCol, libhighs), Cint, (Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}), highs, col, col_vector, col_num_nz, col_indices)
end

function Highs_getBasisSolve(highs, rhs, solution_vector, solution_num_nz, solution_indices)
    ccall((:Highs_getBasisSolve, libhighs), Cint, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}), highs, rhs, solution_vector, solution_num_nz, solution_indices)
end

function Highs_getBasisTransposeSolve(highs, rhs, solution_vector, solution_nz, solution_indices)
    ccall((:Highs_getBasisTransposeSolve, libhighs), Cint, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}), highs, rhs, solution_vector, solution_nz, solution_indices)
end

function Highs_getReducedRow(highs, row::Cint, row_vector, row_num_nz, row_indices)
    ccall((:Highs_getReducedRow, libhighs), Cint, (Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}), highs, row, row_vector, row_num_nz, row_indices)
end

function Highs_getReducedColumn(highs, col::Cint, col_vector, col_num_nz, col_indices)
    ccall((:Highs_getReducedColumn, libhighs), Cint, (Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}), highs, col, col_vector, col_num_nz, col_indices)
end
