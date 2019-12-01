using Test

import HiGHS

function test_model()
    cc = [1.0, -2.0]
    cl = [0.0, 0.0]
    cu = [10.0, 10.0]
    ru = [2.0, 1.0]
    rl = [0.0, 0.0]
    astart = [0, 2, 4]
    aindex = [0, 1, 0, 1]
    avalue = [1.0, 2.0, 1.0, 3.0]
    return (cc, cl, cu, rl, ru, astart, aindex, avalue)
end

@testset "Direct C call" begin
    (colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue) = test_model()
    n_col = convert(Cint, size(colcost, 1))
    n_row = convert(Cint, size(rowlower, 1))
    n_nz = convert(Cint, size(aindex, 1))

    colcost = convert(Array{Cdouble}, colcost)
    collower = convert(Array{Cdouble}, collower)
    colupper = convert(Array{Cdouble}, colupper)

    rowlower = convert(Array{Cdouble}, rowlower)
    rowupper = convert(Array{Cdouble}, rowupper)
    matstart = convert(Array{Cint}, astart)
    matindex = convert(Array{Cint}, aindex)
    matvalue = convert(Array{Cdouble}, avalue)

    # solution = HighsSolution(Array{Cdouble, 1}(undef, n_col), Array{Cdouble, 1}(undef, n_col), Array{Cdouble, 1}(undef, n_row),  Array{Cdouble, 1}(undef, n_row))
    # basis = HighsBasis(Array{Cint, 1}(undef, n_col), Array{Cint, 1}(undef, n_row))
    colvalue, coldual = (Array{Cdouble, 1}(undef, n_col), Array{Cdouble, 1}(undef, n_col))
    rowvalue, rowdual = (Array{Cdouble, 1}(undef, n_row),  Array{Cdouble, 1}(undef, n_row))
    colbasisstatus, rowbasisstatus = (Array{Cint, 1}(undef, n_col), Array{Cint, 1}(undef, n_row))

    modelstatus = Ref{Cint}(42)
    status = HiGHS.CWrapper.Highs_call(n_col, n_row, n_nz, colcost, collower, colupper, rowlower, rowupper, matstart, matindex, matvalue, colvalue, coldual, rowvalue, rowdual, colbasisstatus, rowbasisstatus, modelstatus)
    @test status == 0
    @test modelstatus[] == 11 # optimal
end

@testset "Managed HiGHS" begin
    @test_nowarn finalize(HiGHS.ManagedHiGHS())
    managed_h = HiGHS.ManagedHiGHS()
    @test HiGHS.free_highs(managed_h)
    @test !HiGHS.free_highs(managed_h)
end
