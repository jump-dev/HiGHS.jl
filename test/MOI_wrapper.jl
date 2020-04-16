import MathOptInterface
import HiGHS

const MOI = MathOptInterface
const MOIT = MOI.Test

const CONFIG = MOIT.TestConfig()

@testset "SolverName" begin
    @test MOI.get(HiGHS.Optimizer(), MOI.SolverName()) == "HiGHS"
end

@testset "MOI variable count and empty" begin
    o = HiGHS.Optimizer()
    x1 = MOI.add_variable(o)
    @test x1.value == 0
    (x2, _) = MOI.add_constrained_variable(o, MOI.Interval(0, 1))
    @test x2.value == 1
    @test MOI.get(o, MOI.NumberOfVariables()) == 2
    MOI.empty!(o)
    @test MOI.get(o, MOI.NumberOfVariables()) == 0
end

@testset "Objective function and value in box constraints" begin
    o = HiGHS.Optimizer()
    (x, _) = MOI.add_constrained_variable(o, MOI.Interval(-3.0, 6.0))
    HiGHS.CWrapper.Highs_changeColCost(o.model.inner, Cint(x.value), 1.0)
    MOI.set(o, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    @test MOI.get(o, MOI.ResultCount()) == 0
    MOI.optimize!(o)
    @test MOI.get(o, MOI.ResultCount()) == 1
    @test MOI.get(o, MOI.ObjectiveValue()) ≈ -3
    MOI.empty!(o)
    @test MOI.get(o, MOI.ResultCount()) == 0
    (x, _) = MOI.add_constrained_variable(o, MOI.Interval(-3.0, 6.0))
    MOI.set(o, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    HiGHS.CWrapper.Highs_changeColCost(o.model.inner, Cint(x.value), 2.0)
    MOI.optimize!(o)
    # BUG in HiGHS
    # see https://github.com/ERGO-Code/HiGHS/issues/316
    # @test MOI.get(o, MOI.ObjectiveValue()) ≈ 2 * 6
    obj_func = MOI.get(o, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    @test obj_func ≈ MOI.ScalarAffineFunction([
            MOI.ScalarAffineTerm(2.0, x),
        ], 0.0,
    )

    MOI.empty!(o)
    (x1, _) = MOI.add_constrained_variable(o, MOI.Interval(-3.0, 6.0))
    (x2, _) = MOI.add_constrained_variable(o, MOI.Interval(1.0, 2.0))
    MOI.set(o, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    HiGHS.CWrapper.Highs_changeColCost(o.model.inner, Cint(x1.value), 2.0)
    HiGHS.CWrapper.Highs_changeColCost(o.model.inner, Cint(x2.value), -1.0)
    F = MOI.get(o, MOI.ObjectiveFunctionType())
    @test F <: MOI.ScalarAffineFunction{Float64}
    obj_func = MOI.get(o, MOI.ObjectiveFunction{F}())
    @test obj_func ≈ MOI.ScalarAffineFunction([
            MOI.ScalarAffineTerm(2.0, x1),
            MOI.ScalarAffineTerm(-1.0, x2),
        ], 0.0,
    )
    MOI.set(o, MOI.ObjectiveFunction{F}(), obj_func)
    obj_func = MOI.get(o, MOI.ObjectiveFunction{F}())
    @test obj_func ≈ MOI.ScalarAffineFunction([
            MOI.ScalarAffineTerm(2.0, x1),
            MOI.ScalarAffineTerm(-1.0, x2),
        ], 0.0,
    )
    @test all(MOI.get(o, MOI.ListOfVariableIndices()) .== [x1, x2])
    # add constraint variable equivalent to add constraint
    MOI.empty!(o)
    x = MOI.add_variable(o)
    _ = MOI.add_constraint(o, MOI.SingleVariable(x), MOI.Interval(-3.0, 6.0))
    HiGHS.CWrapper.Highs_changeColCost(o.model.inner, Cint(x.value), 1.0)
    MOI.set(o, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    @test MOI.get(o, MOI.ResultCount()) == 0
    MOI.optimize!(o)
    @test MOI.get(o, MOI.ResultCount()) == 1
    @test MOI.get(o, MOI.ObjectiveValue()) ≈ -3
end

@testset "Linear constraints" begin
    # max x1 + 2x2
    # st 0 <= x{1,2} <= 5
    # 0 <= x1 + x2 <= 7.5
    o = HiGHS.Optimizer()
    (x1, _) = MOI.add_constrained_variable(o, MOI.Interval(0.0, 5.0))
    (x2, _) = MOI.add_constrained_variable(o, MOI.Interval(0.0, 5.0))
    MOI.set(o, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(o, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(
            [
                MOI.ScalarAffineTerm(1.0, x1),
                MOI.ScalarAffineTerm(2.0, x2),
            ], 0.0,
        )
    )
    @test MOI.get(o, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64}}()) == 0
    MOI.add_constraint(o,
        MOI.ScalarAffineFunction(
            [
                MOI.ScalarAffineTerm(1.0, x1),
                MOI.ScalarAffineTerm(1.0, x2),
            ], 0.0,
        ), MOI.Interval(0.0, 7.5),
    )
    @test MOI.get(o, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64}}()) == 1
    MOI.optimize!(o)
    @test MOI.get(o, MOI.ObjectiveValue()) ≈ 12.5
    @test MOI.get(o, MOI.SimplexIterations()) > 0
end

@testset "Variable names" begin
    MOIT.variablenames(HiGHS.Optimizer(), CONFIG)
end
