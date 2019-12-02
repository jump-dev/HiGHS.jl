import MathOptInterface
const MOI = MathOptInterface

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

@testset "Objective function and value" begin
    o = HiGHS.Optimizer()
    (x, _) = MOI.add_constrained_variable(o, MOI.Interval(-3.0, 6.0))
    HiGHS.CWrapper.Highs_changeColCost(o.model.inner, Cint(x.value), 1.0)
    MOI.set(o, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.optimize!(o)
    @test MOI.get(o, MOI.ObjectiveValue()) ≈ -3
    MOI.empty!(o)
    (x, _) = MOI.add_constrained_variable(o, MOI.Interval(-3.0, 6.0))
    MOI.set(o, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    HiGHS.CWrapper.Highs_changeColCost(o.model.inner, Cint(x.value), 2.0)
    MOI.optimize!(o)
    @test MOI.get(o, MOI.ObjectiveValue()) ≈ 2 * 6
    obj_func = MOI.get(o, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    @test obj_func ≈ MOI.ScalarAffineFunction([
            MOI.ScalarAffineTerm(2.0, x),
        ], 0.0,
    )
end
