import MathOptInterface
import HiGHS

const MOI = MathOptInterface
const MOIT = MOI.Test

const CONFIG = MOIT.TestConfig()

@testset "Attributes" begin
    o = HiGHS.Optimizer()
    @test MOI.supports(o, MOI.SolverName())
    @test MOI.get(o, MOI.SolverName()) == "HiGHS"
    @test MOI.get(o, MOI.TimeLimitSec()) > 10000
    MOI.set(o, MOI.TimeLimitSec(), 500)
    @test MOI.get(o, MOI.TimeLimitSec()) == 500.0
    @test MOI.supports(o, MOI.RawSolver())
    @test MOI.get(o, MOI.RawSolver()) == o.model
end

@testset "MOI variable count and empty" begin
    o = HiGHS.Optimizer()
    x1 = MOI.add_variable(o)
    @test x1.value == 0
    @test MOI.supports_constraint(o, MOI.SingleVariable(x1), MOI.Interval(0, 1))
    (x2, _) = MOI.add_constrained_variable(o, MOI.Interval(0, 1))
    @test x2.value == 1
    @test MOI.get(o, MOI.NumberOfVariables()) == 2
    MOI.empty!(o)
    @test MOI.get(o, MOI.NumberOfVariables()) == 0
end

@testset "Objective function and value in box constraints" begin
    o = HiGHS.Optimizer()
    @testset "Getting objective value" begin
        (x, _) = MOI.add_constrained_variable(o, MOI.Interval(-3.0, 6.0))
        HiGHS.CWrapper.Highs_changeColCost(o.model.inner, Cint(x.value), 1.0)
        @test MOI.get(o, MOI.ObjectiveSense()) == MOI.FEASIBILITY_SENSE
        MOI.set(o, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        @test MOI.get(o, MOI.ObjectiveSense()) == MOI.MIN_SENSE
        @test MOI.get(o, MOI.ResultCount()) == 0
        MOI.optimize!(o)
        @test MOI.get(o, MOI.ResultCount()) == 1
        @test MOI.get(o, MOI.ObjectiveValue()) ≈ -3
    end
    @testset "Max in box" begin
        MOI.empty!(o)
        @test MOI.get(o, MOI.ResultCount()) == 0
        (x, _) = MOI.add_constrained_variable(o, MOI.Interval(-3.0, 6.0))
        MOI.set(o, MOI.ObjectiveSense(), MOI.MAX_SENSE)
        HiGHS.CWrapper.Highs_changeColCost(o.model.inner, Cint(x.value), 2.0)
        MOI.optimize!(o)
        # BUG in HiGHS
        # see https://github.com/ERGO-Code/HiGHS/issues/316
        @test_broken MOI.get(o, MOI.ObjectiveValue()) ≈ 2 * 6
        obj_func = MOI.get(o, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
        @test obj_func ≈ MOI.ScalarAffineFunction([
                MOI.ScalarAffineTerm(2.0, x),
            ], 0.0,
        )
    end
    @testset "Objective function obtained from model corresponds" begin
        MOI.empty!(o)
        (x1, _) = MOI.add_constrained_variable(o, MOI.Interval(-3.0, 6.0))
        (x2, _) = MOI.add_constrained_variable(o, MOI.Interval(1.0, 2.0))
        MOI.set(o, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        HiGHS.CWrapper.Highs_changeColCost(o.model.inner, Cint(x1.value), 2.0)
        HiGHS.CWrapper.Highs_changeColCost(o.model.inner, Cint(x2.value), -1.0)
        F = MOI.get(o, MOI.ObjectiveFunctionType())
        @test F <: MOI.ScalarAffineFunction{Float64}
        obj_func = MOI.get(o, MOI.ObjectiveFunction{F}())
        @test MOI.supports(o, MOI.ObjectiveFunction{F}())
        @test all(MOI.get(o, MOI.ListOfVariableIndices()) .== [x1, x2])        
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
        obj_func.terms[1] = MOI.ScalarAffineTerm(3.0, x1)
        MOI.set(o, MOI.ObjectiveFunction{F}(), obj_func)
        obj_func = MOI.get(o, MOI.ObjectiveFunction{F}())
        @test obj_func ≈ MOI.ScalarAffineFunction([
                MOI.ScalarAffineTerm(3.0, x1),
                MOI.ScalarAffineTerm(-1.0, x2),
            ], 0.0,
        )
    end
    
    @testset "Constrained variable equivalent to add constraint" begin
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

    @testset "Constant in objective function" begin
        MOI.empty!(o)
        x = MOI.add_variable(o)
        _ = MOI.add_constraint(o, MOI.SingleVariable(x), MOI.Interval(-3.0, 6.0))
        MOI.set(o, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        obj_func = MOI.ScalarAffineFunction(
            [MOI.ScalarAffineTerm(1.0, x)], 3.0,
        )
        MOI.set(o, MOI.ObjectiveFunction{typeof(obj_func)}(), obj_func)
        MOI.optimize!(o)
        @test MOI.get(o, MOI.ResultCount()) == 1
        @test MOI.get(o, MOI.ObjectiveValue()) ≈ 0
        obj_func = MOI.get(o, MOI.ObjectiveFunction{typeof(obj_func)}())
        @test MOI.constant(obj_func) ≈ 3        
        MOI.set(o, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
        obj_func = MOI.get(o, MOI.ObjectiveFunction{typeof(obj_func)}())
        @test MOI.constant(obj_func) ≈ 0
        @test isempty(obj_func.terms)
    end
end

@testset "Linear constraints" begin
    # max x1 + 2x2
    # st 0 <= x{1,2} <= 5
    # 0 <= x1 + x2 <= 7.5
    o = HiGHS.Optimizer()
    (x1, _) = MOI.add_constrained_variable(o, MOI.Interval(0.0, 5.0))
    (x2, _) = MOI.add_constrained_variable(o, MOI.Interval(0.0, 5.0))
    MOI.set(o, MOI.ObjectiveSense(), MOI.MAX_SENSE)    
    func = MOI.ScalarAffineFunction(
        [
            MOI.ScalarAffineTerm(1.0, x1),
            MOI.ScalarAffineTerm(2.0, x2),
        ], 0.0,
    )
    @test MOI.supports_constraint(o, func, MOI.Interval(0, 1))
    MOI.set(o, MOI.ObjectiveFunction{typeof(func)}(), func)    
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
    @test MOI.get(o, MOI.BarrierIterations()) == 0
end

@testset "Variable names" begin
    o = HiGHS.Optimizer()
    MOIT.variablenames(o, CONFIG)
    MOI.empty!(o)
    y = MOI.add_variable(o)
    MOI.set(o, MOI.VariableName(), y, "y")
    y2 = MOI.get(o, MOI.VariableIndex, "y")
    @test y == y2
    @test MOI.get(o, MOI.VariableIndex, "y0") === nothing
end

@testset "HiGHS custom options" begin
    o = HiGHS.Optimizer()
    @test MOI.supports(o, MOI.RawParameter("solver"))
    @test MOI.get(o, MOI.RawParameter("solver")) == "choose"
    MOI.set(o, MOI.RawParameter("solver"), "simplex")
    @test MOI.get(o, MOI.RawParameter("solver")) == "simplex"
    @test MOI.get(o, MOI.RawParameter("message_level")) == 4
    MOI.set(o, MOI.RawParameter("message_level"), 1)
    @test MOI.get(o, MOI.RawParameter("message_level")) == 1
    @test MOI.get(o, MOI.RawParameter("time_limit")) > 1000
    MOI.set(o, MOI.RawParameter("time_limit"), 1000.0)
    @test MOI.get(o, MOI.RawParameter("time_limit")) == 1000.0
    # unsupported test
    @test_throws ArgumentError MOI.get(o, MOI.RawParameter("wrong_param"))
end

@testset "Model empty" begin
    o = HiGHS.Optimizer()
    @test MOI.is_empty(o)
    MOI.add_variable(o)
    @test !MOI.is_empty(o)
    MOI.empty!(o)
    MOI.set(
        o,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction{Float64}([], 0.0)
    )
    @test MOI.is_empty(o)
    MOI.set(
        o,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction{Float64}([], 3.0)
    )
    @test !MOI.is_empty(o)
    MOI.set(o, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
    @test MOI.is_empty(o)
    @test MOI.get(o, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}()) ≈ MOI.ScalarAffineFunction{Float64}([], 0.0)
    x = MOI.add_variable(o)
    MOI.set(
        o,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm(1.0, x)], 0.0)
    )
    @test_throws ErrorException MOI.optimize!(o)
end