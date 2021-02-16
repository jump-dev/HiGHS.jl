module TestMOIHighs

import HiGHS
import MathOptInterface
using Test

const MOI = MathOptInterface

const OPTIMIZER = MOI.Bridges.full_bridge_optimizer(HiGHS.Optimizer(), Float64)
MOI.set(OPTIMIZER, MOI.Silent(), true)

const CONFIG = MOI.Test.TestConfig(
    basis = true,

    # TODO(odow): add infeasibility certificates.
    infeas_certificates = false,
)

function test_basic_constraint_tests()
    return MOI.Test.basic_constraint_tests(OPTIMIZER, CONFIG)
end

function test_unittest()
    return MOI.Test.unittest(
        OPTIMIZER,
        CONFIG,
        String[
            # TODO(odow):
            # Add support for MOI.ObjectiveBound.
            "solve_objbound_edge_cases",
            # Add support for MOI.NumberOfThreads.
            "number_threads",

            # These are excluded because HiGHS does not support them.
            "delete_soc_variables",
            "solve_integer_edge_cases",
            "solve_qcp_edge_cases",
            "solve_qp_edge_cases",
            "solve_zero_one_with_bounds_1",
            "solve_zero_one_with_bounds_2",
            "solve_zero_one_with_bounds_3",
        ],
    )
end

function test_modificationtest()
    return MOI.Test.modificationtest(OPTIMIZER, CONFIG)
end

function test_contlineartest()
    return MOI.Test.contlineartest(
        OPTIMIZER,
        CONFIG,
        String[
            # Upstream segfault. Reported: https://github.com/ERGO-Code/HiGHS/issues/448
            "linear8b",

            # VariablePrimalStart not supported.
            "partial_start",
        ],
    )
end

function test_lintest()
    return MOI.Test.lintest(OPTIMIZER, CONFIG)
end

function test_SolverName()
    @test MOI.get(OPTIMIZER, MOI.SolverName()) == "HiGHS"
end

function test_default_objective()
    return MOI.Test.default_objective_test(OPTIMIZER)
end

function test_default_status_test()
    return MOI.Test.default_status_test(OPTIMIZER)
end

function test_nametest()
    return MOI.Test.nametest(OPTIMIZER)
end

function test_validtest()
    return MOI.Test.validtest(OPTIMIZER)
end

function test_emptytest()
    return MOI.Test.emptytest(OPTIMIZER)
end

function test_orderedindicestest()
    return MOI.Test.orderedindicestest(OPTIMIZER)
end

function test_copytest()
    return MOI.Test.copytest(
        OPTIMIZER,
        MOI.Bridges.full_bridge_optimizer(HiGHS.Optimizer(), Float64),
    )
end

function test_scalar_function_constant_not_zero()
    return MOI.Test.scalar_function_constant_not_zero(OPTIMIZER)
end

function test_supports_constrainttest()
    # supports_constrainttest needs VectorOfVariables-in-Zeros,
    # MOI.Test.supports_constrainttest(HiGHS.Optimizer(), Float64, Float32)
    # but supports_constrainttest is broken via bridges:
    MOI.empty!(OPTIMIZER)
    MOI.add_variable(OPTIMIZER)
    @test MOI.supports_constraint(
        OPTIMIZER,
        MOI.SingleVariable,
        MOI.EqualTo{Float64},
    )
    @test MOI.supports_constraint(
        OPTIMIZER,
        MOI.ScalarAffineFunction{Float64},
        MOI.EqualTo{Float64},
    )
    # This test is broken for some reason:
    @test_broken !MOI.supports_constraint(
        OPTIMIZER,
        MOI.ScalarAffineFunction{Int},
        MOI.EqualTo{Float64},
    )
    @test !MOI.supports_constraint(
        OPTIMIZER,
        MOI.ScalarAffineFunction{Int},
        MOI.EqualTo{Int},
    )
    @test !MOI.supports_constraint(
        OPTIMIZER,
        MOI.SingleVariable,
        MOI.EqualTo{Int},
    )
    @test MOI.supports_constraint(OPTIMIZER, MOI.VectorOfVariables, MOI.Zeros)
    @test !MOI.supports_constraint(
        OPTIMIZER,
        MOI.VectorOfVariables,
        MOI.EqualTo{Float64},
    )
    @test !MOI.supports_constraint(OPTIMIZER, MOI.SingleVariable, MOI.Zeros)
    @test !MOI.supports_constraint(
        OPTIMIZER,
        MOI.VectorOfVariables,
        MOI.Test.UnknownVectorSet,
    )
end

function test_set_lower_bound_twice()
    return MOI.Test.set_lower_bound_twice(HiGHS.Optimizer(), Float64)
end

function test_set_upper_bound_twice()
    return MOI.Test.set_upper_bound_twice(HiGHS.Optimizer(), Float64)
end

function test_Attributes()
    o = HiGHS.Optimizer()
    @test MOI.get(o, MOI.SolverName()) == "HiGHS"
    @test MOI.get(o, MOI.TimeLimitSec()) > 10000
    MOI.set(o, MOI.TimeLimitSec(), 500)
    @test MOI.get(o, MOI.TimeLimitSec()) == 500.0
    @test MOI.get(o, MOI.RawSolver()) == o
end

function test_MOI_variable_count_and_empty()
    o = HiGHS.Optimizer()
    @test MOI.get(o, MOI.NumberOfVariables()) == 0
    x1 = MOI.add_variable(o)
    @test MOI.get(o, MOI.NumberOfVariables()) == 1
    @test MOI.supports_constraint(o, MOI.SingleVariable, MOI.Interval{Float64})
    x2, _ = MOI.add_constrained_variable(o, MOI.Interval(0.0, 1.0))
    @test MOI.get(o, MOI.NumberOfVariables()) == 2
    MOI.empty!(o)
    @test MOI.get(o, MOI.NumberOfVariables()) == 0
end

function test_Getting_objective_value()
    o = HiGHS.Optimizer()
    MOI.set(o, MOI.Silent(), true)
    (x, _) = MOI.add_constrained_variable(o, MOI.Interval(-3.0, 6.0))
    MOI.set(
        o,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x)], 0.0),
    )
    @test MOI.get(o, MOI.ObjectiveSense()) == MOI.FEASIBILITY_SENSE
    MOI.set(o, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    @test MOI.get(o, MOI.ObjectiveSense()) == MOI.MIN_SENSE
    @test MOI.get(o, MOI.ResultCount()) == 0
    MOI.optimize!(o)
    @test MOI.get(o, MOI.ResultCount()) == 1
    @test MOI.get(o, MOI.ObjectiveValue()) ≈ -3
end

function test_Max_in_box()
    o = HiGHS.Optimizer()
    MOI.set(o, MOI.Silent(), true)
    @test MOI.get(o, MOI.ResultCount()) == 0
    (x, _) = MOI.add_constrained_variable(o, MOI.Interval(-3.0, 6.0))
    MOI.set(o, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(
        o,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(2.0, x)], 0.0),
    )
    MOI.optimize!(o)
    @test MOI.get(o, MOI.ObjectiveValue()) ≈ 2 * 6
    obj_func =
        MOI.get(o, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    @test MOI.get(o, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test obj_func ≈
          MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(2.0, x)], 0.0)
end

function test_Objective_function_obtained_from_model_corresponds()
    o = HiGHS.Optimizer()
    (x1, _) = MOI.add_constrained_variable(o, MOI.Interval(-3.0, 6.0))
    (x2, _) = MOI.add_constrained_variable(o, MOI.Interval(1.0, 2.0))
    MOI.set(o, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(
        o,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(
            MOI.ScalarAffineTerm.([2.0, -1.0], [x1, x2]),
            0.0,
        ),
    )
    F = MOI.get(o, MOI.ObjectiveFunctionType())
    @test F <: MOI.ScalarAffineFunction{Float64}
    obj_func = MOI.get(o, MOI.ObjectiveFunction{F}())
    @test MOI.supports(o, MOI.ObjectiveFunction{F}())
    @test all(MOI.get(o, MOI.ListOfVariableIndices()) .== [x1, x2])
    @test obj_func ≈ MOI.ScalarAffineFunction(
        [MOI.ScalarAffineTerm(2.0, x1), MOI.ScalarAffineTerm(-1.0, x2)],
        0.0,
    )
    MOI.set(o, MOI.ObjectiveFunction{F}(), obj_func)
    obj_func = MOI.get(o, MOI.ObjectiveFunction{F}())
    @test obj_func ≈ MOI.ScalarAffineFunction(
        [MOI.ScalarAffineTerm(2.0, x1), MOI.ScalarAffineTerm(-1.0, x2)],
        0.0,
    )
    obj_func.terms[1] = MOI.ScalarAffineTerm(3.0, x1)
    MOI.set(o, MOI.ObjectiveFunction{F}(), obj_func)
    obj_func = MOI.get(o, MOI.ObjectiveFunction{F}())
    @test obj_func ≈ MOI.ScalarAffineFunction(
        [MOI.ScalarAffineTerm(3.0, x1), MOI.ScalarAffineTerm(-1.0, x2)],
        0.0,
    )
end

function test_Constrained_variable_equivalent_to_add_constraint()
    o = HiGHS.Optimizer()
    MOI.set(o, MOI.Silent(), true)
    x = MOI.add_variable(o)
    _ = MOI.add_constraint(o, MOI.SingleVariable(x), MOI.Interval(-3.0, 6.0))
    MOI.set(
        o,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x)], 0.0),
    )
    MOI.set(o, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    @test MOI.get(o, MOI.ResultCount()) == 0
    MOI.optimize!(o)
    @test MOI.get(o, MOI.ResultCount()) == 1
    @test MOI.get(o, MOI.ObjectiveValue()) ≈ -3
end

function test_Constant_in_objective_function()
    o = HiGHS.Optimizer()
    MOI.set(o, MOI.Silent(), true)
    x = MOI.add_variable(o)
    _ = MOI.add_constraint(o, MOI.SingleVariable(x), MOI.Interval(-3.0, 6.0))
    MOI.set(o, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj_func = MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x)], 3.0)
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

function test_Linear_constraints()
    # max x1 + 2x2
    # st 0 <= x{1,2} <= 5
    # 0 <= x1 + x2 <= 7.5
    o = HiGHS.Optimizer()
    MOI.set(o, MOI.Silent(), true)
    (x1, _) = MOI.add_constrained_variable(o, MOI.Interval(0.0, 5.0))
    (x2, _) = MOI.add_constrained_variable(o, MOI.Interval(0.0, 5.0))
    MOI.set(o, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    func = MOI.ScalarAffineFunction(
        [MOI.ScalarAffineTerm(1.0, x1), MOI.ScalarAffineTerm(2.0, x2)],
        0.0,
    )
    @test MOI.supports_constraint(o, typeof(func), MOI.Interval{Float64})
    MOI.set(o, MOI.ObjectiveFunction{typeof(func)}(), func)
    @test MOI.get(
        o,
        MOI.NumberOfConstraints{
            MOI.ScalarAffineFunction{Float64},
            MOI.Interval{Float64},
        }(),
    ) == 0
    MOI.add_constraint(
        o,
        MOI.ScalarAffineFunction(
            [MOI.ScalarAffineTerm(1.0, x1), MOI.ScalarAffineTerm(1.0, x2)],
            0.0,
        ),
        MOI.Interval(0.0, 7.5),
    )
    @test MOI.get(
        o,
        MOI.NumberOfConstraints{
            MOI.ScalarAffineFunction{Float64},
            MOI.Interval{Float64},
        }(),
    ) == 1
    MOI.optimize!(o)
    @test MOI.get(o, MOI.ObjectiveValue()) ≈ 12.5
    @test MOI.get(o, MOI.SimplexIterations()) > 0
end

function test_Variable_names()
    o = HiGHS.Optimizer()
    MOI.Test.variablenames(o, CONFIG)
    MOI.empty!(o)
    y = MOI.add_variable(o)
    MOI.set(o, MOI.VariableName(), y, "y")
    y2 = MOI.get(o, MOI.VariableIndex, "y")
    @test y == y2
    @test MOI.get(o, MOI.VariableIndex, "y0") === nothing
end

function test_HiGHS_custom_options()
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
    @test MOI.supports(o, MOI.RawParameter("wrong_param")) == false
    @test_throws MOI.UnsupportedAttribute MOI.get(o, MOI.RawParameter("wrong_param"))
    for v in [false, 1, 1.0, "A"]
        @test_throws(
            MOI.UnsupportedAttribute,
            MOI.set(o, MOI.RawParameter("wrong_param"), v)
        )
    end
end

function test_Model_empty()
    o = HiGHS.Optimizer()
    @test MOI.is_empty(o)
    MOI.add_variable(o)
    @test !MOI.is_empty(o)
    MOI.empty!(o)
    MOI.set(
        o,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction{Float64}([], 0.0),
    )
    @test MOI.is_empty(o)
    MOI.set(
        o,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction{Float64}([], 3.0),
    )
    @test !MOI.is_empty(o)
    MOI.set(o, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
    @test MOI.is_empty(o)
    @test MOI.get(
        o,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
    ) ≈ MOI.ScalarAffineFunction{Float64}([], 0.0)
    x = MOI.add_variable(o)
    return MOI.set(
        o,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm(1.0, x)], 0.0),
    )
end

function test_show()
    model = HiGHS.Optimizer()
    @test sprint(show, model) == "A HiGHS model with 0 columns and 0 rows."
    return
end

function test_options()
    model = HiGHS.Optimizer()
    options = [
        "write_solution_to_file",   # Bool
        "simplex_strategy",         # Cint
        "time_limit",               # Cdouble
        "presolve",                 # String
    ]
    for key in options
        v = MOI.get(model, MOI.RawParameter(key))
        MOI.set(model, MOI.RawParameter(key), v)
        v2 = MOI.get(model, MOI.RawParameter(key))
        @test v == v2
    end
    return
end

function test_option_invalid()
    model = HiGHS.Optimizer()
    @test_throws(
        ErrorException(
            "Encountered an error in HiGHS: Check the log for details."
        ),
        MOI.set(model, MOI.RawParameter("time_limit"), -1.0),
    )
    return
end

function test_option_invalid()
    model = HiGHS.Optimizer()
    MOI.supports(model, MOI.RawParameter(:presolve)) == false
    @test_throws MOI.UnsupportedAttribute MOI.get(model, MOI.RawParameter(:presolve))
    return
end

function test_option_unknown_option()
    model = HiGHS.Optimizer()
    err = ErrorException(
        "Encountered an error in HiGHS: Check the log for details."
    )
    @test_throws err MOI.set(model, MOI.RawParameter("write_solution_to_file"), 1)
    @test_throws err MOI.set(model, MOI.RawParameter("simplex_strategy"), "on")
    @test_throws err MOI.set(model, MOI.RawParameter("time_limit"), 1)
    @test_throws err MOI.set(model, MOI.RawParameter("presolve"), 1)
    return
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith(string(name), "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

end

TestMOIHighs.runtests()
