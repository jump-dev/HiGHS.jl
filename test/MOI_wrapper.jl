module TestMOIHighs

import HiGHS
using Test

const MOI = HiGHS.MOI

function test_basic_constraint_tests(model, config)
    return MOI.Test.basic_constraint_tests(model, config)
end

function test_unittest(model, config)
    return MOI.Test.unittest(
        model,
        config,
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

function test_modificationtest(model, config)
    return MOI.Test.modificationtest(model, config)
end

function test_contlineartest(model, config)
    excludes = [
        # VariablePrimalStart not supported.
        "partial_start",
    ]
    if MOI.get(model, MOI.RawParameter("solver")) == "ipm"
        # TODO(odow): investigate
        push!(excludes, "linear8b")
        push!(excludes, "linear8c")
    end
    return MOI.Test.contlineartest(model, config, excludes)
end

function test_lintest(model, config)
    return MOI.Test.lintest(model, config)
end

function test_SolverName(model, ::Any)
    @test MOI.get(model, MOI.SolverName()) == "HiGHS"
end

function test_default_objective(model, ::Any)
    return MOI.Test.default_objective_test(model)
end

function test_default_status_test(model, ::Any)
    return MOI.Test.default_status_test(model)
end

function test_nametest(model, ::Any)
    return MOI.Test.nametest(model)
end

function test_validtest(model, ::Any)
    return MOI.Test.validtest(model)
end

function test_emptytest(model, ::Any)
    return MOI.Test.emptytest(model)
end

function test_orderedindicestest(model, ::Any)
    return MOI.Test.orderedindicestest(model)
end

function test_copytest(model, ::Any)
    return MOI.Test.copytest(
        model,
        MOI.Bridges.full_bridge_optimizer(HiGHS.Optimizer(), Float64),
    )
end

function test_scalar_function_constant_not_zero(model, ::Any)
    return MOI.Test.scalar_function_constant_not_zero(model)
end

function test_supports_constrainttest(model, ::Any)
    # supports_constrainttest needs VectorOfVariables-in-Zeros,
    # MOI.Test.supports_constrainttest(HiGHS.Optimizer(), Float64, Float32)
    # but supports_constrainttest is broken via bridges:
    MOI.empty!(model)
    MOI.add_variable(model)
    @test MOI.supports_constraint(
        model,
        MOI.SingleVariable,
        MOI.EqualTo{Float64},
    )
    @test MOI.supports_constraint(
        model,
        MOI.ScalarAffineFunction{Float64},
        MOI.EqualTo{Float64},
    )
    # This test is broken for some reason:
    @test_broken !MOI.supports_constraint(
        model,
        MOI.ScalarAffineFunction{Int},
        MOI.EqualTo{Float64},
    )
    @test !MOI.supports_constraint(
        model,
        MOI.ScalarAffineFunction{Int},
        MOI.EqualTo{Int},
    )
    @test !MOI.supports_constraint(model, MOI.SingleVariable, MOI.EqualTo{Int})
    @test MOI.supports_constraint(model, MOI.VectorOfVariables, MOI.Zeros)
    @test !MOI.supports_constraint(
        model,
        MOI.VectorOfVariables,
        MOI.EqualTo{Float64},
    )
    @test !MOI.supports_constraint(model, MOI.SingleVariable, MOI.Zeros)
    @test !MOI.supports_constraint(
        model,
        MOI.VectorOfVariables,
        MOI.Test.UnknownVectorSet,
    )
end

function test_set_lower_bound_twice(::Any, ::Any)
    return MOI.Test.set_lower_bound_twice(HiGHS.Optimizer(), Float64)
end

function test_set_upper_bound_twice(::Any, ::Any)
    return MOI.Test.set_upper_bound_twice(HiGHS.Optimizer(), Float64)
end

function test_Attributes(::Any, ::Any)
    model = HiGHS.Optimizer()
    @test MOI.get(model, MOI.SolverName()) == "HiGHS"
    @test MOI.get(model, MOI.TimeLimitSec()) > 10000
    MOI.set(model, MOI.TimeLimitSec(), 500)
    @test MOI.get(model, MOI.TimeLimitSec()) == 500.0
    @test MOI.get(model, MOI.RawSolver()) == model
end

function test_MOI_variable_count_and_empty(model, ::Any)
    MOI.empty!(model)
    @test MOI.get(model, MOI.NumberOfVariables()) == 0
    x1 = MOI.add_variable(model)
    @test MOI.get(model, MOI.NumberOfVariables()) == 1
    @test MOI.supports_constraint(
        model,
        MOI.SingleVariable,
        MOI.Interval{Float64},
    )
    x2, _ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0))
    @test MOI.get(model, MOI.NumberOfVariables()) == 2
    MOI.empty!(model)
    @test MOI.get(model, MOI.NumberOfVariables()) == 0
end

function test_Getting_objective_value(model, ::Any)
    MOI.empty!(model)
    MOI.set(model, MOI.Silent(), true)
    (x, _) = MOI.add_constrained_variable(model, MOI.Interval(-3.0, 6.0))
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x)], 0.0),
    )
    @test MOI.get(model, MOI.ObjectiveSense()) == MOI.FEASIBILITY_SENSE
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    @test MOI.get(model, MOI.ObjectiveSense()) == MOI.MIN_SENSE
    @test MOI.get(model, MOI.ResultCount()) == 0
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ResultCount()) == 1
    @test MOI.get(model, MOI.ObjectiveValue()) ≈ -3
end

function test_Max_in_box(model, ::Any)
    MOI.empty!(model)
    MOI.set(model, MOI.Silent(), true)
    @test MOI.get(model, MOI.ResultCount()) == 0
    (x, _) = MOI.add_constrained_variable(model, MOI.Interval(-3.0, 6.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(2.0, x)], 0.0),
    )
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ObjectiveValue()) ≈ 2 * 6
    obj_func = MOI.get(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
    )
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test obj_func ≈
          MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(2.0, x)], 0.0)
end

function test_Objective_function_obtained_from_model_corresponds(model, ::Any)
    MOI.empty!(model)
    (x1, _) = MOI.add_constrained_variable(model, MOI.Interval(-3.0, 6.0))
    (x2, _) = MOI.add_constrained_variable(model, MOI.Interval(1.0, 2.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(
            MOI.ScalarAffineTerm.([2.0, -1.0], [x1, x2]),
            0.0,
        ),
    )
    F = MOI.get(model, MOI.ObjectiveFunctionType())
    @test F <: MOI.ScalarAffineFunction{Float64}
    obj_func = MOI.get(model, MOI.ObjectiveFunction{F}())
    @test MOI.supports(model, MOI.ObjectiveFunction{F}())
    @test all(MOI.get(model, MOI.ListOfVariableIndices()) .== [x1, x2])
    @test obj_func ≈ MOI.ScalarAffineFunction(
        [MOI.ScalarAffineTerm(2.0, x1), MOI.ScalarAffineTerm(-1.0, x2)],
        0.0,
    )
    MOI.set(model, MOI.ObjectiveFunction{F}(), obj_func)
    obj_func = MOI.get(model, MOI.ObjectiveFunction{F}())
    @test obj_func ≈ MOI.ScalarAffineFunction(
        [MOI.ScalarAffineTerm(2.0, x1), MOI.ScalarAffineTerm(-1.0, x2)],
        0.0,
    )
    obj_func.terms[1] = MOI.ScalarAffineTerm(3.0, x1)
    MOI.set(model, MOI.ObjectiveFunction{F}(), obj_func)
    obj_func = MOI.get(model, MOI.ObjectiveFunction{F}())
    @test obj_func ≈ MOI.ScalarAffineFunction(
        [MOI.ScalarAffineTerm(3.0, x1), MOI.ScalarAffineTerm(-1.0, x2)],
        0.0,
    )
end

function test_Constrained_variable_equivalent_to_add_constraint(model, ::Any)
    MOI.empty!(model)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    _ = MOI.add_constraint(
        model,
        MOI.SingleVariable(x),
        MOI.Interval(-3.0, 6.0),
    )
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x)], 0.0),
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    @test MOI.get(model, MOI.ResultCount()) == 0
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ResultCount()) == 1
    @test MOI.get(model, MOI.ObjectiveValue()) ≈ -3
end

function test_Constant_in_objective_function(model, ::Any)
    MOI.empty!(model)
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    _ = MOI.add_constraint(
        model,
        MOI.SingleVariable(x),
        MOI.Interval(-3.0, 6.0),
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj_func = MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x)], 3.0)
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_func)}(), obj_func)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ResultCount()) == 1
    @test MOI.get(model, MOI.ObjectiveValue()) ≈ 0
    obj_func = MOI.get(model, MOI.ObjectiveFunction{typeof(obj_func)}())
    @test MOI.constant(obj_func) ≈ 3
    MOI.set(model, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
    obj_func = MOI.get(model, MOI.ObjectiveFunction{typeof(obj_func)}())
    @test MOI.constant(obj_func) ≈ 0
    @test isempty(obj_func.terms)
end

function test_Linear_constraints(model, ::Any)
    MOI.empty!(model)
    # max x1 + 2x2
    # st 0 <= x{1,2} <= 5
    # 0 <= x1 + x2 <= 7.5
    MOI.set(model, MOI.Silent(), true)
    (x1, _) = MOI.add_constrained_variable(model, MOI.Interval(0.0, 5.0))
    (x2, _) = MOI.add_constrained_variable(model, MOI.Interval(0.0, 5.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    func = MOI.ScalarAffineFunction(
        [MOI.ScalarAffineTerm(1.0, x1), MOI.ScalarAffineTerm(2.0, x2)],
        0.0,
    )
    @test MOI.supports_constraint(model, typeof(func), MOI.Interval{Float64})
    MOI.set(model, MOI.ObjectiveFunction{typeof(func)}(), func)
    @test MOI.get(
        model,
        MOI.NumberOfConstraints{
            MOI.ScalarAffineFunction{Float64},
            MOI.Interval{Float64},
        }(),
    ) == 0
    MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(
            [MOI.ScalarAffineTerm(1.0, x1), MOI.ScalarAffineTerm(1.0, x2)],
            0.0,
        ),
        MOI.Interval(0.0, 7.5),
    )
    @test MOI.get(
        model,
        MOI.NumberOfConstraints{
            MOI.ScalarAffineFunction{Float64},
            MOI.Interval{Float64},
        }(),
    ) == 1
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ObjectiveValue()) ≈ 12.5
end

function test_Variable_names(model, config)
    MOI.empty!(model)
    MOI.Test.variablenames(model, config)
    MOI.empty!(model)
    y = MOI.add_variable(model)
    MOI.set(model, MOI.VariableName(), y, "y")
    y2 = MOI.get(model, MOI.VariableIndex, "y")
    @test y == y2
    @test MOI.get(model, MOI.VariableIndex, "y0") === nothing
end

function test_HiGHS_custom_options(::Any, ::Any)
    model = HiGHS.Optimizer()
    @test MOI.supports(model, MOI.RawParameter("solver"))
    @test MOI.get(model, MOI.RawParameter("solver")) == "choose"
    MOI.set(model, MOI.RawParameter("solver"), "simplex")
    @test MOI.get(model, MOI.RawParameter("solver")) == "simplex"
    @test MOI.get(model, MOI.RawParameter("output_flag")) == true
    MOI.set(model, MOI.RawParameter("output_flag"), false)
    @test MOI.get(model, MOI.RawParameter("output_flag")) == false
    @test MOI.get(model, MOI.RawParameter("time_limit")) > 1000
    MOI.set(model, MOI.RawParameter("time_limit"), 1000.0)
    @test MOI.get(model, MOI.RawParameter("time_limit")) == 1000.0
    # unsupported test
    @test MOI.supports(model, MOI.RawParameter("wrong_param")) == false
    @test_throws MOI.UnsupportedAttribute MOI.get(
        model,
        MOI.RawParameter("wrong_param"),
    )
    for v in [false, 1, 1.0, "A"]
        @test_throws(
            MOI.UnsupportedAttribute,
            MOI.set(model, MOI.RawParameter("wrong_param"), v)
        )
    end
end

function test_Model_empty(model, ::Any)
    MOI.empty!(model)
    @test MOI.is_empty(model)
    MOI.add_variable(model)
    @test !MOI.is_empty(model)
    MOI.empty!(model)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction{Float64}([], 0.0),
    )
    @test MOI.is_empty(model)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction{Float64}([], 3.0),
    )
    @test !MOI.is_empty(model)
    MOI.set(model, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
    @test MOI.is_empty(model)
    @test MOI.get(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
    ) ≈ MOI.ScalarAffineFunction{Float64}([], 0.0)
    x = MOI.add_variable(model)
    return MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction{Float64}([MOI.ScalarAffineTerm(1.0, x)], 0.0),
    )
end

function test_show(::Any, ::Any)
    model = HiGHS.Optimizer()
    @test sprint(show, model) == "A HiGHS model with 0 columns and 0 rows."
    return
end

function test_options(model, ::Any)
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

function test_option_invalid(model, ::Any)
    @test_throws(
        ErrorException(
            "Encountered an error in HiGHS: Check the log for details.",
        ),
        MOI.set(model, MOI.RawParameter("time_limit"), -1.0),
    )
    return
end

function test_option_invalid(model, ::Any)
    MOI.supports(model, MOI.RawParameter(:presolve)) == false
    @test_throws MOI.UnsupportedAttribute MOI.get(
        model,
        MOI.RawParameter(:presolve),
    )
    return
end

function test_option_unknown_option(model, ::Any)
    err = ErrorException(
        "Encountered an error in HiGHS: Check the log for details.",
    )
    @test_throws err MOI.set(
        model,
        MOI.RawParameter("write_solution_to_file"),
        1,
    )
    @test_throws err MOI.set(model, MOI.RawParameter("simplex_strategy"), "on")
    @test_throws err MOI.set(model, MOI.RawParameter("time_limit"), 1)
    @test_throws err MOI.set(model, MOI.RawParameter("presolve"), 1)
    return
end

function test_copy_to(::Any, ::Any)
    src = MOI.Utilities.Model{Float64}()
    MOI.Utilities.loadfromstring!(
        src,
        """
variables: w, x, y, z
minobjective: w + x + y + z
c1: w <= 1.0
c2: w >= 0.5
c3: x in MathOptInterface.Interval(1.0, 2.0)
c4: y == 2.0
c5: y + z >= 4.5
c6: x + y <= 3.0
c7: w + x == 1.5
c8: w + x in MathOptInterface.Interval(1.0, 2.0)
""",
    )
    dest = HiGHS.Optimizer()
    MOI.copy_to(dest, src)
    @test MOI.get(dest, MOI.NumberOfVariables()) == 4
    list = MOI.get(dest, MOI.ListOfConstraints())
    @test length(list) == 8
    for S in (
        MOI.GreaterThan{Float64},
        MOI.LessThan{Float64},
        MOI.EqualTo{Float64},
        MOI.Interval{Float64},
    )
        @test (MOI.SingleVariable, S) in list
        @test (MOI.ScalarAffineFunction{Float64}, S) in list
    end
    MOI.optimize!(dest)
    @test MOI.get(dest, MOI.ObjectiveValue()) == 0.5 + 1.0 + 2.0 + 2.5
    return
end

function optimizer(; solver::String = "simplex")
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    # Turn off presolve for simplex to get infeas_certificates.
    presolve = solver == "simplex" ? "off" : "on"
    MOI.set(model, MOI.RawParameter("presolve"), presolve)
    MOI.set(model, MOI.RawParameter("solver"), solver)
    return MOI.Bridges.full_bridge_optimizer(model, Float64)
end

function runtests()
    config = Dict(
        "simplex" => MOI.Test.TestConfig(basis = true),
        "ipm" =>
            MOI.Test.TestConfig(basis = true, infeas_certificates = false),
    )
    @testset "$(solver)" for solver in ["simplex", "ipm"]
        model = optimizer(solver = solver)
        for name in names(@__MODULE__; all = true)
            if startswith(string(name), "test_")
                @testset "$(name)" begin
                    getfield(@__MODULE__, name)(model, config[solver])
                end
            end
        end
    end
end

end

TestMOIHighs.runtests()
