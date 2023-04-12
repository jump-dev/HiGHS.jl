# Copyright (c) 2019 Mathieu Besançon, Oscar Dowson, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestMOIHighs

import HiGHS
using Test

const MOI = HiGHS.MOI

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$name" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

const _EXPLICIT_METHOD_FAILURES = [
    "test_objective_qp_ObjectiveFunction_edge_cases",
    "test_objective_qp_ObjectiveFunction_zero_ofdiag",
    "test_quadratic_duplicate_terms",
    "test_quadratic_integration",
    "test_quadratic_nonhomogeneous",
    "test_linear_Semicontinuous_integration",
    "test_linear_Semiinteger_integration",
]

function test_runtests()
    model = MOI.Bridges.full_bridge_optimizer(HiGHS.Optimizer(), Float64)
    MOI.set(model, MOI.Silent(), true)
    # Slightly loosen tolerances, particularly for QP tests
    MOI.Test.runtests(model, MOI.Test.Config(; atol = 1e-7))
    return
end

function test_runtests_cache()
    model = MOI.Bridges.full_bridge_optimizer(
        MOI.Utilities.CachingOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
            HiGHS.Optimizer(),
        ),
        Float64,
    )
    MOI.set(model, MOI.Silent(), true)
    # Slightly loosen tolerances, particularly for QP tests
    MOI.Test.runtests(model, MOI.Test.Config(; atol = 1e-7))
    return
end

function test_runtests_simplex()
    model = MOI.Bridges.full_bridge_optimizer(HiGHS.Optimizer(), Float64)
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, MOI.RawOptimizerAttribute("solver"), "simplex")
    for presolve in ("on", "off")
        MOI.set(model, MOI.RawOptimizerAttribute("presolve"), presolve)
        MOI.Test.runtests(
            model,
            MOI.Test.Config();
            exclude = _EXPLICIT_METHOD_FAILURES,
        )
    end
    return
end

function test_runtests_ipm()
    model = MOI.Bridges.full_bridge_optimizer(HiGHS.Optimizer(), Float64)
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, MOI.RawOptimizerAttribute("solver"), "ipm")
    MOI.Test.runtests(
        model,
        MOI.Test.Config();
        exclude = _EXPLICIT_METHOD_FAILURES,
    )
    return
end

function test_runtests_ipm_no_presolve()
    model = MOI.Bridges.full_bridge_optimizer(HiGHS.Optimizer(), Float64)
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, MOI.RawOptimizerAttribute("solver"), "ipm")
    MOI.set(model, MOI.RawOptimizerAttribute("presolve"), "off")
    MOI.Test.runtests(
        model,
        MOI.Test.Config(),
        exclude = String[
            # Termination status is OTHER_ERROR
            "test_conic_linear_INFEASIBLE",
            "test_conic_linear_INFEASIBLE_2",
            _EXPLICIT_METHOD_FAILURES...,
        ],
    )
    return
end

function test_SolverName()
    @test MOI.get(HiGHS.Optimizer(), MOI.SolverName()) == "HiGHS"
    return
end

function test_attributes()
    model = HiGHS.Optimizer()
    @test MOI.get(model, MOI.SolverName()) == "HiGHS"
    @test MOI.get(model, MOI.TimeLimitSec()) > 10000
    MOI.set(model, MOI.TimeLimitSec(), 500)
    @test MOI.get(model, MOI.TimeLimitSec()) == 500.0
    @test MOI.get(model, MOI.RawSolver()) == model
    return
end

function test_MOI_variable_count_and_empty()
    model = HiGHS.Optimizer()
    @test MOI.get(model, MOI.NumberOfVariables()) == 0
    x1 = MOI.add_variable(model)
    @test MOI.get(model, MOI.NumberOfVariables()) == 1
    @test MOI.supports_constraint(
        model,
        MOI.VariableIndex,
        MOI.Interval{Float64},
    )
    x2, _ = MOI.add_constrained_variable(model, MOI.Interval(0.0, 1.0))
    @test MOI.get(model, MOI.NumberOfVariables()) == 2
    MOI.empty!(model)
    @test MOI.get(model, MOI.NumberOfVariables()) == 0
end

function test_HiGHS_custom_options()
    model = HiGHS.Optimizer()
    @test MOI.supports(model, MOI.RawOptimizerAttribute("solver"))
    @test MOI.get(model, MOI.RawOptimizerAttribute("solver")) == "choose"
    MOI.set(model, MOI.RawOptimizerAttribute("solver"), "simplex")
    @test MOI.get(model, MOI.RawOptimizerAttribute("solver")) == "simplex"
    @test MOI.get(model, MOI.RawOptimizerAttribute("output_flag")) == true
    MOI.set(model, MOI.RawOptimizerAttribute("output_flag"), false)
    @test MOI.get(model, MOI.RawOptimizerAttribute("output_flag")) == false
    @test MOI.get(model, MOI.RawOptimizerAttribute("time_limit")) > 1000
    MOI.set(model, MOI.RawOptimizerAttribute("time_limit"), 1000.0)
    @test MOI.get(model, MOI.RawOptimizerAttribute("time_limit")) == 1000.0
    # unsupported test
    @test MOI.supports(model, MOI.RawOptimizerAttribute("wrong_param")) == false
    @test_throws(
        MOI.UnsupportedAttribute,
        MOI.get(model, MOI.RawOptimizerAttribute("wrong_param")),
    )
    for v in [false, 1, 1.0, "A"]
        @test_throws(
            MOI.UnsupportedAttribute,
            MOI.set(model, MOI.RawOptimizerAttribute("wrong_param"), v)
        )
    end
    return
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
        v = MOI.get(model, MOI.RawOptimizerAttribute(key))
        MOI.set(model, MOI.RawOptimizerAttribute(key), v)
        v2 = MOI.get(model, MOI.RawOptimizerAttribute(key))
        @test v == v2
    end
    return
end

function test_option_invalid()
    model = HiGHS.Optimizer()
    @test_throws(
        ErrorException(
            "Encountered an error in HiGHS: Check the log for details.",
        ),
        MOI.set(model, MOI.RawOptimizerAttribute("time_limit"), -1.0),
    )
    return
end

function test_option_unknown_option()
    model = HiGHS.Optimizer()
    param = MOI.RawOptimizerAttribute("write_solution_to_file")
    err = MOI.SetAttributeNotAllowed(
        param,
        "\n\nInvalid value `1::Int64` for option \"write_solution_to_file\", expected a value of type `Bool`.\n\n",
    )
    @test_throws(err, MOI.set(model, param, 1))
    param = MOI.RawOptimizerAttribute("simplex_strategy")
    err = MOI.SetAttributeNotAllowed(
        param,
        "\n\nInvalid value `on::String` for option \"simplex_strategy\", expected a value of type `Integer`.\n\n",
    )
    @test_throws(err, MOI.set(model, param, "on"))
    param = MOI.RawOptimizerAttribute("time_limit")
    err = MOI.SetAttributeNotAllowed(
        param,
        "\n\nInvalid value `1::Int64` for option \"time_limit\", expected a value of type `AbstractFloat`.\n\n",
    )
    @test_throws err MOI.set(model, param, 1)
    param = MOI.RawOptimizerAttribute("presolve")
    err = MOI.SetAttributeNotAllowed(
        param,
        "\n\nInvalid value `1::Int64` for option \"presolve\", expected a value of type `String`.\n\n",
    )
    @test_throws err MOI.set(model, param, 1)
    return
end

function test_copy_to()
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
    list = MOI.get(dest, MOI.ListOfConstraintTypesPresent())
    @test length(list) == 8
    for S in (
        MOI.GreaterThan{Float64},
        MOI.LessThan{Float64},
        MOI.EqualTo{Float64},
        MOI.Interval{Float64},
    )
        @test (MOI.VariableIndex, S) in list
        @test (MOI.ScalarAffineFunction{Float64}, S) in list
    end
    MOI.optimize!(dest)
    @test MOI.get(dest, MOI.ObjectiveValue()) == 0.5 + 1.0 + 2.0 + 2.5
    return
end

function _knapsack_model(; mip::Bool, solver::String)
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.RawOptimizerAttribute("presolve"), "off")
    MOI.set(model, MOI.RawOptimizerAttribute("solver"), solver)
    N = 30
    x = MOI.add_variables(model, N)
    if mip
        MOI.add_constraints(model, x, MOI.ZeroOne())
    end
    item_weights, item_values = rand(N), rand(N)
    MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(item_weights, x), 0.0),
        MOI.LessThan(10.0),
    )
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(item_values, x), 0.0),
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    return model
end

function test_RelativeGap()
    model = _knapsack_model(mip = true, solver = "choose")
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.RelativeGap()), 0.0; atol = 1e-14)
    return
end

function test_SimplexIterations_BarrierIterations()
    model = _knapsack_model(mip = false, solver = "simplex")
    @test MOI.get(model, MOI.SimplexIterations()) == 0
    @test MOI.get(model, MOI.BarrierIterations()) == 0
    MOI.optimize!(model)
    @test MOI.get(model, MOI.SimplexIterations()) > 0
    @test MOI.get(model, MOI.BarrierIterations()) == 0
    model = _knapsack_model(mip = false, solver = "ipm")
    MOI.optimize!(model)
    # Not == 0 because HiGHS will use Simplex to clean-up occasionally
    @test MOI.get(model, MOI.SimplexIterations()) >= 0
    @test MOI.get(model, MOI.BarrierIterations()) > 0
end

function test_NodeCount()
    model = _knapsack_model(mip = true, solver = "choose")
    @test MOI.get(model, MOI.NodeCount()) == 0
    MOI.optimize!(model)
    @test MOI.get(model, MOI.NodeCount()) > 0
    return
end

function test_option_nothing()
    model = HiGHS.Optimizer()
    @test_throws(
        MOI.SetAttributeNotAllowed,
        MOI.set(model, MOI.RawOptimizerAttribute("presolve"), nothing),
    )
    return
end

function test_copy_to_names()
    dest = HiGHS.Optimizer()
    src = MOI.Utilities.Model{Float64}()
    MOI.Utilities.loadfromstring!(src, "variables: x\nc: 2.0 * x <= 1.0")
    _ = MOI.copy_to(dest, src)
    @test MOI.get(dest, MOI.VariableIndex, "x") isa MOI.VariableIndex
    F, S = MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}
    @test MOI.get(dest, MOI.ConstraintIndex, "c") isa MOI.ConstraintIndex{F,S}
    return
end

function test_copy_to_sets()
    for (s_set, set) in [
        "Semicontinuous(1.0, 2.0)" => MOI.Semicontinuous(1.0, 2.0),
        "Semiinteger(2.0, 3.0)" => MOI.Semiinteger(2.0, 3.0),
        "ZeroOne()" => MOI.ZeroOne(),
        "Integer()" => MOI.Integer(),
    ]
        dest = HiGHS.Optimizer()
        src = MOI.Utilities.Model{Float64}()
        MOI.Utilities.loadfromstring!(src, "variables: x\nc: x in $s_set")
        _ = MOI.copy_to(dest, src)
        x = MOI.get(dest, MOI.VariableIndex, "x")
        ci = MOI.ConstraintIndex{MOI.VariableIndex,typeof(set)}(x.value)
        @test MOI.get(dest, MOI.ConstraintSet(), ci) == set
        @test MOI.get(dest, MOI.ListOfConstraintTypesPresent()) ==
              [(MOI.VariableIndex, typeof(set))]
        if set isa MOI.Semicontinuous || set isa MOI.Semiinteger
            @test_throws(
                MOI.UpperBoundAlreadySet{typeof(set),MOI.LessThan{Float64}},
                MOI.add_constraint(dest, x, MOI.LessThan(1.0)),
            )
            @test_throws(
                MOI.LowerBoundAlreadySet{typeof(set),MOI.GreaterThan{Float64}},
                MOI.add_constraint(dest, x, MOI.GreaterThan(1.0)),
            )
        end
    end
    return
end

function test_delete_vector()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 5)
    MOI.add_constraint.(model, x, MOI.GreaterThan(0.0))
    c = MOI.add_constraint.(model, 1.0 .* x, MOI.GreaterThan.(1.0:5.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    f = sum(Float64(i) * x[i] for i in 1:5)
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    F, S = MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}
    @test MOI.get(model, MOI.NumberOfConstraints{F,S}()) == 5
    MOI.delete(model, c)
    @test MOI.get(model, MOI.NumberOfConstraints{F,S}()) == 0
    c = MOI.add_constraint.(model, 1.0 .* x, MOI.GreaterThan.(1.0:5.0))
    MOI.delete(model, [c[2], c[4]])
    @test MOI.get(model, MOI.NumberOfConstraints{F,S}()) == 3
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ObjectiveValue()) ≈ 1^2 + 3^2 + 5^2
    MOI.set.(model, MOI.ConstraintSet(), c[1:2:5], MOI.GreaterThan(1.0))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ObjectiveValue()) ≈ 1 + 3 + 5
    return
end

function test_option_type()
    for x in ["1", 1.0, 1, true]
        k = HiGHS._highs_option_type(x)
        T = HiGHS._type_for_highs_option(k)
        @test x isa T
    end
end

function test_quadratic_sets_objective()
    model = HiGHS.Optimizer()
    MOI.Utilities.loadfromstring!(
        model,
        """
        variables: x
        minobjective: 1.0 * x * x
        """,
    )
    attr = MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}()
    @test attr in MOI.get(model, MOI.ListOfModelAttributesSet())
    return
end

function test_dual_issue_157()
    model = HiGHS.Optimizer()
    x, y = MOI.add_variables(model, 2)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    f = 2.0 * x * x + 1.0 * x * y + 1.0 * y * y + x + y
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    c = MOI.add_constraint(model, 1.0 * x + y, MOI.LessThan(-1.0))
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ConstraintDual(), c), -0.75; atol = 1e-4)
    @test_throws(
        MOI.GetAttributeNotAllowed,
        MOI.get(model, MOI.ConstraintBasisStatus(), c),
    )
    @test_throws(
        MOI.GetAttributeNotAllowed,
        MOI.get(model, MOI.VariableBasisStatus(), x),
    )
    return
end

end

TestMOIHighs.runtests()
