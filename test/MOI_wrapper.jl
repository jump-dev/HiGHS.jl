# Copyright (c) 2019 Mathieu Besançon, Oscar Dowson, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestMOIHighs

using Test

import HiGHS
import MathOptInterface as MOI

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

function test_runtests()
    model = MOI.Bridges.full_bridge_optimizer(HiGHS.Optimizer(), Float64)
    MOI.set(model, MOI.Silent(), true)
    # Turn presolve off so that we generate infeasibility certificates. This is
    # a temporary work-around until we fix this upstream in HiGHS.
    MOI.set(model, MOI.RawOptimizerAttribute("presolve"), "off")
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

function test_SolverName()
    @test MOI.get(HiGHS.Optimizer(), MOI.SolverName()) == "HiGHS"
    return
end

function test_attributes()
    model = HiGHS.Optimizer()
    @test MOI.get(model, MOI.SolverName()) == "HiGHS"
    @test MOI.get(model, MOI.TimeLimitSec()) === nothing
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
    @test 0 <= MOI.get(model, MOI.RelativeGap()) <= 1e-4
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

function test_attribute_TimeLimitSec()
    model = HiGHS.Optimizer()
    @test MOI.supports(model, MOI.TimeLimitSec())
    @test MOI.get(model, MOI.TimeLimitSec()) === nothing
    MOI.set(model, MOI.TimeLimitSec(), 0.0)
    @test MOI.get(model, MOI.TimeLimitSec()) == 0.0
    MOI.set(model, MOI.TimeLimitSec(), nothing)
    @test MOI.get(model, MOI.TimeLimitSec()) === nothing
    MOI.set(model, MOI.TimeLimitSec(), 1.0)
    @test MOI.get(model, MOI.TimeLimitSec()) == 1.0
    return
end

function test_copy_to_bug_172()
    model = MOI.Utilities.Model{Float64}()
    x = MOI.add_variable(model)
    F = MOI.ScalarAffineFunction{Float64}
    c1 = MOI.add_constraint(model, 2.0 * x, MOI.GreaterThan(0.0))
    c2 = MOI.add_constraint(model, zero(F), MOI.GreaterThan(0.0))
    c3 = MOI.add_constraint(model, 1.0 * x, MOI.EqualTo(1.0))
    h = HiGHS.Optimizer()
    MOI.set(h, MOI.Silent(), true)
    index_map = MOI.copy_to(h, model)
    y = index_map[x]
    @test MOI.get(h, MOI.ConstraintFunction(), index_map[c1]) ≈ 2.0 * y
    @test MOI.get(h, MOI.ConstraintFunction(), index_map[c2]) ≈ zero(F)
    @test MOI.get(h, MOI.ConstraintFunction(), index_map[c3]) ≈ 1.0 * y
    @test MOI.get(h, MOI.ConstraintSet(), index_map[c1]) == MOI.GreaterThan(0.0)
    @test MOI.get(h, MOI.ConstraintSet(), index_map[c2]) == MOI.GreaterThan(0.0)
    @test MOI.get(h, MOI.ConstraintSet(), index_map[c3]) == MOI.EqualTo(1.0)
    MOI.optimize!(h)
    @test MOI.get(h, MOI.TerminationStatus()) == MOI.OPTIMAL
    return
end

function test_relax_integrality_after_solve()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    MOI.add_constraint(model, x, MOI.LessThan(2.0))
    c = MOI.add_constraint(model, x, MOI.ZeroOne())
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    f = 1.0 * x
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.VariablePrimal(), x), 1.0; atol = 1e-6)
    MOI.delete(model, c)
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.VariablePrimal(), x), 2.0; atol = 1e-6)
    return
end

function test_quadratic_modification_from_affine()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    MOI.add_constraint(model, x, MOI.GreaterThan(2.0))
    f = 2.0 * x + 1.0
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    attr = MOI.ObjectiveFunction{typeof(f)}()
    MOI.set(model, attr, f)
    MOI.optimize!(model)
    @test isapprox(MOI.get(model, MOI.ObjectiveValue()), 5, atol = 1e-5)
    F = MOI.ScalarQuadraticFunction{Float64}
    attr = MOI.ObjectiveFunction{F}()
    MOI.modify(model, attr, MOI.ScalarQuadraticCoefficientChange(x, x, 3.0))
    MOI.optimize!(model)
    @test isapprox(MOI.get(model, MOI.ObjectiveValue()), 11, atol = 1e-5)
    return
end

function test_quadratic_off_diagonal_modification()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 2)
    MOI.add_constraint.(model, x, MOI.GreaterThan.([2.0, 3.0]))
    f = 4.0 * x[1] * x[1] + 2.0 * x[1] * x[2] + 2.0 * x[2] * x[2]
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    attr = MOI.ObjectiveFunction{typeof(f)}()
    MOI.set(model, attr, f)
    MOI.optimize!(model)
    a = MOI.get(model, MOI.VariablePrimal(), x)
    y = 0.5 * a' * [8 2; 2 4] * a
    @test isapprox(MOI.get(model, MOI.ObjectiveValue()), y, atol = 1e-5)
    MOI.modify(
        model,
        attr,
        MOI.ScalarQuadraticCoefficientChange(x[1], x[2], -1.0),
    )
    MOI.optimize!(model)
    a = MOI.get(model, MOI.VariablePrimal(), x)
    y = 0.5 * a' * [8 -1; -1 4] * a
    @test isapprox(MOI.get(model, MOI.ObjectiveValue()), y, atol = 1e-5)
    return
end

function test_quadratic_diagonal_modification()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    MOI.add_constraint(model, x, MOI.GreaterThan(2.0))
    f = 3.0 * x * x + 2.0 * x + 1.0
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    attr = MOI.ObjectiveFunction{typeof(f)}()
    MOI.set(model, attr, f)
    MOI.optimize!(model)
    @test isapprox(MOI.get(model, MOI.ObjectiveValue()), 17, atol = 1e-5)
    MOI.modify(model, attr, MOI.ScalarConstantChange(2.0))
    MOI.optimize!(model)
    @test isapprox(MOI.get(model, MOI.ObjectiveValue()), 18, atol = 1e-5)
    MOI.modify(model, attr, MOI.ScalarCoefficientChange(x, 3.0))
    MOI.optimize!(model)
    @test isapprox(MOI.get(model, MOI.ObjectiveValue()), 20, atol = 1e-5)
    MOI.modify(model, attr, MOI.ScalarQuadraticCoefficientChange(x, x, 8.0))
    MOI.optimize!(model)
    @test isapprox(MOI.get(model, MOI.ObjectiveValue()), 24, atol = 1e-5)
    return
end

function test_change_col_bounds_by_set_dimension_mismatch()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 3)
    c = MOI.add_constraint.(model, x, MOI.GreaterThan.(1.0:3.0))
    @test_throws(
        DimensionMismatch,
        MOI.set(model, MOI.ConstraintSet(), c, MOI.GreaterThan.([4.0, 5.0])),
    )
    return
end

function test_change_col_bounds_by_set_invalid()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    c = MOI.add_constraint(model, x, MOI.GreaterThan(0.0))
    c_invalid = typeof(c)(-123456)
    sets = MOI.GreaterThan.(1.0:2.0)
    @test_throws(
        MOI.InvalidIndex(c_invalid),
        MOI.set(model, MOI.ConstraintSet(), [c, c_invalid], sets),
    )
    return
end

function test_change_col_bounds_by_set_greater_than()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 3)
    c = MOI.add_constraint.(model, x, MOI.GreaterThan.(1.0:3.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    f = 1.0 * x[1] + x[2] + x[3]
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 6; atol = 1e-6)
    MOI.set(
        model,
        MOI.ConstraintSet(),
        [c[1], c[3]],
        MOI.GreaterThan.([4.0, 5.0]),
    )
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 11; atol = 1e-6)
    @test MOI.get(model, MOI.ConstraintSet(), c) ==
          MOI.GreaterThan.([4.0, 2.0, 5.0])
    return
end

function test_change_col_bounds_by_set_less_than()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 3)
    c = MOI.add_constraint.(model, x, MOI.LessThan.(1.0:3.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    f = 1.0 * x[1] + x[2] + x[3]
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 6; atol = 1e-6)
    MOI.set(model, MOI.ConstraintSet(), [c[1], c[3]], MOI.LessThan.([4.0, 5.0]))
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 11; atol = 1e-6)
    @test MOI.get(model, MOI.ConstraintSet(), c) ==
          MOI.LessThan.([4.0, 2.0, 5.0])
    return
end

function test_change_col_bounds_by_set_equal_to()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 3)
    c = MOI.add_constraint.(model, x, MOI.EqualTo.(1.0:3.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    f = 1.0 * x[1] + x[2] + x[3]
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 6; atol = 1e-6)
    MOI.set(model, MOI.ConstraintSet(), [c[1], c[3]], MOI.EqualTo.([4.0, 5.0]))
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 11; atol = 1e-6)
    @test MOI.get(model, MOI.ConstraintSet(), c) ==
          MOI.EqualTo.([4.0, 2.0, 5.0])
    return
end

function test_change_row_bounds_by_set_dimension_mismatch()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 3)
    c = MOI.add_constraint.(model, 1.0 .* x, MOI.GreaterThan.(1.0:3.0))
    @test_throws(
        DimensionMismatch,
        MOI.set(model, MOI.ConstraintSet(), c, MOI.GreaterThan.([4.0, 5.0])),
    )
    return
end

function test_change_row_bounds_by_set_invalid()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    c = MOI.add_constraint(model, 1.0 .* x, MOI.GreaterThan(0.0))
    c_invalid = typeof(c)(-123456)
    sets = MOI.GreaterThan.(1.0:2.0)
    @test_throws(
        MOI.InvalidIndex(c_invalid),
        MOI.set(model, MOI.ConstraintSet(), [c, c_invalid], sets),
    )
    return
end

function test_change_row_bounds_by_set_greater_than()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 3)
    c = MOI.add_constraint.(model, 1.0 .* x, MOI.GreaterThan.(1.0:3.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    f = 1.0 * x[1] + x[2] + x[3]
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 6; atol = 1e-6)
    MOI.set(
        model,
        MOI.ConstraintSet(),
        [c[1], c[3]],
        MOI.GreaterThan.([4.0, 5.0]),
    )
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 11; atol = 1e-6)
    @test MOI.get(model, MOI.ConstraintSet(), c) ==
          MOI.GreaterThan.([4.0, 2.0, 5.0])
    return
end

function test_change_row_bounds_by_set_less_than()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 3)
    c = MOI.add_constraint.(model, 1.0 .* x, MOI.LessThan.(1.0:3.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    f = 1.0 * x[1] + x[2] + x[3]
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 6; atol = 1e-6)
    MOI.set(model, MOI.ConstraintSet(), [c[1], c[3]], MOI.LessThan.([4.0, 5.0]))
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 11; atol = 1e-6)
    @test MOI.get(model, MOI.ConstraintSet(), c) ==
          MOI.LessThan.([4.0, 2.0, 5.0])
    return
end

function test_change_row_bounds_by_set_equal_to()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 3)
    c = MOI.add_constraint.(model, 1.0 .* x, MOI.EqualTo.(1.0:3.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    f = 1.0 * x[1] + x[2] + x[3]
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 6; atol = 1e-6)
    MOI.set(model, MOI.ConstraintSet(), [c[1], c[3]], MOI.EqualTo.([4.0, 5.0]))
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 11; atol = 1e-6)
    @test MOI.get(model, MOI.ConstraintSet(), c) ==
          MOI.EqualTo.([4.0, 2.0, 5.0])
    return
end

function test_set_affine_after_quadratic()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    MOI.add_constraint(model, x, MOI.GreaterThan(0.25))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    f = (x - 1.0)^2 + x + 1.0
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.VariablePrimal(), x), 0.5; atol = 1e-5)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 1.75; atol = 1e-5)
    g = 2.0 * x + 3.0
    MOI.set(model, MOI.ObjectiveFunction{typeof(g)}(), g)
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.VariablePrimal(), x), 0.25; atol = 1e-5)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 3.5; atol = 1e-5)
    return
end

function test_write_mps_gz()
    model = HiGHS.Optimizer()
    x = MOI.add_variable(model)
    filename = joinpath(mktempdir(), "test.mps.gz")
    @test HiGHS.Highs_writeModel(model, filename) == HiGHS.kHighsStatusWarning
    @test isfile(filename)
    return
end

function test_get_empty_objective_function()
    model = HiGHS.Optimizer()
    F = MOI.ScalarAffineFunction{Float64}
    f = MOI.get(model, MOI.ObjectiveFunction{F}())
    @test f ≈ zero(F)
    return
end

function test_is_valid_variable_bound()
    model = HiGHS.Optimizer()
    for S in (
        MOI.LessThan{Float64},
        MOI.GreaterThan{Float64},
        MOI.EqualTo{Float64},
        MOI.Interval{Float64},
        MOI.Semicontinuous{Float64},
        MOI.Semiinteger{Float64},
        MOI.ZeroOne,
        MOI.Integer,
    )
        ci = MOI.ConstraintIndex{MOI.VariableIndex,S}(1)
        @test_throws MOI.InvalidIndex HiGHS.column(model, ci)
        @test !MOI.is_valid(model, ci)
    end
    return
end

function test_set_function_variable()
    model = HiGHS.Optimizer()
    x, ci = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
    y = MOI.add_variable(model)
    @test_throws(
        MOI.SettingVariableIndexNotAllowed,
        MOI.set(model, MOI.ConstraintFunction(), ci, y),
    )
    return
end

function test_throw_if_existing_bound()
    model = HiGHS.Optimizer()
    x, ci = MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
    MOI.add_constraint(model, x, MOI.LessThan(1.0))
    @test_throws(
        MOI.LowerBoundAlreadySet,
        MOI.add_constraint(model, x, MOI.GreaterThan(1.0)),
    )
    @test_throws(
        MOI.UpperBoundAlreadySet,
        MOI.add_constraint(model, x, MOI.LessThan(0.5)),
    )
    return
end

function test_delete_double_bound_less_than()
    model = HiGHS.Optimizer()
    x = MOI.add_variable(model)
    _ = MOI.add_constraint(model, x, MOI.GreaterThan(0.0))
    ci = MOI.add_constraint(model, x, MOI.LessThan(1.0))
    @test MOI.is_valid(model, ci)
    MOI.delete(model, ci)
    @test !MOI.is_valid(model, ci)
    return
end

function test_delete_double_bound_greater_than()
    model = HiGHS.Optimizer()
    x = MOI.add_variable(model)
    ci = MOI.add_constraint(model, x, MOI.GreaterThan(0.0))
    _ = MOI.add_constraint(model, x, MOI.LessThan(1.0))
    @test MOI.is_valid(model, ci)
    MOI.delete(model, ci)
    @test !MOI.is_valid(model, ci)
    return
end

function test_set_constraint_function_constant_not_zero()
    model = HiGHS.Optimizer()
    x = MOI.add_variable(model)
    ci = MOI.add_constraint(model, 1.0 * x, MOI.GreaterThan(0.0))
    @test_throws(
        MOI.ScalarFunctionConstantNotZero,
        MOI.set(model, MOI.ConstraintFunction(), ci, 2.0 * x + 1.0),
    )
    return
end

function test_delete_integrality()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x, _ = MOI.add_constrained_variable(model, MOI.LessThan(1.5))
    y, _ = MOI.add_constrained_variable(model, MOI.LessThan(1.5))
    ci = MOI.add_constraint(model, y, MOI.Integer())
    f = 1.0 * x + 1.0 * y
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 2.5; atol = 1e-6)
    @test MOI.is_valid(model, ci)
    MOI.delete(model, ci)
    @test !MOI.is_valid(model, ci)
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 3.0; atol = 1e-6)
    return
end

function test_copy_to_unsupported_constraint()
    src = MOI.Utilities.Model{Float64}()
    x = MOI.add_variable(src)
    f = MOI.ScalarNonlinearFunction(:log, Any[x])
    MOI.add_constraint(src, f, MOI.EqualTo(0.0))
    model = HiGHS.Optimizer()
    @test_throws(
        MOI.UnsupportedConstraint{typeof(f),MOI.EqualTo{Float64}},
        MOI.copy_to(model, src),
    )
    return
end

function test_delete_constraint_twice()
    model = HiGHS.Optimizer()
    x = MOI.add_variable(model)
    ci = MOI.add_constraint(model, 1.0 * x, MOI.EqualTo(1.0))
    @test_throws(
        ErrorException(
            "Encountered an error in HiGHS (Status -1). Check the log " *
            "for details.",
        ),
        MOI.delete(model, [ci, ci]),
    )
    return
end

function test_RawStatusStringOptimizeNotCalled()
    model = HiGHS.Optimizer()
    @test MOI.get(model, MOI.RawStatusString()) == "OPTIMIZE_NOT_CALLED"
    @test MOI.get(model, MOI.ResultCount()) == 0
    return
end

function test_nonbasic_equality_constraint()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    ci = MOI.add_constraint(model, 1.0 * x, MOI.EqualTo(1.0))
    MOI.optimize!(model)
    @test MOI.get(model, MOI.ConstraintBasisStatus(), ci) == MOI.NONBASIC
    return
end

function test_variable_basis_status_zero()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariableBasisStatus(), x) == MOI.SUPER_BASIC
    return
end

function test_optimize_errored()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variable(model)
    MOI.set(model, MOI.NumberOfThreads(), 123_456)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.OTHER_ERROR
    @test MOI.get(model, MOI.RawStatusString()) ==
          "There was an error calling optimize!"
    return
end

function test_infeasible_point()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 2)
    MOI.add_constraint.(model, x, MOI.GreaterThan(0.0))
    ci = MOI.add_constraint(model, x[1] + 2.0 * x[2], MOI.LessThan(-1.0))
    MOI.set(model, MOI.RawOptimizerAttribute("presolve"), "off")
    MOI.set(model, MOI.RawOptimizerAttribute("solver"), "ipm")
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.INFEASIBLE
    @test MOI.get(model, MOI.PrimalStatus()) == MOI.INFEASIBLE_POINT
    dual_stat = MOI.get(model, MOI.DualStatus())
    @test dual_stat in (MOI.INFEASIBLE_POINT, MOI.NO_SOLUTION)
    @test MOI.get(model, MOI.ResultCount()) == 1
    @test MOI.get(model, MOI.VariablePrimal(), x) isa Vector{Float64}
    @test MOI.get(model, MOI.ConstraintDual(), ci) isa Float64
    return
end

function test_DualObjectiveValue_int()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x, _ = MOI.add_constrained_variable(model, MOI.ZeroOne())
    MOI.optimize!(model)
    @test isnan(MOI.get(model, MOI.DualObjectiveValue()))
    return
end

function test_continuous_objective_bound()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, MOI.RawOptimizerAttribute("solver"), "ipm")
    MOI.set(model, MOI.RawOptimizerAttribute("run_crossover"), "off")
    MOI.set(model, MOI.RawOptimizerAttribute("presolve"), "off")
    x = MOI.add_variables(model, 3)
    c = MOI.add_constraint.(model, 1.0 .* x, MOI.EqualTo.(1.0:3.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    f = 1.0 * x[1] + x[2] + x[3]
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    primal = MOI.get(model, MOI.ObjectiveValue())
    dual = MOI.get(model, MOI.DualObjectiveValue())
    @test MOI.get(model, MOI.ObjectiveBound()) == dual
    @test 0 < MOI.get(model, MOI.RelativeGap()) <= 1e-6
    return
end

function test_callback_interrupt()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.RawOptimizerAttribute("solver"), "ipm")
    MOI.set(model, MOI.RawOptimizerAttribute("presolve"), "off")
    x = MOI.add_variables(model, 3)
    c = MOI.add_constraint.(model, 1.0 .* x, MOI.EqualTo.(1.0:3.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    f = 1.0 * x[1] + x[2] + x[3]
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    callback_types = Cint[]
    function user_callback(
        callback_type::Cint,
        ::Ptr{Cchar},
        ::HiGHS.HighsCallbackDataOut,
    )::Cint
        push!(callback_types, callback_type)
        return 1
    end
    MOI.set(model, HiGHS.CallbackFunction(), user_callback)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.INTERRUPTED
    @test HiGHS.kHighsCallbackIpmInterrupt in callback_types
    return
end

function test_active_bound()
    for ((l, x, u, d), result) in [
        # Primal exists. Pick closest bound lower
        (0.0, 0.0, 1.0, 2.0) => 0.0,
        (0.0, 0.4, 1.0, 2.0) => 0.0,
        (0.0, 0.4, 1.0, -2.0) => 0.0,  # incorrect d but doesn't matter
        # Primal exists. Pick closest bound upper
        (0.0, 1.0, 1.0, -2.0) => 1.0,
        (0.0, 0.6, 1.0, -2.0) => 1.0,
        (0.0, 0.6, 1.0, 2.0) => 1.0,  # incorrect d but doesn't matter
        (-Inf, 0.0, Inf, 0.0) => 0.0,
        # It's a ray. Choose based on sign
        (0.0, NaN, 1.0, 2.0) => 0.0,
        (0.0, NaN, 1.0, 1e-10) => 0.0,
        (0.0, NaN, 1.0, -2.0) => 1.0,
        (0.0, NaN, 1.0, -1e-10) => 1.0,
        (-Inf, NaN, Inf, 0.0) => 0.0,
        # It's a one-sided ray
        (0.0, NaN, Inf, 2.0) => 0.0,
        (0.0, NaN, Inf, -1e-10) => 0.0,
        (-Inf, NaN, 1.0, 1e-10) => 1.0,
        (-Inf, NaN, 1.0, -2.0) => 1.0,
    ]
        @test HiGHS._active_bound(l, x, u, d) == result
    end
    return
end

function test_add_constrained_variable_tuple()
    F = MOI.VariableIndex
    model = HiGHS.Optimizer()
    set = (MOI.GreaterThan(0.0), MOI.LessThan(1.0))
    @test MOI.supports_add_constrained_variable(model, typeof(set))
    x, (c_l, c_u) = MOI.add_constrained_variable(model, set)
    @test c_l == MOI.ConstraintIndex{F,MOI.GreaterThan{Float64}}(x.value)
    @test c_u == MOI.ConstraintIndex{F,MOI.LessThan{Float64}}(x.value)
    @test MOI.get(model, MOI.ConstraintFunction(), c_l) == x
    @test MOI.get(model, MOI.ConstraintSet(), c_l) == set[1]
    @test MOI.get(model, MOI.ConstraintFunction(), c_u) == x
    @test MOI.get(model, MOI.ConstraintSet(), c_u) == set[2]
    return
end

function test_dual_objective_value_infeasible()
    model = HiGHS.Optimizer()
    x = MOI.add_variable(model)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    f = 1.0 * x
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBLE_POINT
    @test MOI.get(model, MOI.DualObjectiveValue()) == 0.0
    return
end

function test_ObjectiveFunction_VectorAffineFunction_to_ScalarAffineFunction()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 2)
    MOI.add_constraint(model, x[1], MOI.Interval(1.0, 2.0))
    MOI.add_constraint(model, x[2], MOI.Interval(3.0, 4.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    f = MOI.Utilities.operate(vcat, Float64, 1.0 * x[1], 1.0 * x[2])
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 4.0)
    g = 1.0 * x[1] - 2.0 * x[2]
    MOI.set(model, MOI.ObjectiveFunction{typeof(g)}(), g)
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), -7.0; atol = 1e-6)
    @test ≈(MOI.get(model, MOI.VariablePrimal(), x), [1.0, 4.0]; atol = 1e-6)
    @test MOI.get(model, MOI.ObjectiveFunctionType()) ==
          MOI.ScalarAffineFunction{Float64}
    @test model.multi_objective === nothing
    return
end

function test_ObjectiveFunction_VectorAffineFunction_to_ScalarQuadraticFunction()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 2)
    MOI.add_constraint(model, x[1], MOI.Interval(1.0, 2.0))
    MOI.add_constraint(model, x[2], MOI.Interval(3.0, 4.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    f = MOI.Utilities.operate(vcat, Float64, 1.0 * x[1], 1.0 * x[2])
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 4.0)
    g = 1.0 * x[1] + (1.0 * x[2] * x[2] - 8.0 * x[2] + 16.0)
    MOI.set(model, MOI.ObjectiveFunction{typeof(g)}(), g)
    MOI.optimize!(model)
    @test ≈(MOI.get(model, MOI.ObjectiveValue()), 1.0; atol = 1e-6)
    @test ≈(MOI.get(model, MOI.VariablePrimal(), x), [1.0, 4.0]; atol = 1e-6)
    @test MOI.get(model, MOI.ObjectiveFunctionType()) ==
          MOI.ScalarQuadraticFunction{Float64}
    @test model.multi_objective === nothing
    return
end

function test_write_to_file()
    dir = mktempdir()
    model = HiGHS.Optimizer()
    x, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(1.0))
    MOI.set(model, MOI.VariableName(), x, "x")
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    f = 1.0 * x + 2.0
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    filename = joinpath(dir, "model.mps")
    MOI.write_to_file(model, filename)
    contents = read(filename, String)
    @test occursin("ENDATA", contents)
    filename = joinpath(dir, "model.lp")
    MOI.write_to_file(model, filename)
    contents = read(filename, String)
    @test occursin("obj:", contents)
    @test occursin("bounds", contents)
    return
end

function test_row_type()
    @test HiGHS._row_type(MOI.GreaterThan(0.0)) == HiGHS._ROW_TYPE_GREATERTHAN
    @test HiGHS._row_type(MOI.LessThan(0.0)) == HiGHS._ROW_TYPE_LESSTHAN
    @test HiGHS._row_type(MOI.EqualTo(0.0)) == HiGHS._ROW_TYPE_EQUAL_TO
    @test HiGHS._row_type(MOI.Interval(0.0, 1.0)) == HiGHS._ROW_TYPE_INTERVAL
    return
end

end  # module

TestMOIHighs.runtests()
