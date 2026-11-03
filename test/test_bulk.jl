# Copyright (c) 2019 Mathieu Besançon, Oscar Dowson, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestBulkHiGHS

using Test
import HiGHS
import MathOptInterface as MOI

const F = MOI.ScalarAffineFunction{Float64}
const SETS = (
    MOI.GreaterThan{Float64},
    MOI.LessThan{Float64},
    MOI.EqualTo{Float64},
    MOI.Interval{Float64},
)

function optimizer()
    model = HiGHS.Optimizer()
    MOI.set(model, MOI.Silent(), true)
    return model
end

function native_columns(model)
    n = HiGHS.Highs_getNumCol(model)
    cost, lower, upper = zeros(n), zeros(n), zeros(n)
    num_col, nnz = Ref{HiGHS.HighsInt}(0), Ref{HiGHS.HighsInt}(0)
    ret = HiGHS.Highs_getColsByRange(
        model,
        0,
        n - 1,
        num_col,
        cost,
        lower,
        upper,
        nnz,
        C_NULL,
        C_NULL,
        C_NULL,
    )
    @test ret == HiGHS.kHighsStatusOk
    @test num_col[] == n
    return cost, lower, upper
end

set(::Type{MOI.GreaterThan{Float64}}, x) = MOI.GreaterThan(Float64(x))
set(::Type{MOI.LessThan{Float64}}, x) = MOI.LessThan(Float64(x))
set(::Type{MOI.EqualTo{Float64}}, x) = MOI.EqualTo(Float64(x))
function set(::Type{MOI.Interval{Float64}}, x)
    return MOI.Interval(Float64(x), Float64(x + 0.5))
end

function test_add_variables()
    model = optimizer()
    first = MOI.add_variable(model)
    MOI.set(model, MOI.VariableName(), first, "first")
    @test MOI.get(model, MOI.VariableIndex, "first") == first
    variables = MOI.add_variables(model, 3)
    @test variables == MOI.VariableIndex.([2, 3, 4])
    @test all(MOI.is_valid.(model, variables))
    @test HiGHS.column.(model, variables) == [1, 2, 3]
    MOI.delete(model, variables[1])
    MOI.delete(model, variables[3])
    added = MOI.add_variables(model, 2)
    @test added == MOI.VariableIndex.([5, 6])
    @test HiGHS.column.(model, added) == [2, 3]
    scalar = MOI.add_variable(model)
    @test scalar == MOI.VariableIndex(7)
    @test HiGHS.column(model, scalar) == 4
    @test MOI.get(model, MOI.ListOfVariableIndices()) ==
          [first, variables[2], added..., scalar]
    @test MOI.get(model, MOI.NumberOfVariables()) == 5
    @test HiGHS.Highs_getNumCol(model) == 5
    @test MOI.get(model, MOI.VariableIndex, "first") == first
    MOI.set(model, MOI.VariableName(), added[1], "added")
    @test MOI.get(model, MOI.VariableIndex, "added") == added[1]
    @test MOI.add_variables(model, 0) == MOI.VariableIndex[]
    @test eltype(MOI.add_variables(model, 0)) == MOI.VariableIndex
    @test_throws ArgumentError MOI.add_variables(model, -1)
    @test MOI.get(model, MOI.NumberOfVariables()) == 5
    cost, lower, upper = native_columns(model)
    @test cost == zeros(5)
    @test lower == fill(-Inf, 5)
    @test upper == fill(Inf, 5)
    MOI.empty!(model)
    @test MOI.add_variables(model, 2) == MOI.VariableIndex.([1, 2])
    return
end

function test_add_bounds()
    for S in SETS
        model = optimizer()
        deleted = MOI.add_variable(model)
        MOI.delete(model, deleted)
        variables = MOI.add_variables(model, 3)
        opposite = if S == MOI.GreaterThan{Float64}
            MOI.add_constraint.(model, variables, MOI.LessThan(8.0))
        elseif S == MOI.LessThan{Float64}
            MOI.add_constraint.(model, variables, MOI.GreaterThan(0.0))
        else
            nothing
        end
        order = [3, 1, 2]
        sets = [set(S, i) for i in 1:3]
        indices = MOI.add_constraints(model, variables[order], sets[order])
        @test indices == MOI.ConstraintIndex{MOI.VariableIndex,S}.(
            getfield.(variables[order], :value),
        )
        @test MOI.get(model, MOI.ConstraintSet(), indices) == sets[order]
        @test MOI.get(model, MOI.ConstraintFunction(), indices) ==
              variables[order]
        @test all(MOI.is_valid.(model, indices))
        @test MOI.get(model, MOI.NumberOfConstraints{MOI.VariableIndex,S}()) ==
              3
        @test MOI.get(
            model,
            MOI.ListOfConstraintIndices{MOI.VariableIndex,S}(),
        ) == indices[[2, 3, 1]]
        if opposite !== nothing
            @test all(MOI.is_valid.(model, opposite))
        end
        _, lower, upper = native_columns(model)
        @test lower == (S == MOI.LessThan{Float64} ? zeros(3) : [1.0, 2.0, 3.0])
        @test upper == (
            S == MOI.GreaterThan{Float64} ? fill(8.0, 3) :
            S == MOI.Interval{Float64} ? [1.5, 2.5, 3.5] : [1.0, 2.0, 3.0]
        )
        sense = S == MOI.LessThan{Float64} ? MOI.MAX_SENSE : MOI.MIN_SENSE
        MOI.set(model, MOI.ObjectiveSense(), sense)
        objective = sum(1.0 * v for v in variables)
        MOI.set(model, MOI.ObjectiveFunction{F}(), objective)
        MOI.optimize!(model)
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
        @test MOI.get(model, MOI.ObjectiveValue()) ≈ 6.0
        @test MOI.get(model, MOI.VariablePrimal(), variables) ≈ [1.0, 2.0, 3.0]
        expected_dual = S == MOI.LessThan{Float64} ? -ones(3) : ones(3)
        @test MOI.get(model, MOI.ConstraintDual(), indices) ≈ expected_dual
        MOI.delete(model, indices[1])
        @test !MOI.is_valid(model, indices[1])
        scalar = MOI.add_constraint(model, variables[3], sets[3])
        @test MOI.is_valid(model, scalar)
        @test MOI.get(model, MOI.ConstraintSet(), scalar) == sets[3]
        @test isempty(MOI.add_constraints(model, MOI.VariableIndex[], S[]))
    end
    return
end

function test_bound_column_ranges()
    for S in SETS,
        deleted in (false, true),
        order in ([2, 3, 4], [2, 4, 6], [4, 2, 3], [3])

        model = optimizer()
        variables = MOI.add_variables(model, 7)
        if deleted
            MOI.delete(model, variables[3])
            deleteat!(variables, 3)
        end
        n = length(variables)
        lower, upper = fill(-Inf, n), fill(Inf, n)
        if S == MOI.GreaterThan{Float64}
            MOI.add_constraints(model, variables, fill(MOI.LessThan(100.0), n))
            fill!(upper, 100.0)
        elseif S == MOI.LessThan{Float64}
            MOI.add_constraints(
                model,
                variables,
                fill(MOI.GreaterThan(-100.0), n),
            )
            fill!(lower, -100.0)
        end
        sets = [set(S, i) for i in eachindex(order)]
        indices = MOI.add_constraints(model, variables[order], sets)
        @test all(MOI.is_valid.(model, indices))
        @test MOI.get(model, MOI.ConstraintSet(), indices) == sets
        for (column, bound) in zip(order, sets)
            if S != MOI.LessThan{Float64}
                lower[column] =
                    S == MOI.EqualTo{Float64} ? bound.value : bound.lower
            end
            if S != MOI.GreaterThan{Float64}
                upper[column] =
                    S == MOI.EqualTo{Float64} ? bound.value : bound.upper
            end
        end
        cost, actual_lower, actual_upper = native_columns(model)
        @test cost == zeros(n)
        @test actual_lower == lower
        @test actual_upper == upper
    end
    return
end

function test_bounds_reject_without_partial_updates()
    for S in SETS
        model = optimizer()
        variables = MOI.add_variables(model, 3)
        before = native_columns(model)
        invalid = MOI.VariableIndex(-1)
        @test_throws MOI.InvalidIndex(invalid) MOI.add_constraints(
            model,
            [variables[1], invalid],
            fill(set(S, 1), 2),
        )
        @test_throws DimensionMismatch MOI.add_constraints(
            model,
            variables,
            [set(S, 1)],
        )
        for duplicate in (
            [variables[1], variables[1]],
            [variables[3], variables[1], variables[3]],
        )
            error =
                S == MOI.LessThan{Float64} ?
                MOI.UpperBoundAlreadySet{S,S}(duplicate[1]) :
                MOI.LowerBoundAlreadySet{S,S}(duplicate[1])
            @test_throws error MOI.add_constraints(
                model,
                duplicate,
                fill(set(S, 1), length(duplicate)),
            )
        end
        @test MOI.get(model, MOI.NumberOfConstraints{MOI.VariableIndex,S}()) ==
              0
        @test native_columns(model) == before
        MOI.add_constraint(model, variables[2], set(S, 1))
        before = native_columns(model)
        error =
            S == MOI.LessThan{Float64} ?
            MOI.UpperBoundAlreadySet{S,S}(variables[2]) :
            MOI.LowerBoundAlreadySet{S,S}(variables[2])
        @test_throws error MOI.add_constraints(
            model,
            variables[1:2],
            fill(set(S, 2), 2),
        )
        @test native_columns(model) == before
        @test MOI.get(model, MOI.NumberOfConstraints{MOI.VariableIndex,S}()) ==
              1
    end
    for (old, new, E) in (
        (MOI.Interval(0.0, 3.0), MOI.EqualTo(1.0), MOI.LowerBoundAlreadySet),
        (MOI.EqualTo(1.0), MOI.Interval(0.0, 3.0), MOI.LowerBoundAlreadySet),
        (MOI.LessThan(3.0), MOI.EqualTo(1.0), MOI.UpperBoundAlreadySet),
        (
            MOI.Semicontinuous(1.0, 3.0),
            MOI.GreaterThan(0.0),
            MOI.LowerBoundAlreadySet,
        ),
        (
            MOI.Semiinteger(1.0, 3.0),
            MOI.LessThan(4.0),
            MOI.UpperBoundAlreadySet,
        ),
    )
        model = optimizer()
        x = MOI.add_variables(model, 2)
        MOI.add_constraint(model, x[2], old)
        before = native_columns(model)
        @test_throws E MOI.add_constraints(model, x, fill(new, 2))
        @test native_columns(model) == before
    end
    model = optimizer()
    x = MOI.add_variables(model, 2)
    before = native_columns(model)
    @test_throws ErrorException MOI.add_constraints(
        model,
        x,
        MOI.GreaterThan.([1.0, NaN]),
    )
    @test native_columns(model) == before
    @test isempty(MOI.get(model, MOI.ListOfConstraintTypesPresent()))
    return
end

function test_binary_bounds_and_mixed_sets()
    model = optimizer()
    x = MOI.add_variables(model, 2)
    binary = MOI.add_constraint.(model, x, MOI.ZeroOne())
    lower = MOI.add_constraints(model, x, MOI.GreaterThan.([-1.0, 0.0]))
    upper = MOI.add_constraints(model, reverse(x), MOI.LessThan.([0.5, 2.0]))
    @test all(MOI.is_valid.(model, binary))
    @test MOI.get(model, MOI.ConstraintSet(), lower) ==
          MOI.GreaterThan.([-1.0, 0.0])
    @test MOI.get(model, MOI.ConstraintSet(), upper) ==
          MOI.LessThan.([0.5, 2.0])
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{F}(), 1.0 * x[1] + x[2])
    MOI.optimize!(model)
    @test MOI.get(model, MOI.VariablePrimal(), x) ≈ [1.0, 0.0]
    model = optimizer()
    x = MOI.add_variable(model)
    sets = Union{MOI.GreaterThan{Float64},MOI.LessThan{Float64}}[
        MOI.GreaterThan(1.0),
        MOI.LessThan(2.0),
    ]
    indices = MOI.add_constraints(model, [x, x], sets)
    @test typeof(indices[1]) ==
          MOI.ConstraintIndex{MOI.VariableIndex,typeof(sets[1])}
    @test typeof(indices[2]) ==
          MOI.ConstraintIndex{MOI.VariableIndex,typeof(sets[2])}
    @test all(MOI.is_valid.(model, indices))
    @test MOI.get(model, MOI.ConstraintSet(), indices) == sets
    return
end

function test_affine_rows()
    for S in SETS
        model = optimizer()
        original = MOI.add_variables(model, 4)
        MOI.delete(model, original[2])
        x = original[[1, 3, 4]]
        MOI.add_constraints(model, x, fill(MOI.Interval(0.0, 2.0), 3))
        deleted = MOI.add_constraint(model, 1.0 * x[1], MOI.LessThan(20.0))
        existing = MOI.add_constraint(model, 1.0 * x[3], MOI.LessThan(20.0))
        MOI.delete(model, deleted)
        MOI.set(model, MOI.ConstraintName(), existing, "existing")
        @test MOI.get(model, MOI.ConstraintIndex, "existing") == existing
        functions = [
            F(
                [
                    MOI.ScalarAffineTerm(2.0, x[3]),
                    MOI.ScalarAffineTerm(3.0, x[1]),
                    MOI.ScalarAffineTerm(-1.0, x[3]),
                ],
                0.0,
            ),
            F(
                [
                    MOI.ScalarAffineTerm(1.0, x[2]),
                    MOI.ScalarAffineTerm(-1.0, x[2]),
                ],
                0.0,
            ),
            1.0 * x[2],
            F(MOI.ScalarAffineTerm{Float64}[], 0.0),
        ]
        before = deepcopy(functions)
        sets =
            S == MOI.Interval{Float64} ?
            [
                MOI.Interval(3.0, 4.0),
                MOI.Interval(0.0, 0.0),
                MOI.Interval(1.0, 2.0),
                MOI.Interval(0.0, 0.0),
            ] : [set(S, 3.0), set(S, 0.0), set(S, 1.0), set(S, 0.0)]
        indices = MOI.add_constraints(model, functions, sets)
        @test getfield.(indices, :value) == [3, 4, 5, 6]
        @test HiGHS.row.(model, indices) == [1, 2, 3, 4]
        @test HiGHS.Highs_getNumRow(model) == 5
        @test all(MOI.is_valid.(model, indices))
        @test MOI.get(model, MOI.ConstraintSet(), indices) == sets
        for i in eachindex(functions)
            @test MOI.get(model, MOI.ConstraintFunction(), indices[i]) ≈
                  before[i]
            @test functions[i].terms == before[i].terms
        end
        @test MOI.get(model, MOI.ConstraintIndex, "existing") == existing
        MOI.set(model, MOI.ConstraintName(), indices[1], "bulk")
        @test MOI.get(model, MOI.ConstraintIndex, "bulk") == indices[1]
        empty!(functions[1].terms)
        @test MOI.get(model, MOI.ConstraintFunction(), indices[1]) ≈ before[1]
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(model, MOI.ObjectiveFunction{F}(), sum(1.0 * v for v in x))
        MOI.optimize!(model)
        @test MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
        @test MOI.get(model, MOI.ObjectiveValue()) ≈
              (S == MOI.LessThan{Float64} ? 0.0 : 2.0)
        @test MOI.get(model, MOI.ConstraintPrimal(), indices[2]) ≈ 0.0
        MOI.delete(model, indices[2])
        @test !MOI.is_valid(model, indices[2])
        @test HiGHS.row(model, indices[4]) == 3
        added = MOI.add_constraint(model, 1.0 * x[1], MOI.LessThan(2.0))
        @test HiGHS.row(model, added) == 4
        @test isempty(MOI.add_constraints(model, F[], S[]))
    end
    return
end

function test_affine_rows_reject_without_partial_updates()
    model = optimizer()
    x = MOI.add_variables(model, 2)
    f = 1.0 * x[1] + x[2]
    existing = MOI.add_constraint(model, f, MOI.LessThan(3.0))
    @test_throws DimensionMismatch MOI.add_constraints(
        model,
        [f],
        MOI.EqualTo{Float64}[],
    )
    @test_throws MOI.ScalarFunctionConstantNotZero MOI.add_constraints(
        model,
        [f, F(copy(f.terms), 1.0)],
        fill(MOI.EqualTo(0.0), 2),
    )
    invalid = MOI.VariableIndex(-1)
    @test_throws MOI.InvalidIndex(invalid) MOI.add_constraints(
        model,
        [f, 1.0 * invalid],
        fill(MOI.LessThan(0.0), 2),
    )
    @test_throws ErrorException MOI.add_constraints(
        model,
        [f, 1e30 * x[1]],
        fill(MOI.LessThan(0.0), 2),
    )
    @test HiGHS.Highs_getNumRow(model) == 1
    @test length(model.affine_constraint_info) == 1
    @test MOI.get(model, MOI.ConstraintFunction(), existing) ≈ f
    added = MOI.add_constraints(model, [f], [MOI.EqualTo(2.0)])
    @test only(added).value == 2
    @test HiGHS.row(model, only(added)) == 1
    sets = Union{MOI.LessThan{Float64},MOI.GreaterThan{Float64}}[
        MOI.LessThan(3.0),
        MOI.GreaterThan(1.0),
    ]
    mixed = MOI.add_constraints(model, [f, f], sets)
    @test all(MOI.is_valid.(model, mixed))
    @test MOI.get(model, MOI.ConstraintSet(), mixed) == sets
    @test typeof(mixed[1]) == MOI.ConstraintIndex{F,typeof(sets[1])}
    @test typeof(mixed[2]) == MOI.ConstraintIndex{F,typeof(sets[2])}
    return
end

function runtests()
    is_test(name) = startswith("$name", "test_")
    @testset "$name" for name in filter(is_test, names(@__MODULE__; all = true))
        getfield(@__MODULE__, name)()
    end
    return
end

end # module TestBulkHiGHS

TestBulkHiGHS.runtests()
