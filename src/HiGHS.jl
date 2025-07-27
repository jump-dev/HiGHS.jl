# Copyright (c) 2019 Mathieu Besançon, Oscar Dowson, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module HiGHS

import HiGHS_jll: libhighs
import MathOptInterface as MOI
import MathOptIIS
import SparseArrays

include("gen/libhighs.jl")
include("MOI_wrapper.jl")

# HiGHS exports all `Highs_xxx` symbols. If you don't want all of these symbols
# in your environment, then use `import HiGHS` instead of `using HiGHS`.

for sym in names(@__MODULE__, all = true)
    if startswith("$sym", "Highs_") || startswith("$sym", "kHighs")
        @eval export $sym
    end
end

import PrecompileTools

PrecompileTools.@setup_workload begin
    PrecompileTools.@compile_workload begin
        let
            model = MOI.Utilities.CachingOptimizer(
                MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
                MOI.instantiate(HiGHS.Optimizer; with_bridge_type = Float64),
            )
            MOI.set(model, MOI.Silent(), true)
            x = MOI.add_variables(model, 3)
            MOI.supports(model, MOI.VariableName(), typeof(x[1]))
            MOI.set(model, MOI.VariableName(), x[1], "x1")
            MOI.set(model, MOI.VariablePrimalStart(), x[1], 0.0)
            MOI.add_constraint(model, x[1], MOI.ZeroOne())
            MOI.add_constraint(model, x[2], MOI.Integer())
            for F in (MOI.VariableIndex, MOI.ScalarAffineFunction{Float64})
                MOI.supports_constraint(model, F, MOI.GreaterThan{Float64})
                MOI.supports_constraint(model, F, MOI.LessThan{Float64})
                MOI.supports_constraint(model, F, MOI.EqualTo{Float64})
            end
            MOI.supports_constraint(model, MOI.VariableIndex, MOI.ZeroOne)
            MOI.supports_constraint(model, MOI.VariableIndex, MOI.Integer)
            MOI.add_constraint(model, x[1], MOI.GreaterThan(0.0))
            MOI.add_constraint(model, x[2], MOI.LessThan(0.0))
            MOI.add_constraint(model, x[3], MOI.EqualTo(0.0))
            MOI.add_constrained_variable(model, MOI.GreaterThan(0.0))
            MOI.add_constrained_variable(model, MOI.LessThan(0.0))
            MOI.add_constrained_variable(model, MOI.EqualTo(0.0))
            MOI.add_constrained_variable(model, MOI.Integer())
            MOI.add_constrained_variable(model, MOI.ZeroOne())
            set = (MOI.GreaterThan(0.0), MOI.LessThan(0.0))
            MOI.supports_add_constrained_variable(model, typeof(set))
            MOI.add_constrained_variable(model, set)
            f = 1.0 * x[1] + x[2] + x[3]
            c1 = MOI.add_constraint(model, f, MOI.GreaterThan(0.0))
            MOI.set(model, MOI.ConstraintName(), c1, "c1")
            MOI.supports(model, MOI.ConstraintName(), typeof(c1))
            MOI.add_constraint(model, f, MOI.LessThan(0.0))
            MOI.add_constraint(model, f, MOI.EqualTo(0.0))
            y, _ = MOI.add_constrained_variables(model, MOI.Nonnegatives(2))
            MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
            MOI.supports(model, MOI.ObjectiveFunction{typeof(f)}())
            MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
            MOI.optimize!(model)
            MOI.get(model, MOI.TerminationStatus())
            MOI.get(model, MOI.PrimalStatus())
            MOI.get(model, MOI.DualStatus())
            MOI.get(model, MOI.VariablePrimal(), x)
        end
    end
end

end # module HiGHS
