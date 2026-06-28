# Copyright (c) 2019 Mathieu Besançon, Oscar Dowson, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestRunTests

using Test

import HiGHS
import MathOptInterface as MOI

function runtests()
    is_test(name) = startswith("$name", "test_")
    @testset "$name" for name in filter(is_test, names(@__MODULE__; all = true))
        getfield(@__MODULE__, name)()
    end
    return
end

function test_runtests()
    model = MOI.Bridges.full_bridge_optimizer(HiGHS.Optimizer(), Float64)
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, MOI.RawOptimizerAttribute("parallel"), "on")
    # Slightly loosen tolerances, particularly for QP tests
    MOI.Test.runtests(model, MOI.Test.Config(; atol = 1e-7))
    return
end

end  # TestRunTests

TestRunTests.runtests()
