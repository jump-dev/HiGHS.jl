# Copyright (c) 2019 Mathieu Besançon, Oscar Dowson, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestRunTestsCache

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

end  # module

TestRunTestsCache.runtests()
