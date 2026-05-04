# Copyright (c) 2019 Mathieu Besançon, Oscar Dowson, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

import HiGHS
import ParallelTestRunner
import Test

is_test_file(f) = startswith(f, "test_") && endswith(f, ".jl")

testsuite = Dict{String,Expr}()
for (root, dirs, files) in walkdir(@__DIR__)
    for file in joinpath.(root, filter(is_test_file, files))
        testsuite[file] = :(include($file))
    end
end

if "--parallel=false" in ARGS
    Test.@testset "$file" for file in keys(testsuite)
        include(file)
    end
else
    ParallelTestRunner.runtests(HiGHS, ARGS; testsuite)
end
