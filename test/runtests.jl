# Copyright (c) 2019 Mathieu Besançon, Oscar Dowson, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

import HiGHS
import ParallelTestRunner
import Test

push!(ARGS, "--jobs=4")

is_test_file(f) = startswith(f, "test_") && endswith(f, ".jl")

testsuite = Dict{String,Expr}()
for (root, dirs, files) in walkdir(@__DIR__)
    for file in joinpath.(root, filter(is_test_file, files))
        testsuite[file] = :(include($file))
    end
end

ParallelTestRunner.runtests(HiGHS, ARGS; testsuite)
