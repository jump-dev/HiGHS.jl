# Copyright (c) 2019 Mathieu Besançon, Oscar Dowson, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using Test

@testset "$(file)" for file in readdir(@__DIR__)
    if file == "runtests.jl" || !endswith(file, ".jl")
        continue
    end
    include(file)
end
