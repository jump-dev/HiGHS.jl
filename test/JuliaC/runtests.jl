# Copyright (c) 2019 Mathieu Besançon, Oscar Dowson, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using Test

import JuliaC

function compile(output_dir)
    image_recipe = JuliaC.ImageRecipe(
        output_type = "--output-exe",
        file = joinpath(@__DIR__, "MyApp"),
        trim_mode = "no",
        add_ccallables = false,
        verbose = true,
    )
    link_recipe = JuliaC.LinkRecipe(;
        image_recipe,
        outname = joinpath(output_dir, "MyApp"),
    )
    bundle_recipe = JuliaC.BundleRecipe(; link_recipe, output_dir)
    JuliaC.compile_products(image_recipe)
    JuliaC.link_products(link_recipe)
    JuliaC.bundle_products(bundle_recipe)
    return
end

@testset "JuliaC" begin
    output_dir = mktempdir()
    compile(output_dir)
    app = joinpath(output_dir, "bin", "MyApp")
    output = sprint(io -> run(pipeline(`$app`; stdout = io)))
    @test occursin("HiGHS", output)
    @test occursin("Optimal", output)
end
