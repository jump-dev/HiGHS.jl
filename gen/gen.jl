using Clang.Generators
import HiGHS_jll

highs = joinpath(HiGHS_jll.artifact_dir, "include")
c_api = joinpath(highs, "interfaces", "highs_c_api.h")

build!(
    create_context(
        [c_api, joinpath(highs, "util", "HighsInt.h")],
        vcat(get_default_args(), "-I$highs"),
        load_options(joinpath(@__DIR__, "generate.toml")),
    ),
)

open(joinpath(@__DIR__, "..", "src", "gen", "libhighs.jl"), "a") do io
    for line in readlines(c_api)
        m = match(r"const HighsInt kHighs([a-zA-Z]+) = (-?[0-9]+);", line)
        if m === nothing
            continue
        end
        println(io, "const kHighs$(m[1]) = HighsInt($(m[2]))")
    end
end
