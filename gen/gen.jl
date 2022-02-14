using Clang.Generators
import HiGHS_jll

highs = joinpath(HiGHS_jll.artifact_dir, "include")
build!(
    create_context(
        [
            joinpath(highs, "interfaces", "highs_c_api.h"),
            joinpath(highs, "util", "HighsInt.h"),
        ],
        vcat(get_default_args(), "-I$highs"),
        load_options(joinpath(@__DIR__, "generate.toml")),
    ),
)
