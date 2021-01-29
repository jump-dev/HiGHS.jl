import Clang
import HiGHS_jll

const LIBHIGHS_HEADERS = [
    joinpath(HiGHS_jll.artifact_dir, "include", "interfaces", "highs_c_api.h"),
]

const GEN_DIR = joinpath(@__DIR__, "..", "src", "wrapper")
if !isdir(GEN_DIR)
    mkdir(GEN_DIR)
end

wc = Clang.init(
    headers = LIBHIGHS_HEADERS,
    output_file = joinpath(GEN_DIR, "libhighs_api.jl"),
    common_file = joinpath(GEN_DIR, "libhighs_common.jl"),
    header_wrapped = (root, current) -> root == current,
    header_library = x -> "libhighs",
    clang_diagnostics = true,
)

run(wc)

rm(joinpath(GEN_DIR, "LibTemplate.jl"))
