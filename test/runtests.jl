using HiGHS_jll

@show run(`uname -m`)
@show run(`which valgrind`)
@show run(`ls /snap/bin/valgrind`)
@show run(`valgrind ls`)

exe = HiGHS_jll.highs_path
@show run(`valgrind $(exe) fail-on-32-bit.mps`)

# valgrind = addenv(`valgrind`, HiGHS_jll.JLLWrappers.JLLWrappers.LIBPATH_env=>HiGHS_jll.LIBPATH[]);
# HiGHS_jll.highs() do exe
#     run(`valgrind $(exe) fail-on-32-bit.mps`)
# end

# using Test

# @testset "$(file)" for file in readdir(@__DIR__)
#     if file == "runtests.jl" || !endswith(file, ".jl")
#         continue
#     end
#     include(file)
# end
