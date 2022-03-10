using HiGHS_jll

@show run(`uname -m`)
@show run(`which valgrind`)
@show run(`valgrind ls`)

exe = HiGHS_jll.highs_path
run(`$(exe) --parallel=off --presolve=off fail-on-32-bit.mps`)
run(`$(exe) --parallel=off fail-on-32-bit.mps`)
run(`$(exe) fail-on-32-bit.mps`)
@show run(`valgrind --leak-check=full $(exe) --parallel=off --presolve=off fail-on-32-bit.mps`)
# @show run(`export VALGRIND_LIB="/usr/local/lib/valgrind"; valgrind --leak-check=full $(exe) fail-on-32-bit.mps`)

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
