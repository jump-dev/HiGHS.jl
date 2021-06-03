if get(ENV, "GITHUB_ACTIONS", "") == "true"
    import Pkg
    Pkg.add(Pkg.PackageSpec(name = "MathOptInterface", rev = "master"))
end

using Test

@testset "$(file)" for file in readdir(@__DIR__)
    if file == "runtests.jl" || !endswith(file, ".jl")
        continue
    end
    include(file)
end
