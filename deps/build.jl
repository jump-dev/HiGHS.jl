using BinaryProvider

# Parse some basic command-line arguments
const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))

# These are the two binary objects we care about
const products = Product[LibraryProduct(prefix, "libhighs", :libhighs)]

# Download binaries from hosted location
bin_prefix = "https://github.com/matbesancon/HiGHSBuilder/releases/download/v0.1.0"

download_info = Dict(
    Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(:gcc4))  => ("$bin_prefix/highs.v0.1.0.x86_64-linux-gnu-gcc4.tar.gz", "a0df6e2899f2ba5aadfbd5d865136abf31e41f768156038bca72f859dbed5373"),
    Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(:gcc7))  => ("$bin_prefix/highs.v0.1.0.x86_64-linux-gnu-gcc7.tar.gz", "9b2d37d473126f3394a131e718f22114a892effc0e23c7471a3a242c6122a539"),
    Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(:gcc8))  => ("$bin_prefix/highs.v0.1.0.x86_64-linux-gnu-gcc8.tar.gz", "277ddbc459941a1cd4ce052c16008ccc65c379cdee34377eda70f1ecfe1d48ac"),
)



# First, check to see if we're all satisfied
if any(!satisfied(p; verbose=verbose) for p in products)
    try
        # Download and install binaries
        url, tarball_hash = choose_download(download_info)
        install(url, tarball_hash; prefix=prefix, force=true, verbose=true)
    catch e
        if typeof(e) <: ArgumentError
            error("Your platform $(Sys.MACHINE) is not supported by this package!")
        else
            rethrow(e)
        end
    end

    # Finally, write out a deps.jl file
    write_deps_file(joinpath(@__DIR__, "deps.jl"), products)
end
