using BinaryProvider

# Parse some basic command-line arguments
const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))

# These are the two binary objects we care about
const products = Product[LibraryProduct(prefix, "libhighs", :libhighs)]

# Download binaries from hosted location
const bin_prefix = "https://github.com/matbesancon/HiGHSBuilder/releases/download/v0.1.1"

const download_info = Dict(
    Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(:gcc4))  => ("$bin_prefix/highs.v0.1.1.x86_64-linux-gnu-gcc4.tar.gz", "efc83bcd2052b46b8148a866fcafb88d161f0a0cd5735090586f220e62c28035"),
    Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(:gcc7))  => ("$bin_prefix/highs.v0.1.1.x86_64-linux-gnu-gcc7.tar.gz", "6b3f2f2d1405e3f73ed0702607966f9b86fb2e0b0edee1cae6f9b9c9552e09eb"),
    Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(:gcc8))  => ("$bin_prefix/highs.v0.1.1.x86_64-linux-gnu-gcc8.tar.gz", "97f275fb5dbe729fd29c3d74b9a3c6cc572a8393345943328f3a5de8a267c057"),
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
