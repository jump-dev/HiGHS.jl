import Clang

const HIGHS_DIR = haskey(ENV, "HIGHS_PATH") ? ENV["HIGHS_PATH"] : joinpath(@__DIR__, "../deps/usr")

const header_file = joinpath(HIGHS_DIR, "src/interfaces/highs_c_api.h")

const LIB_HEADERS = [header_file]

const ctx = Clang.DefaultContext()

Clang.parse_headers!(ctx, LIB_HEADERS,
    includes=[Clang.CLANG_INCLUDE],
)

ctx.libname = "libhighs"
ctx.options["is_function_strictly_typed"] = true
ctx.options["is_struct_mutable"] = false

const api_file = joinpath(@__DIR__, "../src/wrapper", "$(ctx.libname)_api.jl")

open(api_file, "w") do f
    for trans_unit in ctx.trans_units
        root_cursor = Clang.getcursor(trans_unit)
        push!(ctx.cursor_stack, root_cursor)
        header = Clang.spelling(root_cursor)
        @info "wrapping header: $header ..."
        # loop over all of the child cursors and wrap them, if appropriate.
        ctx.children = Clang.children(root_cursor)
        for (i, child) in enumerate(ctx.children)
           child_name = Clang.name(child)
           child_header = Clang.filename(child)
           ctx.children_index = i
           # choose which cursor to wrap
           startswith(child_name, "__") && continue  # skip compiler definitions
           child_name in keys(ctx.common_buffer) && continue  # already wrapped
           child_header != header && continue  # skip if cursor filename is not in the headers to be wrapped
           Clang.wrap!(ctx, child)
        end
        @info "writing $(api_file)"
        println(f, "# Julia wrapper for header: $(basename(header))")
        println(f, "# Automatically generated using Clang.jl\n")
        Clang.print_buffer(f, ctx.api_buffer)
        empty!(ctx.api_buffer)  # clean up api_buffer for the next header
    end
end

const common_file = joinpath(@__DIR__, "../src/wrapper", "$(ctx.libname)_common.jl")

open(common_file, "w") do f
    println(f, "# Automatically generated using Clang.jl\n")
    Clang.print_buffer(f, Clang.dump_to_buffer(ctx.common_buffer))
end
