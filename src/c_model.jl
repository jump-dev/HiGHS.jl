
mutable struct ManagedHiGHS
    inner::Ptr{Cvoid}

    function ManagedHiGHS()
        mghs = new(
            CWrapper.Highs_create()
        )
        # register the memory cleanup function
        finalizer(mghs) do m
            success = free_highs(m)
            if !success
                @warn "Memory free failure, possible leak."
            end
        end
    end
end

"""
Release references and free memory.
Returns a boolean indicating the success of the operation.
"""
function free_highs(mhgs::ManagedHiGHS)
    # Avoid double-free (SCIP will set the pointers to NULL).
    if mhgs.inner == C_NULL
        return false
    end
    # only mscip.scip is GC-protected during ccall!
    GC.@preserve mhgs begin
        CWrapper.Highs_destroy(mhgs.inner)
    end
    mhgs.inner = C_NULL
    return true
end
