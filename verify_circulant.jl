#!/usr/bin/env julia

include("problem_nf_matrix_reduction.jl")

# Check first few number fields
println("Verifying circulant transformation for first 3 fields:\n")

for i in 1:3
    nf = NUMBER_FIELDS[i]

    println("Field $i: $(nf.label)")
    println("  Original u = ($(nf.u1), $(nf.u2), $(nf.u3))")
    println("  Circulant matrix R:")
    for row in eachrow(nf.R)
        println("    ", row)
    end
    println("  det(R) = $(det(nf.R))")

    # Verify v = R * u
    u = [nf.u1, nf.u2, nf.u3]
    v_computed = nf.R * u
    v_stored = [nf.v1, nf.v2, nf.v3]

    println("  Computed v = R*u = ($(v_computed[1]), $(v_computed[2]), $(v_computed[3]))")
    println("  Stored v = ($(nf.v1), $(nf.v2), $(nf.v3))")
    println("  Match: ", isapprox(v_computed, v_stored, atol=1e-10))
    println()
end
