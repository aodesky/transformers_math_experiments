#!/usr/bin/env julia

"""
Generate initial dataset for number field matrix reduction problem
Creates a file with format: field_idx:matrix,features...
"""

include("problem_nf_matrix_reduction.jl")

function main()
    output_file = ARGS[1]
    num_samples = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 3000

    println("Generating $num_samples initial samples...")

    open(output_file, "w") do f
        for i in 1:num_samples
            # Generate a random starting point
            obj = empty_starting_point()

            # Compute some basic features of the current state
            field_idx, A = parse_object(obj)
            nf = NUMBER_FIELDS[field_idx]
            M = construct_M_matrix(nf)
            AM = Float64.(A) * M

            # Features: matrix norm, discriminant, index, and first few entries of AM
            norm_val = matrix_max_norm(AM)
            disc = nf.disc_abs
            index = nf.index

            # Flatten AM and take some entries as features
            AM_flat = vec(AM')

            # Write: object,norm,disc,index,AM_entries...
            write(f, obj)
            write(f, ",$norm_val")
            write(f, ",$disc")
            write(f, ",$index")

            # Add first 6 entries of AM as additional features
            for val in AM_flat[1:min(6, length(AM_flat))]
                write(f, ",$(round(val, digits=2))")
            end

            write(f, "\n")

            if i % 500 == 0
                println("  Generated $i samples...")
            end
        end
    end

    println("Generated $num_samples samples to $output_file")
end

if length(ARGS) < 1
    println("Usage: julia generate_nf_input.jl <output_file> [num_samples]")
    println("Example: julia generate_nf_input.jl nf_matrix_input.txt 3000")
    exit(1)
end

main()
