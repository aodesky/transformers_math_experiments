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

            # Just write the object (it already contains A, disc, index, M)
            # No need to write features - those are computed later
            write(f, obj)
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
