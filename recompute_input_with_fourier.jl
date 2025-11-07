#!/usr/bin/env julia

"""
Recompute elliptic_input.txt with Fourier coefficients (ellap values)

Reads: elliptic_input.txt
Writes: elliptic_input_with_fourier_coefficients.txt

Format: num/den,a_2,a_3,a_5,...,a_97
"""

include("problem_mestre_rank12.jl")

function main()
    input_file = "elliptic_input.txt"
    output_file = "elliptic_input_with_fourier_coefficients.txt"

    if !isfile(input_file)
        println("Error: $input_file not found")
        return
    end

    println("="^80)
    println("Recomputing elliptic_input.txt with Fourier coefficients")
    println("="^80)
    println()
    println("Reading from: $input_file")
    println("Writing to:   $output_file")
    println()

    # Read input rationals
    rationals = String[]
    open(input_file, "r") do f
        for line in eachline(f)
            line = strip(line)
            if !isempty(line)
                push!(rationals, line)
            end
        end
    end

    n_total = length(rationals)
    println("Found $n_total rationals to process")
    println()

    # Process and write output
    n_processed = 0
    n_errors = 0

    open(output_file, "w") do out
        for (i, rational_str) in enumerate(rationals)
            try
                # Parse rational
                num, den = parse_rational(rational_str)

                # Compute ellap values
                ellap_values = compute_ellap_batch(num, den)

                # Write: rational,a_2,a_3,...,a_97
                write(out, rational_str)
                for ap in ellap_values
                    write(out, ",$ap")
                end
                write(out, "\n")

                n_processed += 1

                # Progress indicator
                if i % 100 == 0
                    println("Processed $i/$n_total ($(round(100*i/n_total, digits=1))%)")
                end
            catch e
                println("Error processing $rational_str: $e")
                n_errors += 1
            end
        end
    end

    println()
    println("="^80)
    println("COMPLETE")
    println("="^80)
    println("Processed: $n_processed")
    println("Errors:    $n_errors")
    println("Output:    $output_file")
    println()
end

main()
