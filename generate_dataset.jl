include("problem_mestre_rank12.jl")
using Random

"""
Generate a dataset of elliptic curves from Mestre's family

Parameters:
- target_size: number of valid curves to generate
- max_height: maximum height (max(|num|, |den|)) for rational points
"""
function generate_dataset(target_size::Int, max_height::Int)
    println("Generating dataset with $target_size curves (height < $max_height)...")
    println("Using $(Threads.nthreads()) threads")
    println("="^80)

    # Pre-generate candidate rationals (more than needed to account for filtering)
    Random.seed!(42)  # For reproducibility
    n_candidates = target_size * 3  # Generate 3x more candidates to account for filtering

    candidates = Vector{Tuple{Int64, Int64}}()
    while length(candidates) < n_candidates
        num = rand(-max_height+1:max_height-1)
        den = rand(1:max_height-1)

        if den != 0
            num, den = simplify_rational(num, den)
            push!(candidates, (num, den))
        end
    end

    println("Generated $(length(candidates)) candidate rationals")
    println()

    # Process candidates in parallel
    dataset = Vector{String}()
    dataset_lock = ReentrantLock()

    skipped_disc_zero = Threads.Atomic{Int}(0)
    skipped_conductor_fail = Threads.Atomic{Int}(0)
    processed = Threads.Atomic{Int}(0)

    Threads.@threads for (num, den) in candidates
        # Stop if we have enough curves
        if length(dataset) >= target_size
            break
        end

        # Check discriminant is non-zero
        disc = eval_discriminant(num, den)
        if abs(disc) < 1e-10
            Threads.atomic_add!(skipped_disc_zero, 1)
            continue
        end

        # Compute ap values (skip conductor check for speed)
        ap_values = compute_ellap_batch(num, den)

        # Check if all zeros (indicates error)
        if all(x -> x == 0, ap_values)
            Threads.atomic_add!(skipped_conductor_fail, 1)
            continue
        end

        # Create entry: t,a_2,a_3,a_5,...,a_97
        t_str = rational_to_string(num, den)
        ap_str = join(ap_values, ",")
        entry = "$t_str,$ap_str"

        # Thread-safe append
        lock(dataset_lock) do
            if length(dataset) < target_size
                push!(dataset, entry)

                # Progress update every 100 entries
                if length(dataset) % 100 == 0
                    println("Progress: $(length(dataset))/$target_size curves generated")
                end
            end
        end

        Threads.atomic_add!(processed, 1)
    end

    println()
    println("="^80)
    println("Dataset generation complete!")
    println("Total curves: $(length(dataset))")
    println("Total processed: $(processed[])")
    println("Skipped (disc=0): $(skipped_disc_zero[])")
    println("Skipped (ap computation fail): $(skipped_conductor_fail[])")
    println("="^80)

    return dataset
end

# Generate dataset
println("Starting dataset generation...")
println()

dataset = generate_dataset(3000, 100)

# Write to file
output_file = "elliptic_input_with_fourier_coefficients.txt"
println()
println("Writing to $output_file...")
open(output_file, "w") do f
    for entry in dataset
        println(f, entry)
    end
end

println("✓ Done! Dataset written to $output_file")
println("File size: ", stat(output_file).size ÷ 1024, " KB")
