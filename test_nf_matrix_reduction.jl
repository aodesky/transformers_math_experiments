#!/usr/bin/env julia

"""
Test script for the Number Field Matrix Reduction problem
"""

include("problem_nf_matrix_reduction.jl")

println("\n" * "="^80)
println("Testing Number Field Matrix Reduction Problem")
println("="^80)

# Test 1: Generate an empty starting point
println("\n[Test 1] Generate empty starting point")
obj1 = empty_starting_point()
println("Starting point: $obj1")

# Parse it back
field_idx, A = parse_object(obj1)
println("Field index: $field_idx")
println("Matrix A:\n$A")

# Test 2: Compute reward for starting point
println("\n[Test 2] Compute reward for starting point")
reward = reward_calc(obj1)
println("Reward (negative norm): $reward")

# Show the actual norm
nf = NUMBER_FIELDS[field_idx]
M = construct_M_matrix(nf)
println("M matrix:\n$M")
AM = Float64.(A) * M
println("AM matrix:\n$AM")
println("||AM||_∞ = $(matrix_max_norm(AM))")

# Test 3: Perform local search
println("\n[Test 3] Perform local search")
candidates = greedy_search_from_startpoint(nothing, obj1)
println("Generated $(length(candidates)) candidates")

# Evaluate each candidate
println("\nEvaluating candidates:")
for (i, cand) in enumerate(candidates[1:min(5, length(candidates))])
    r = reward_calc(cand)
    field_idx2, A2 = parse_object(cand)
    AM2 = Float64.(A2) * M
    norm_val = matrix_max_norm(AM2)
    println("  Candidate $i: reward = $r, ||AM||_∞ = $norm_val")
end

# Test 4: Test random invertible matrix generation
println("\n[Test 4] Test random invertible matrix generation")
for i in 1:5
    R = generate_random_invertible_matrix()
    det_val = det(R)
    println("  Matrix $i: det = $det_val")
    if abs(det_val) ≈ 1.0
        println("    ✓ Invertible (det ≈ ±1)")
    else
        println("    ✗ ERROR: Not invertible!")
    end
end

# Test 5: Run several iterations of improvement
println("\n[Test 5] Run improvement iterations")
function run_improvement_test()
    current = empty_starting_point()
    current_reward = reward_calc(current)
    println("Initial: reward = $current_reward")

    for iter in 1:10
        candidates = greedy_search_from_startpoint(nothing, current)

        # Find best candidate
        best_cand = current
        best_reward = current_reward

        for cand in candidates
            r = reward_calc(cand)
            if r > best_reward
                best_cand = cand
                best_reward = r
            end
        end

        if best_reward > current_reward
            current = best_cand
            current_reward = best_reward
            println("Iteration $iter: improved to reward = $current_reward (norm = $(-current_reward))")
        else
            println("Iteration $iter: no improvement")
        end
    end

    field_idx, A_final = parse_object(current)
    nf = NUMBER_FIELDS[field_idx]
    M = construct_M_matrix(nf)
    AM_final = Float64.(A_final) * M
    println("\nFinal matrix A:\n$A_final")
    println("Final AM:\n$AM_final")
    println("Final ||AM||_∞ = $(matrix_max_norm(AM_final))")
end

run_improvement_test()

println("\n" * "="^80)
println("All tests completed!")
println("="^80 * "\n")
