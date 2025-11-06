"""
Test script for the elliptic curve conductor minimization problem
"""

include("constants.jl")
include("problem_elliptic_simple.jl")

println("="^80)
println("Testing Elliptic Curve Problem")
println("="^80)

# Test 1: Parsing rationals
println("\n1. Testing rational parsing:")
test_vals = ["1/2", "3/4", "-5/2", "7"]
for s in test_vals
    num, den = parse_rational(s)
    println("  $s -> numerator=$num, denominator=$den")
    println("  Back to string: $(rational_to_string(num, den))")
end

# Test 2: Discriminant computation
println("\n2. Testing discriminant D(t) = -16(4t³ + 27):")
test_rationals = [(1, 1), (0, 1), (2, 1), (1, 2), (-1, 1)]
for (num, den) in test_rationals
    disc = eval_discriminant(num, den)
    println("  t=$num/$den: D(t) = $disc")
end

# Test 3: Discriminant derivative
println("\n3. Testing discriminant derivative D'(t) = -192t²:")
for (num, den) in test_rationals
    deriv = eval_discriminant_derivative(num, den)
    println("  t=$num/$den: D'(t) = $deriv")
end

# Test 4: Check if GP/Pari is available
println("\n4. Testing GP/Pari availability:")
try
    result = read(pipeline(`gp --version`), String)
    println("  GP/Pari found:")
    println("  ", split(result, "\n")[1])
catch e
    println("  ERROR: GP/Pari not found! Please install GP/Pari.")
    println("  Error: ", e)
end

# Test 5: Conductor computation (if GP/Pari works)
println("\n5. Testing conductor computation:")
test_simple = [(1, 1), (2, 1), (0, 1)]
for (num, den) in test_simple
    disc = eval_discriminant(num, den)
    if disc != 0
        try
            cond = compute_conductor(num, den)
            println("  t=$num/$den: conductor = $cond")
        catch e
            println("  t=$num/$den: ERROR computing conductor: ", e)
        end
    else
        println("  t=$num/$den: SKIP (discriminant = 0)")
    end
end

# Test 6: Gradient descent step
println("\n6. Testing gradient descent:")
obj = "1/1"
println("  Starting point: $obj")
new_objs = greedy_search_from_startpoint(nothing, obj)
println("  After gradient descent: $(new_objs[1])")

# Test 7: Reward computation
println("\n7. Testing reward function:")
test_objs = ["1/1", "2/1", "1/2"]
for obj in test_objs
    reward = reward_calc(obj)
    println("  $obj: reward = $reward")
end

# Test 8: Empty starting point
println("\n8. Testing empty starting point:")
start = empty_starting_point()
println("  Empty starting point: $start")

println("\n" * "="^80)
println("Testing complete!")
println("="^80)
