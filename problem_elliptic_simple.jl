include("constants.jl")

"""
Simple Elliptic Curve Conductor Minimization Problem

Elliptic curve family: E(t): y² = x³ + t·x + 1
- D(t) = discriminant = -16(4t³ + 27)
- D'(t) = -192t²
- C(t) = conductor of E(t) (computed via GP/Pari)

Goal: Minimize C(t) where t is a rational number
Constraint: D(t) ≠ 0

Local search: Single step of gradient descent on D(t)
Learning rate: 1/50
"""

# Helper functions for rational arithmetic
function parse_rational(s::String)::Tuple{Int64, Int64}
    """Parse a rational string 'num/den' into (numerator, denominator)"""
    if occursin("/", s)
        parts = split(s, "/")
        num = parse(Int64, parts[1])
        den = parse(Int64, parts[2])
        return (num, den)
    else
        # If no slash, treat as integer
        num = parse(Int64, s)
        return (num, 1)
    end
end

function rational_to_string(num::Int64, den::Int64)::String
    """Convert rational to string, simplifying if needed"""
    if den == 1
        return string(num)
    else
        return string(num) * "/" * string(den)
    end
end

function simplify_rational(num::Int64, den::Int64)::Tuple{Int64, Int64}
    """Simplify a rational number using GCD"""
    if den < 0
        num = -num
        den = -den
    end
    if num == 0
        return (0, 1)
    end
    g = gcd(abs(num), abs(den))
    return (div(num, g), div(den, g))
end

function eval_discriminant(num::Int64, den::Int64)::Rational{Int64}
    """
    Evaluate D(t) = -16(4t³ + 27) at t = num/den
    Returns the discriminant as a rational number
    """
    # t = num/den
    # t³ = num³/den³
    # 4t³ = 4*num³/den³
    # 4t³ + 27 = (4*num³ + 27*den³)/den³
    # -16(4t³ + 27) = -16*(4*num³ + 27*den³)/den³

    numer = -16 * (4 * num^3 + 27 * den^3)
    denom = den^3

    return numer // denom
end

function eval_discriminant_derivative(num::Int64, den::Int64)::Rational{Int64}
    """
    Evaluate D'(t) = -192t² at t = num/den
    Returns the derivative as a rational number
    """
    # t² = num²/den²
    # -192t² = -192*num²/den²

    numer = -192 * num^2
    denom = den^2

    return numer // denom
end

function compute_conductor(num::Int64, den::Int64)::Int64
    """
    Compute the conductor of E(t): y² = x³ + t·x + 1
    where t = num/den, using GP/Pari

    The curve in Weierstrass form [a1,a2,a3,a4,a6] is [0,0,0,t,1]
    """
    t_val = num // den  # Rational number

    # Create GP/Pari script to compute conductor
    # For curve y² = x³ + a4*x + a6, we have [0,0,0,a4,a6]
    gp_script = """
    t = $num/$den;
    E = ellinit([0, 0, 0, t, 1]);
    if (E == 0, print("INVALID"), red = ellglobalred(E); print(red[1]));
    quit();
    """

    # Write script to temp file
    temp_file = tempname() * ".gp"
    open(temp_file, "w") do f
        write(f, gp_script)
    end

    try
        # Run GP/Pari
        result = read(pipeline(`gp -q $temp_file`), String)
        rm(temp_file)

        result = strip(result)
        if result == "INVALID"
            return typemax(Int64)  # Return large value for invalid curves
        end

        return parse(Int64, result)
    catch e
        # If GP/Pari fails, return a large conductor
        println("Warning: GP/Pari failed for t=$num/$den: ", e)
        if isfile(temp_file)
            rm(temp_file)
        end
        return typemax(Int64)
    end
end

function greedy_search_from_startpoint(db, obj::OBJ_TYPE)::Vector{OBJ_TYPE}
    """
    Perform one step of gradient descent on D(t)

    Input: obj is a string "num/den" representing rational t
    Output: Vector containing the new rational after gradient descent

    Gradient descent rule:
    - If D(t) > 0: t_new = t - learning_rate * D'(t)
    - If D(t) < 0: t_new = t + learning_rate * D'(t)
    - If D(t) = 0: invalid, return original
    """

    # Parse the input rational
    local num, den
    try
        num, den = parse_rational(obj)
    catch
        # If parsing fails, return empty starting point
        return [empty_starting_point()]
    end

    # Evaluate discriminant
    disc = eval_discriminant(num, den)

    # Check if discriminant is zero (invalid curve)
    if disc == 0
        # Return a random nearby rational
        new_num, new_den = simplify_rational(num + rand(-5:5), den + max(1, rand(-2:2)))
        return [rational_to_string(new_num, new_den)]
    end

    # Evaluate derivative
    disc_deriv = eval_discriminant_derivative(num, den)

    # Learning rate: 1/50
    learning_rate = 1 // 50

    # Gradient descent step (accounting for sign of D)
    # t_new = t - sign(D) * learning_rate * D'(t)
    t_current = num // den

    if disc > 0
        t_new = t_current - learning_rate * disc_deriv
    else  # disc < 0
        t_new = t_current + learning_rate * disc_deriv
    end

    # Convert back to simplified rational
    new_num = numerator(t_new)
    new_den = denominator(t_new)
    new_num, new_den = simplify_rational(Int64(new_num), Int64(new_den))

    # Verify the new point has non-zero discriminant
    new_disc = eval_discriminant(new_num, new_den)
    if new_disc == 0
        # If we landed on invalid point, try a small perturbation
        new_num += 1
        new_num, new_den = simplify_rational(new_num, new_den)
    end

    return [rational_to_string(new_num, new_den)]
end

function reward_calc(obj::OBJ_TYPE)::REWARD_TYPE
    """
    Compute the reward = -log(conductor(t))
    (Negative because we want to minimize log(conductor), but the framework maximizes reward)
    """
    try
        num, den = parse_rational(obj)

        # Check discriminant is non-zero
        disc = eval_discriminant(num, den)
        if disc == 0
            return Float32(-1e9)  # Very bad reward for invalid curves
        end

        conductor = compute_conductor(num, den)

        if conductor <= 0 || conductor == typemax(Int64)
            return Float32(-1e9)  # Invalid conductor
        end

        # Return negative log of conductor (we want to minimize it)
        return Float32(-log(Float64(conductor)))
    catch e
        println("Error computing reward for $obj: ", e)
        return Float32(-1e9)
    end
end

function empty_starting_point()::OBJ_TYPE
    """
    Initial starting point: t = 1/1
    """
    return "1/1"
end
