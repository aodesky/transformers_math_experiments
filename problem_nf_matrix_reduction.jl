include("constants.jl")
using LinearAlgebra

"""
Number Field Matrix Reduction Problem

Dataset: Degree 3 Galois number fields with Minkowski embeddings
Goal: Minimize ||AM||_∞ where:
  - M is the 3×3 matrix [[1, u1, u1²], [1, u2, u2²], [1, u3, u3²]]
  - u1, u2, u3 are the three Minkowski embeddings
  - A is a 3×3 invertible matrix with integer entries
  - ||·||_∞ is the max absolute value of matrix coefficients

Local search: Multiply A on the left by random 3×3 invertible integer matrices
"""

# Load the number field dataset
const NF_DATA_PATH = joinpath(dirname(@__FILE__), "..", "nf_deg_3_galois.csv")

struct NumberField
    label::String
    disc_abs::Int
    index::Int
    u1::Float64
    u2::Float64
    u3::Float64
    v1::Float64  # Randomized embedding: v = R * u
    v2::Float64
    v3::Float64
    R::Matrix{Int64}  # The random transformation matrix used (circulant form)
end

function generate_circulant_matrix()::Matrix{Int64}
    """
    Generate a random invertible circulant matrix of the form:
    R = [[a, b, c],
         [c, a, b],
         [b, c, a]]

    For this matrix, det(R) = a³ + b³ + c³ - 3abc
    We need det(R) ≠ 0 for invertibility
    """
    max_attempts = 100
    for _ in 1:max_attempts
        a = rand(-5:5)
        b = rand(-5:5)
        c = rand(-5:5)

        # Compute determinant: a³ + b³ + c³ - 3abc
        det_val = a^3 + b^3 + c^3 - 3*a*b*c

        if det_val != 0
            R = [a b c;
                 c a b;
                 b c a]
            return R
        end
    end

    # Fallback: return identity if we can't find an invertible one
    # (this shouldn't happen in practice)
    return Matrix{Int64}(I, 3, 3)
end

function load_csv_data(filepath::String)::Vector{NumberField}
    """
    Simple CSV parser for the number field dataset.
    For each field, generates a random circulant matrix R and computes v = R * u
    """
    fields = NumberField[]

    open(filepath, "r") do f
        # Skip header line
        readline(f)

        # Read data lines
        for line in eachline(f)
            parts = split(line, ",")
            if length(parts) != 6
                continue
            end

            try
                label = parts[1]
                disc_abs = parse(Int, parts[2])
                index = parse(Int, parts[3])
                u1 = parse(Float64, parts[4])
                u2 = parse(Float64, parts[5])
                u3 = parse(Float64, parts[6])

                # Generate random circulant matrix R
                R = generate_circulant_matrix()

                # Compute v = R * u
                u = [u1, u2, u3]
                v = R * u

                push!(fields, NumberField(
                    label, disc_abs, index,
                    u1, u2, u3,  # Keep original embeddings for reference
                    v[1], v[2], v[3],  # Randomized embeddings
                    R  # Store the transformation matrix
                ))
            catch e
                println("Warning: Failed to parse line: $line")
            end
        end
    end

    return fields
end

# Load dataset at module initialization
const NUMBER_FIELDS = let
    if !isfile(NF_DATA_PATH)
        error("Dataset not found at $NF_DATA_PATH")
    end

    fields = load_csv_data(NF_DATA_PATH)
    println("Loaded $(length(fields)) number fields from dataset")
    fields
end

function construct_M_matrix(nf::NumberField)::Matrix{Float64}
    """
    Construct the 3×3 matrix M = [[1, v1, v1²], [1, v2, v2²], [1, v3, v3²]]
    where v = R * u are the randomized embeddings
    """
    v1, v2, v3 = nf.v1, nf.v2, nf.v3
    return [
        1.0  v1     v1^2;
        1.0  v2     v2^2;
        1.0  v3     v3^2
    ]
end

function parse_object(obj::OBJ_TYPE)::Tuple{Int, Matrix{Int64}}
    """
    Parse object string "field_idx:a11,a12,a13,a21,a22,a23,a31,a32,a33"
    Returns (field_index, A_matrix)
    """
    parts = split(obj, ":")
    if length(parts) != 2
        error("Invalid object format: $obj")
    end

    field_idx = parse(Int, parts[1])
    matrix_str = parts[2]

    # Parse the 9 matrix entries
    entries = [parse(Int64, s) for s in split(matrix_str, ",")]
    if length(entries) != 9
        error("Expected 9 matrix entries, got $(length(entries))")
    end

    # Construct 3×3 matrix (row-major order)
    A = reshape(entries, 3, 3)'

    return (field_idx, A)
end

function object_to_string(field_idx::Int, A::Matrix{Int64})::OBJ_TYPE
    """
    Convert (field_index, A_matrix) to string representation
    """
    # Flatten matrix in row-major order
    entries = vec(A')
    matrix_str = join(entries, ",")
    return "$field_idx:$matrix_str"
end

function matrix_max_norm(M::Matrix{Float64})::Float64
    """
    Compute max absolute value of matrix entries
    """
    return maximum(abs.(M))
end

function generate_random_invertible_matrix()::Matrix{Int64}
    """
    Generate a random 3×3 invertible matrix with integer entries
    Uses elementary row operations to ensure invertibility
    """
    # Start with identity matrix
    A = Matrix{Int64}(I, 3, 3)

    # Apply random elementary row operations
    num_operations = rand(3:6)

    for _ in 1:num_operations
        op_type = rand(1:3)
        i, j = rand(1:3), rand(1:3)

        if op_type == 1 && i != j
            # Type I: Add a multiple of row j to row i
            multiplier = rand(-3:3)
            A[i, :] .+= multiplier .* A[j, :]
        elseif op_type == 2
            # Type II: Multiply row i by -1 (keeps determinant ±1)
            A[i, :] .*= -1
        elseif op_type == 3 && i != j
            # Type III: Swap rows i and j
            A[i, :], A[j, :] = A[j, :], A[i, :]
        end
    end

    # Verify determinant is ±1
    det_val = round(Int, det(A))
    if abs(det_val) != 1
        # Fallback: return a matrix we know is invertible
        return Matrix{Int64}(I, 3, 3)
    end

    return A
end

function greedy_search_from_startpoint(db, obj::OBJ_TYPE)::Vector{OBJ_TYPE}
    """
    Local search: Generate multiple candidate matrices by multiplying A
    on the left by random 3×3 invertible integer matrices

    Returns several candidates to explore different search directions
    """
    try
        field_idx, A_current = parse_object(obj)

        # Get the number field and its M matrix
        nf = NUMBER_FIELDS[field_idx]
        M = construct_M_matrix(nf)

        # Generate multiple random transformations
        results = Vector{OBJ_TYPE}()
        num_samples = 10

        for _ in 1:num_samples
            # Generate random invertible matrix R
            R = generate_random_invertible_matrix()

            # Compute new A' = R * A
            A_new = R * A_current

            # Add to results
            push!(results, object_to_string(field_idx, A_new))
        end

        return results
    catch e
        println("Error in greedy_search_from_startpoint: $e")
        return [empty_starting_point()]
    end
end

function reward_calc(obj::OBJ_TYPE)::REWARD_TYPE
    """
    Compute reward = -||AM||_∞

    Negative because we want to minimize the norm, but the framework maximizes reward
    """
    try
        field_idx, A = parse_object(obj)

        # Get the number field and construct M
        nf = NUMBER_FIELDS[field_idx]
        M = construct_M_matrix(nf)

        # Compute AM
        AM = Float64.(A) * M

        # Compute max norm
        norm_val = matrix_max_norm(AM)

        # Return negative (we want to minimize)
        return Float32(-norm_val)
    catch e
        println("Error computing reward for $obj: $e")
        return Float32(-1e9)
    end
end

function empty_starting_point()::OBJ_TYPE
    """
    Generate a random starting point:
    - Pick a random number field from the dataset
    - Start with the identity matrix
    """
    field_idx = rand(1:length(NUMBER_FIELDS))
    A_identity = Matrix{Int64}(I, 3, 3)
    return object_to_string(field_idx, A_identity)
end

# Print configuration
println()
println("="^80)
println("Number Field Matrix Reduction Problem")
println("="^80)
println("Dataset: $(length(NUMBER_FIELDS)) degree 3 Galois number fields")
println("Goal: Minimize ||AM||_∞ where M is constructed from Minkowski embeddings")
println("Local search: Multiply A by random 3×3 invertible integer matrices")
println("="^80)
println()
