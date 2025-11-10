# Number Field Matrix Reduction Problem

## Problem Description

**Objective**: Minimize ||AM||∞ (max absolute value of matrix entries) where:
- **M** is a 3×3 matrix constructed from Minkowski embeddings of a degree 3 Galois number field
- **A** is a 3×3 invertible matrix with integer entries (what we optimize)

## How It Works

1. **Load dataset**: 4,679 degree 3 Galois number fields with Minkowski embeddings u = (u1, u2, u3)

2. **Randomize generators**: Apply circulant transformation v = R·u where
   ```
   R = [[a, b, c],
        [c, a, b],
        [b, c, a]]
   ```
   with random integers a, b, c (det(R) ≠ 0)

3. **Construct M** (Vandermonde form using randomized embeddings):
   ```
   M = [[1,  v1,  v1²],
        [1,  v2,  v2²],
        [1,  v3,  v3²]]
   ```

4. **Local search**: Multiply A on the left by random 3×3 invertible integer matrices (det = ±1)

5. **Reward**: `-||AM||∞` (negative because we minimize)

## Quick Start

### Test the implementation:
```bash
julia test_nf_matrix_reduction.jl
```

### Verify circulant randomization:
```bash
julia verify_circulant.jl
```

### Generate input dataset:
```bash
julia generate_nf_input.jl nf_matrix_input.txt 3000
```

### Run local search:
```bash
# Edit search_fc.jl line 18: include("problem_nf_matrix_reduction.jl")
julia search_fc.jl ./output 120 1000 1000 500 -i nf_matrix_input.txt
```

## Encoding

States are encoded as: `"field_idx:a11,a12,a13,a21,a22,a23,a31,a32,a33"`

Example: `"951:1,0,0,0,1,0,0,0,1"` (field 951, identity matrix)

## Example Results

| Field | Initial ||AM||∞ | After 10 Iterations | Reduction |
|-------|-----------------|---------------------|-----------|
| 951   | 131,174.9       | ~1,746.7           | ~98.7%    |
| Random| 8,483.3         | ~1,746.7           | ~79.4%    |

## Files

- `problem_nf_matrix_reduction.jl`: Main problem definition
- `test_nf_matrix_reduction.jl`: Test suite
- `generate_nf_input.jl`: Generate input dataset
- `verify_circulant.jl`: Verify randomization
- `../nf_deg_3_galois.csv`: Number field dataset (4,679 fields)

## Key Functions

- `empty_starting_point()`: Random field + identity matrix
- `greedy_search_from_startpoint(db, obj)`: Multiply A by random matrices
- `reward_calc(obj)`: Compute `-||AM||∞`
- `generate_circulant_matrix()`: Create randomizing transformation R
- `construct_M_matrix(nf)`: Build Vandermonde matrix from v
