# Elliptic Curve Conductor Minimization Problem

This document describes the elliptic curve conductor minimization problem added to the PatternBoost framework.

## Problem Description

### Objective
Minimize the conductor C(t) of an elliptic curve E(t) over rational numbers t, subject to the constraint that the discriminant D(t) ≠ 0.

### Elliptic Curve Family
We use a simplified test family:
```
E(t): y² = x³ + t·x + 1
```

Where:
- **t** is a rational number (represented as `numerator/denominator`)
- **C(t)** is the conductor of E(t) (computed via GP/Pari)
- **D(t)** is the discriminant: `D(t) = -16(4t³ + 27)`
- **D'(t)** is the derivative: `D'(t) = -192t²`

### Constraints
- **Valid curves**: D(t) ≠ 0 (discriminant must be non-zero)
- **Reward function**: `-C(t)` (negative conductor, since we want to minimize)

### Local Search Algorithm

The local search performs **one step of gradient descent** on the discriminant D(t):

1. **Parse** the input rational t = num/den
2. **Evaluate** D(t) and D'(t)
3. **Check validity**: If D(t) = 0, return a perturbed rational
4. **Gradient descent step**:
   - If D(t) > 0: `t_new = t - lr × D'(t)`
   - If D(t) < 0: `t_new = t + lr × D'(t)`
   - Learning rate: `lr = 1/50`
5. **Simplify** the resulting rational
6. **Return** the new rational

### Encoding

Rationals are encoded as strings:
- Format: `"numerator/denominator"` (e.g., `"3/4"`, `"-5/2"`)
- Integers: `"7"` (equivalent to `"7/1"`)
- Initial starting point: `"1/1"`

## Installation

### Prerequisites

1. **Julia 1.8+**
   ```bash
   # macOS with Homebrew
   brew install julia

   # Or download from https://julialang.org/downloads/
   ```

2. **Python 3.10+**
   ```bash
   # Should already be installed
   python3 --version
   ```

3. **GP/Pari** (for conductor computation)
   ```bash
   # macOS with Homebrew
   brew install pari

   # Verify installation
   gp --version
   ```

4. **SageMath** (optional, for Mestre's family)
   ```bash
   # macOS with Homebrew
   brew install --cask sage
   ```

### Julia Packages

Install required Julia packages:

```bash
julia -e 'using Pkg; Pkg.add(["Dictionaries", "Plots", "Combinatorics", "StatsBase", "ArgParse"])'
```

### Python Packages

From the project directory:

```bash
pip3 install torch tokenizers numpy
```

## Running the Problem

### Option 1: Julia Local Search Only (Quick Test)

Run just the local search algorithm without transformer training:

```bash
cd /Users/ao6409/Desktop/projects/patternboost/assets/transformers_math_experiments

# Edit search_fc.jl to uncomment the elliptic problem:
# Line 17: include("problem_elliptic_simple.jl")

julia search_fc.jl ./output 120 1000 1000 500
```

**Parameters:**
- `./output`: Output directory
- `120`: Number of local searches per batch
- `1000`: Number of initial empty objects
- `1000`: Final database size
- `500`: Target database size

**Output:**
- `./output/search_output_1.txt`: Best rationals found
- `./output/plot_1.png`: Score distribution
- `./output/distribution.txt`: Score statistics

### Option 2: Full PatternBoost Loop (With Transformer)

Run the complete algorithm with transformer-based global learning:

```bash
python3 fc_loop.py \
  --dump_path ./output \
  --exp_name elliptic_test \
  --num_initial_empty_objects 1000 \
  --final_database_size 1000 \
  --target_db_size 5000 \
  --max_epochs 10 \
  --max-steps 5000 \
  --cpu true \
  --nb_threads 4
```

Here's the command for a toy run to check the pipeline works:

```bash
python3 fc_loop.py     --dump_path ./output     --exp_name elliptic_gpu_test     --num_initial_empty_objects 20     --nb_local_searches 8     --final_database_size 20     --target_db_size 50     --sample-only 20     --max_epochs 5     --max-steps 20     --nb_threads 8
```

**Key Parameters:**
- `--dump_path`: Output directory
- `--exp_name`: Experiment name
- `--num_initial_empty_objects`: Initial rollouts (first iteration)
- `--final_database_size`: Training set size
- `--max_epochs`: Number of PatternBoost iterations
- `--max-steps`: Training steps per iteration
- `--cpu`: Use CPU instead of GPU
- `--nb_threads`: Number of Julia threads

### Option 3: Test Script

Run the test script to verify everything works:

```bash
julia test_elliptic.jl
```

This will test:
- Rational number parsing
- Discriminant computation
- GP/Pari integration
- Gradient descent
- Reward function

## Example Results

From our testing:

| t     | Conductor | Reward   | Notes                    |
|-------|-----------|----------|--------------------------|
| 0/1   | 36        | -36.0    | Best so far              |
| 1/1   | 496       | -496.0   | Starting point           |
| 2/1   | 472       | -472.0   | Better than starting     |
| 1/2   | 14080     | -14080.0 | Much worse               |
| -71/25| ?         | ?        | After gradient descent   |

## Switching Between Problems

To switch between different problems, edit `search_fc.jl` (lines 14-17):

```julia
# Choose the problem to work on here!

#include("problem_triangle_free.jl")       # Triangle-free graphs
#include("problem_4_cycle_free.jl")        # 4-cycle-free graphs
#include("problem_permanent_avoid_123.jl") # Permanent with 312-avoidance
include("problem_elliptic_simple.jl")      # Elliptic curve conductors ← Active
```

Simply comment out the current problem and uncomment the one you want to run.

## Implementation Details

### File Structure

- `problem_elliptic_simple.jl`: Main problem definition
- `test_elliptic.jl`: Test suite
- `search_fc.jl`: Local search framework
- `fc_loop.py`: Full PatternBoost loop
- `constants.jl`: Global constants

### Key Functions

All defined in `problem_elliptic_simple.jl`:

1. **`empty_starting_point()`**: Returns `"1/1"`
2. **`greedy_search_from_startpoint(db, obj)`**: Gradient descent on discriminant
3. **`reward_calc(obj)`**: Computes `-conductor(t)` via GP/Pari
4. **`compute_conductor(num, den)`**: Calls GP/Pari to compute conductor
5. **`eval_discriminant(num, den)`**: Evaluates D(t)
6. **`eval_discriminant_derivative(num, den)`**: Evaluates D'(t)

### Future Work: Mestre's Family

For more advanced experiments, you can compute the discriminant of Mestre's rank ≥ 12 family:

```bash
sage compute_mestre_discriminant.sage
```

This will output the discriminant polynomial D(t) and its derivative D'(t) for Mestre's family, which can then be hardcoded into a new problem file.

## Troubleshooting

### GP/Pari Not Found
```bash
brew install pari
# Verify: gp --version
```

### Julia Package Errors
```bash
julia -e 'using Pkg; Pkg.update(); Pkg.resolve()'
```

### Threading Issues
If you see threading errors, try running with a single thread:
```bash
julia -t 1 search_fc.jl ./output 120 1000 1000 500
```

### Python Dependencies
```bash
pip3 install --upgrade torch tokenizers numpy
```

## References

- **PatternBoost**: Original framework for mathematical construction search
- **Mestre's Family**: Rank ≥ 12 elliptic curves over ℚ
- **GP/Pari**: Computer algebra system for number theory
- **Conductor**: Invariant measuring arithmetic complexity of elliptic curves

## Contact

For questions about this problem, refer to the main PatternBoost documentation.
