#!/usr/bin/env python3
"""
Convert Mestre dataset to the format expected by problem_elliptic_simple.jl
Input: sign,numerator,denominator,conductor
Output: numerator/denominator (one per line)
"""

input_file = "/Users/ao6409/Desktop/projects/patternboost/assets/mestre_rank12_10_30_100.txt"
output_file = "/Users/ao6409/Desktop/projects/patternboost/assets/transformers_math_experiments/elliptic_input.txt"

with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
    # Skip header
    next(f_in)

    for line in f_in:
        parts = line.strip().split(',')
        if len(parts) == 4:
            sign = int(parts[0])
            numerator = int(parts[1])
            denominator = int(parts[2])

            # Apply sign to numerator
            if sign == -1:
                numerator = -numerator

            # Write in format "num/den"
            f_out.write(f"{numerator}/{denominator}\n")

print(f"Converted {input_file}")
print(f"Output written to {output_file}")
