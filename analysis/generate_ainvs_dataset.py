"""Generate ainvs dataset from rationals

Reads an input file containing rationals and computes
the ainvs (a-invariants) for each rational using Mestre's family.

Output format: rational,a1,a2,a3,a4,a6

Usage:
    sage -python generate_ainvs_dataset.py <input_file>

Output will be written to ainvs_of_<input_filename> in the same directory.
"""

import sys
import os
sys.path.append('..')

from conductor_from_rational import ainvs_from_rational
from sage.all import QQ

def parse_rational(s):
    """Parse rational string like '3/5' or '7' into QQ"""
    s = s.strip()
    if '/' in s:
        parts = s.split('/')
        return QQ(int(parts[0])) / QQ(int(parts[1]))
    else:
        return QQ(int(s))

def main():
    if len(sys.argv) < 2:
        print("Usage: sage -python generate_ainvs_dataset.py <input_file>")
        print("Example: sage -python generate_ainvs_dataset.py ../output/elliptic_gpu_test/le5oyswdvp/search_output_1.txt")
        sys.exit(1)

    input_file = sys.argv[1]

    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)

    # Generate output filename: prepend "ainvs_of_" to the input filename
    input_dir = os.path.dirname(input_file)
    input_basename = os.path.basename(input_file)
    output_basename = f"ainvs_of_{input_basename}"
    output_file = os.path.join(input_dir, output_basename)

    print(f"Reading from {input_file}")
    print(f"Writing to {output_file}")
    print()

    with open(input_file, 'r') as f_in:
        lines = f_in.readlines()

    total = len(lines)
    print(f"Processing {total} rationals...")
    print()

    results = []

    for i, line in enumerate(lines, 1):
        line = line.strip()
        if not line:
            continue

        # Extract rational (first column before comma)
        rational_str = line.split(',')[0]

        try:
            X = parse_rational(rational_str)
            ainvs = ainvs_from_rational(X)

            # Format: rational,a1,a2,a3,a4,a6
            ainvs_str = ','.join(str(a) for a in ainvs)
            result_line = f"{rational_str},{ainvs_str}"
            results.append(result_line)

            if i % 100 == 0:
                print(f"  Processed {i}/{total} rationals...")

        except Exception as e:
            print(f"  ERROR on line {i}: {rational_str} - {e}")

    print()
    print(f"Successfully processed {len(results)}/{total} rationals")
    print(f"Writing to {output_file}...")

    with open(output_file, 'w') as f_out:
        for line in results:
            f_out.write(line + '\n')

    print("Done!")

if __name__ == "__main__":
    main()
