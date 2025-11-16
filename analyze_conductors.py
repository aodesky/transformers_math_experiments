"""Analyze conductors across initial dataset and PatternBoost results"""

print("Script started, importing modules...")
from conductor_from_rational import conductor_from_rational
from sage.all import QQ
import os
print("Imports complete!")

def parse_rational(s):
    """Parse rational string like '3/5' or '7' into QQ"""
    s = s.strip()
    if '/' in s:
        parts = s.split('/')
        return QQ(int(parts[0])) / QQ(int(parts[1]))
    else:
        return QQ(int(s))

def extract_rationals_from_file(filename):
    """Extract rationals (first column) from a file"""
    rationals = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Extract rational (before first comma)
            rational_str = line.split(',')[0]
            rationals.append(rational_str)
    return rationals

def main():
    print("="*80)
    print("ANALYZING LOG CONDUCTORS: Initial Dataset vs PatternBoost Results")
    print("="*80)
    print()

    # Part 1: Initial dataset
    print("PART 1: Initial Dataset (elliptic_input_small.txt)")
    print("-" * 80)

    initial_file = "elliptic_input_small.txt"
    initial_rationals = extract_rationals_from_file(initial_file)

    print(f"Found {len(initial_rationals)} rationals in initial dataset")
    print()

    initial_results = []
    for i, rational_str in enumerate(initial_rationals, 1):
        try:
            X = parse_rational(rational_str)
            log_conductor = conductor_from_rational(X)
            initial_results.append({
                'rational': rational_str,
                'log_conductor': float(log_conductor)
            })
            print(f"  [{i}/{len(initial_rationals)}] t = {rational_str:>10s} : log(C) = {log_conductor:.6f}")
        except Exception as e:
            print(f"  [{i}/{len(initial_rationals)}] t = {rational_str:>10s} : ERROR - {e}")

    print()

    # Part 2: PatternBoost results
    print("PART 2: PatternBoost Search Outputs")
    print("-" * 80)

    output_dir = "./output/elliptic_gpu_test/le5oyswdvp"
    all_pb_rationals = []

    for gen in range(1, 20):  # Try up to generation 20
        filename = f"{output_dir}/search_output_{gen}.txt"
        if not os.path.exists(filename):
            break

        gen_rationals = extract_rationals_from_file(filename)
        print(f"Generation {gen}: {len(gen_rationals)} rationals")
        all_pb_rationals.extend(gen_rationals)

    # Get unique rationals
    unique_pb_rationals = list(set(all_pb_rationals))
    print(f"\nTotal unique rationals across all generations: {len(unique_pb_rationals)}")
    print()

    pb_results = []
    for i, rational_str in enumerate(sorted(unique_pb_rationals), 1):
        try:
            X = parse_rational(rational_str)
            log_conductor = conductor_from_rational(X)
            pb_results.append({
                'rational': rational_str,
                'log_conductor': float(log_conductor)
            })
            print(f"  [{i}/{len(unique_pb_rationals)}] t = {rational_str:>10s} : log(C) = {log_conductor:.6f}")
        except Exception as e:
            print(f"  [{i}/{len(unique_pb_rationals)}] t = {rational_str:>10s} : ERROR - {e}")

    print()
    print("="*80)
    print("SUMMARY")
    print("="*80)
    print()

    # Question 1: Smallest log conductor in initial dataset
    if initial_results:
        initial_results.sort(key=lambda x: x['log_conductor'])
        min_initial = initial_results[0]
        print(f"(1) SMALLEST LOG CONDUCTOR IN STARTING DATASET:")
        print(f"    Rational: t = {min_initial['rational']}")
        print(f"    Log(Conductor): {min_initial['log_conductor']:.6f}")
        print()

    # Question 2: Smallest log conductor after PatternBoost
    if pb_results:
        pb_results.sort(key=lambda x: x['log_conductor'])
        min_pb = pb_results[0]
        print(f"(2) SMALLEST LOG CONDUCTOR AFTER PATTERNBOOST:")
        print(f"    Rational: t = {min_pb['rational']}")
        print(f"    Log(Conductor): {min_pb['log_conductor']:.6f}")
        print()

    # Improvement
    if initial_results and pb_results:
        improvement = min_initial['log_conductor'] - min_pb['log_conductor']
        percent_improvement = (improvement / min_initial['log_conductor']) * 100
        print(f"IMPROVEMENT:")
        print(f"    Reduction in log(C): {improvement:.6f}")
        print(f"    Percentage improvement: {percent_improvement:.2f}%")
        print()

        # Show top 10 from PatternBoost
        print("TOP 10 RATIONALS FROM PATTERNBOOST (by log conductor):")
        print("-" * 80)
        for i, result in enumerate(pb_results[:10], 1):
            print(f"  {i:2d}. t = {result['rational']:>10s} : log(C) = {result['log_conductor']:.6f}")
        print()

    print("="*80)

if __name__ == "__main__":
    main()
