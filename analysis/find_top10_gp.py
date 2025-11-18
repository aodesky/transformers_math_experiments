#!/usr/bin/env python3
"""Find the top 10 best results (smallest conductors) using GP/PARI"""

import subprocess
import os
import sys
from fractions import Fraction
import tempfile

# Mestre's plane cubic coefficients from the paper
a1 = [-26940, 51220, -26940]
a2 = [-1320, 17280, 17280, -1320]
a3 = [-18876, -153828, 301221, -153828, -18776]
a4 = [-1489600, 1489600, 1489600, -1489600]
a5 = [5816880, 8043880, -27463500, 8043880, 5816880]
a6 = [3416160, -24166320, 19202040, 19202040, -24166320, 3416160]
a7 = [-745360, -15468024, 18853764, -138394, 18853764, -15468024, -745360]

def eval_poly(coeffs, t_frac):
    """
    Evaluate polynomial at rational t
    Coefficients are in increasing degree order (constant, t, t^2, ...)
    Returns Fraction
    """
    result = Fraction(0)
    t_power = Fraction(1)
    for c in coeffs:
        result += c * t_power
        t_power *= t_frac
    return result

def compute_conductor_gp(t_frac):
    """
    Compute conductor using GP/PARI
    Returns (conductor, elliptic_curve_string) or (None, None) on error
    """
    try:
        # Evaluate Mestre's polynomials at t
        A1 = eval_poly(a1, t_frac)
        A2 = eval_poly(a2, t_frac)
        A3 = eval_poly(a3, t_frac)
        A4 = eval_poly(a4, t_frac)
        A5 = eval_poly(a5, t_frac)
        A6 = eval_poly(a6, t_frac)
        A7 = eval_poly(a7, t_frac)

        # Create GP/PARI script
        gp_script = f"""
A1={A1.numerator}/{A1.denominator};
A2={A2.numerator}/{A2.denominator};
A3={A3.numerator}/{A3.denominator};
A4={A4.numerator}/{A4.denominator};
A5={A5.numerator}/{A5.denominator};
A6={A6.numerator}/{A6.denominator};
A7={A7.numerator}/{A7.denominator};

\\\\ Affine plane cubic (z=1)
x='x; y='y;
f = y^3 + A1*x^2*y + A2*x*y + A3*y + A4*x^3 + A5*x^2 + A6*x + A7;

\\\\ Convert to Weierstrass form
E = ellinit(ellfromeqn(f));

if (E == 0,
    print("INVALID"),
    print(ellglobalred(E)[1]);
    print(E)
);
quit();
"""

        # Write to temp file and execute
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gp', delete=False) as f:
            temp_file = f.name
            f.write(gp_script)

        try:
            result = subprocess.run(
                ['gp', '-q', '-s', '200000000', temp_file],
                capture_output=True,
                text=True,
                timeout=5
            )

            os.unlink(temp_file)

            lines = result.stdout.strip().split('\n')

            if len(lines) < 2 or lines[0] == "INVALID":
                return None, None

            conductor = int(lines[0])
            elliptic_curve = lines[1]

            return conductor, elliptic_curve

        except (subprocess.TimeoutExpired, ValueError) as e:
            if os.path.exists(temp_file):
                os.unlink(temp_file)
            return None, None

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return None, None

def parse_rational(s):
    """Parse rational string like '3/5' or '7' into Fraction"""
    s = s.strip()
    if '/' in s:
        parts = s.split('/')
        return Fraction(int(parts[0]), int(parts[1]))
    else:
        return Fraction(int(s))

def main():
    # Read all search output files
    run_dir = "./output/elliptic_gpu_test/le5oyswdvp"

    rationals = set()

    print("="*80)
    print("Extracting rationals from search output files...")
    print("="*80)

    for i in range(1, 10):  # Try up to 10 files
        filename = f"{run_dir}/search_output_{i}.txt"
        if not os.path.exists(filename):
            break

        print(f"Reading {filename}")
        with open(filename) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                # Extract rational (before first comma)
                rational_str = line.split(',')[0]
                rationals.add(rational_str)

    print(f"\nFound {len(rationals)} unique rationals")
    print("="*80)
    print("\nComputing conductors using GP/PARI (fast!)...")
    print("="*80)

    results = []

    for i, rational_str in enumerate(sorted(rationals), 1):
        try:
            t_frac = parse_rational(rational_str)
            conductor, curve_str = compute_conductor_gp(t_frac)

            if conductor is not None:
                results.append({
                    'rational': rational_str,
                    'value': t_frac,
                    'conductor': conductor,
                    'curve': curve_str
                })
                print(f"[{i}/{len(rationals)}] t = {rational_str:15s}: conductor = {conductor}")
            else:
                print(f"[{i}/{len(rationals)}] t = {rational_str:15s}: FAILED")

        except Exception as e:
            print(f"[{i}/{len(rationals)}] t = {rational_str:15s}: ERROR - {e}")

    print("\n" + "="*80)
    print(f"Successfully computed {len(results)} conductors")
    print("="*80)

    # Sort by conductor (smallest first)
    results.sort(key=lambda r: r['conductor'])

    # Display top 10
    print("\n" + "="*80)
    print("TOP 10 RESULTS (SMALLEST CONDUCTORS)")
    print("="*80)
    print()

    for i, r in enumerate(results[:10], 1):
        print(f"RANK {i}:")
        print(f"  Rational t = {r['rational']}")
        print(f"  Conductor = {r['conductor']}")
        print(f"  Elliptic Curve: {r['curve']}")
        print()

    # Save to file
    output_file = f"{run_dir}/top10_results.txt"
    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("TOP 10 ELLIPTIC CURVES WITH SMALLEST CONDUCTORS\n")
        f.write("From Mestre's Rank >= 12 Family\n")
        f.write("="*80 + "\n\n")

        for i, r in enumerate(results[:10], 1):
            f.write(f"RANK {i}:\n")
            f.write(f"  Rational parameter: t = {r['rational']}\n")
            f.write(f"  Conductor: {r['conductor']}\n")
            f.write(f"  Elliptic Curve (Weierstrass form):\n")
            f.write(f"    {r['curve']}\n")
            f.write("\n")

        f.write("="*80 + "\n")
        f.write(f"Total curves analyzed: {len(results)}\n")
        f.write("="*80 + "\n")

    print(f"Results saved to {output_file}")
    print()

    if len(results) > 0:
        print("="*80)
        print(f"BEST RESULT: t = {results[0]['rational']} with conductor = {results[0]['conductor']}")
        print("="*80)
    else:
        print("="*80)
        print("ERROR: No valid results found!")
        print("="*80)

if __name__ == "__main__":
    main()
