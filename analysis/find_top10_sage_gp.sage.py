"""Find the top 10 best results using Sage's built-in GP interface"""

from sage.all import QQ, PolynomialRing, gp
import os
import sys

# Mestre's plane cubic coefficients from the paper
R = PolynomialRing(QQ, 't')
t = R.gen()

a1 = [-26940, 51220, -26940]
a2 = [-1320, 17280, 17280, -1320]
a3 = [-18876, -153828, 301221, -153828, -18776]
a4 = [-1489600, 1489600, 1489600, -1489600]
a5 = [5816880, 8043880, -27463500, 8043880, 5816880]
a6 = [3416160, -24166320, 19202040, 19202040, -24166320, 3416160]
a7 = [-745360, -15468024, 18853764, -138394, 18853764, -15468024, -745360]

def compute_conductor_gp(X):
    """
    Compute conductor using Sage's GP interface
    Returns (conductor, elliptic_curve_gp) or (None, None) on error
    """
    try:
        # Evaluate Mestre's polynomials at t=X
        A1 = QQ(R(a1).subs(t=X))
        A2 = QQ(R(a2).subs(t=X))
        A3 = QQ(R(a3).subs(t=X))
        A4 = QQ(R(a4).subs(t=X))
        A5 = QQ(R(a5).subs(t=X))
        A6 = QQ(R(a6).subs(t=X))
        A7 = QQ(R(a7).subs(t=X))

        # Build the plane cubic in GP
        # f = y^3 + A1*x^2*y + A2*x*y + A3*y + A4*x^3 + A5*x^2 + A6*x + A7
        gp.eval("x='x; y='y")

        # Construct the polynomial string
        f_str = f"y^3 + ({A1})*x^2*y + ({A2})*x*y + ({A3})*y + ({A4})*x^3 + ({A5})*x^2 + ({A6})*x + ({A7})"

        # Evaluate in GP to get the polynomial
        f_gp = gp(f_str)

        # Convert to elliptic curve using ellfromeqn
        weierstrass = gp.ellfromeqn(f_gp)

        # Check if ellfromeqn succeeded
        if weierstrass == 0 or weierstrass is None:
            return None, None

        E_gp = gp.ellinit(weierstrass)

        # Check if valid
        if E_gp == 0 or E_gp is None:
            return None, None

        # Get conductor from ellglobalred
        # ellglobalred returns [conductor, [u,r,s,t], c, factorization, local_data]
        try:
            globalred = E_gp.ellglobalred()

            # Check if globalred returned a valid result
            if not globalred or len(globalred) < 1:
                return None, None

            conductor = int(globalred[0])

            # Sanity check
            if conductor <= 0:
                return None, None

            return conductor, E_gp
        except:
            return None, None

    except Exception as e:
        print(f"Error for t={X}: {e}", file=sys.stderr)
        return None, None

def parse_rational(s):
    """Parse rational string like '3/5' or '7' into QQ"""
    s = s.strip()
    if '/' in s:
        parts = s.split('/')
        return QQ(int(parts[0])) / QQ(int(parts[1]))
    else:
        return QQ(int(s))

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
    print("\nComputing conductors using Sage's GP interface...")
    print("="*80)

    results = []

    for i, rational_str in enumerate(sorted(rationals), 1):
        try:
            X = parse_rational(rational_str)
            conductor, E_gp = compute_conductor_gp(X)

            if conductor is not None:
                results.append({
                    'rational': rational_str,
                    'value': X,
                    'conductor': conductor,
                    'curve_gp': E_gp
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
        print(f"  Elliptic Curve (GP format):")
        print(f"    {r['curve_gp']}")
        print()

    # Save to file
    output_file = f"{run_dir}/top10_results_sage.txt"
    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("TOP 10 ELLIPTIC CURVES WITH SMALLEST CONDUCTORS\n")
        f.write("From Mestre's Rank >= 12 Family\n")
        f.write("="*80 + "\n\n")

        for i, r in enumerate(results[:10], 1):
            f.write(f"RANK {i}:\n")
            f.write(f"  Rational parameter: t = {r['rational']}\n")
            f.write(f"  Conductor: {r['conductor']}\n")
            f.write(f"  Elliptic Curve (GP format):\n")
            f.write(f"    {r['curve_gp']}\n")
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

        # Check if in original input
        input_file = "elliptic_input_small.txt"
        if os.path.exists(input_file):
            input_rationals = set()
            with open(input_file) as f:
                for line in f:
                    line = line.strip()
                    if line:
                        input_rationals.add(line.split(',')[0])

            if results[0]['rational'] in input_rationals:
                print("\n⚠️  WARNING: Best result was in the original input dataset")
            else:
                print("\n✓ Best result is a NEW DISCOVERY (not in original input)")
    else:
        print("="*80)
        print("ERROR: No valid results found!")
        print("="*80)

if __name__ == "__main__":
    main()
