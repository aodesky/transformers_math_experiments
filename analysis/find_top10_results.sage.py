"""Find the top 10 best results (smallest conductors) from the GPU run"""

from sage.all import QQ, PolynomialRing, ProjectiveSpace, Curve, Jacobian
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

def compute_curve_and_conductor(X):
    """
    Compute elliptic curve and conductor for rational X
    Returns (E_weierstrass, conductor) or (None, None) on error
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

        # Create plane cubic
        P2 = ProjectiveSpace(QQ, 2, 'xyz')
        x, y, z = P2.gens()

        f = y**3 + A1*x**2*y + A2*x*y*z + A3*y*z**2 + A4*x**3 + A5*x**2*z + A6*x*z**2 + A7*z**3
        C = Curve(f, A=P2)

        # Get Jacobian (this is already the elliptic curve in Weierstrass form)
        E = Jacobian(C)

        # Compute conductor
        conductor = E.conductor()

        return E, conductor
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
    print("\nComputing conductors for all rationals...")
    print("="*80)

    results = []

    for i, rational_str in enumerate(sorted(rationals), 1):
        try:
            X = parse_rational(rational_str)
            E, conductor = compute_curve_and_conductor(X)

            if E is not None and conductor is not None:
                results.append({
                    'rational': rational_str,
                    'value': X,
                    'curve': E,
                    'conductor': conductor
                })
                print(f"[{i}/{len(rationals)}] t = {rational_str}: conductor = {conductor}")
            else:
                print(f"[{i}/{len(rationals)}] t = {rational_str}: FAILED")
        except Exception as e:
            print(f"[{i}/{len(rationals)}] t = {rational_str}: ERROR - {e}")

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
        print(f"  Elliptic Curve:")
        print(f"    {r['curve']}")
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
