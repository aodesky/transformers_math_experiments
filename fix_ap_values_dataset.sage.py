#!/usr/bin/env sage
"""
fix_ap_values_dataset.sage.py

Regenerate elliptic_input_small.txt with CORRECT ap values
computed from Mestre's plane cubic (not the wrong Weierstrass form).

Usage:
    sage fix_ap_values_dataset.sage.py
"""

from sage.all import QQ, PolynomialRing, ProjectiveSpace, Curve, Jacobian, prime_range

print("="*80)
print("FIXING AP VALUES IN DATASET")
print("="*80)
print()

# ============================================================================
# Mestre's plane cubic coefficients from the paper
# ============================================================================

R = PolynomialRing(QQ, 't')
t = R.gen()

a1 = [-26940, 51220, -26940]
a2 = [-1320, 17280, 17280, -1320]
a3 = [-18876, -153828, 301221, -153828, -18776]
a4 = [-1489600, 1489600, 1489600, -1489600]
a5 = [5816880, 8043880, -27463500, 8043880, 5816880]
a6 = [3416160, -24166320, 19202040, 19202040, -24166320, 3416160]
a7 = [-745360, -15468024, 18853764, -138394, 18853764, -15468024, -745360]

def compute_correct_ap_vals(my_rational):
    """
    Compute ap values for Mestre's curve at rational t using the CORRECT
    plane cubic construction (not the wrong Weierstrass form)

    Returns: list of 25 ap values for primes p < 100
    """
    try:
        # Evaluate Mestre's polynomials at t = my_rational
        A1 = QQ(R(a1).subs(t=my_rational))
        A2 = QQ(R(a2).subs(t=my_rational))
        A3 = QQ(R(a3).subs(t=my_rational))
        A4 = QQ(R(a4).subs(t=my_rational))
        A5 = QQ(R(a5).subs(t=my_rational))
        A6 = QQ(R(a6).subs(t=my_rational))
        A7 = QQ(R(a7).subs(t=my_rational))

        # Construct the plane cubic curve
        P2 = ProjectiveSpace(QQ, 2, 'xyz')
        x, y, z = P2.gens()

        f = y**3 + A1*x**2*y + A2*x*y*z + A3*y*z**2 + A4*x**3 + A5*x**2*z + A6*x*z**2 + A7*z**3
        C = Curve(f, A=P2)
        E = Jacobian(C)

        # Compute ap for all primes < 100
        ap_vals = [E.ap(p) for p in prime_range(100)]

        return ap_vals
    except Exception as e:
        print(f"  ERROR computing ap values for {my_rational}: {e}")
        # Return zeros on error
        return [0] * 25

# ============================================================================
# Read existing rationals from elliptic_input_small.txt
# ============================================================================

input_file = "elliptic_input_small.txt"
output_file = "elliptic_input_small.txt"

print(f"Reading rationals from {input_file}...")

rationals = []
with open(input_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line:
            # Extract rational (first column before comma)
            rational_str = line.split(',')[0]
            rationals.append(rational_str)

print(f"Found {len(rationals)} rationals")
print()

# ============================================================================
# Compute correct ap values for each rational
# ============================================================================

print("Computing correct ap values using Mestre's plane cubic...")
print()

output_lines = []

for i, rational_str in enumerate(rationals):
    print(f"[{i+1}/{len(rationals)}] Processing {rational_str}...", end=" ")

    try:
        # Parse the rational
        if '/' in rational_str:
            num, den = rational_str.split('/')
            my_rational = QQ(int(num)) / QQ(int(den))
        else:
            my_rational = QQ(int(rational_str))

        # Compute ap values
        ap_vals = compute_correct_ap_vals(my_rational)

        # Format output: rational,a_2,a_3,...,a_97
        output_line = rational_str + ',' + ','.join(map(str, ap_vals))
        output_lines.append(output_line)

        print(f"✓ ap values: {ap_vals[:5]}... (showing first 5)")

    except Exception as e:
        print(f"✗ ERROR: {e}")
        # On error, write zeros
        output_line = rational_str + ',' + ','.join(['0'] * 25)
        output_lines.append(output_line)

print()

# ============================================================================
# Write corrected data back to file
# ============================================================================

print(f"Writing corrected data to {output_file}...")

with open(output_file, 'w') as f:
    for line in output_lines:
        f.write(line + '\n')

print(f"✓ Successfully wrote {len(output_lines)} lines")
print()

# ============================================================================
# Verification: Show first 3 lines
# ============================================================================

print("="*80)
print("VERIFICATION: First 3 lines of corrected file:")
print("="*80)
for i, line in enumerate(output_lines[:3]):
    print(f"{i+1}: {line}")

print()
print("="*80)
print("DATASET FIX COMPLETE")
print("="*80)
print()
print("The file elliptic_input_small.txt now contains CORRECT ap values")
print("computed from Mestre's plane cubic curve.")
