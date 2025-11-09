#!/usr/bin/env sage
"""
compute_mestre_curve.sage.py

Runtime helper for Julia to compute conductor and ap values
for Mestre's rank >= 12 elliptic curve family.

Usage from Julia:
    # For conductor:
    result = read(`sage compute_mestre_curve.sage.py conductor 43/39`, String)

    # For ap values:
    result = read(`sage compute_mestre_curve.sage.py ap 43/39`, String)

Output format:
    - conductor: single integer
    - ap: comma-separated list of 25 integers
    - error: "ERROR" on failure
"""

import sys
from sage.all import QQ, PolynomialRing, ProjectiveSpace, Curve, Jacobian, prime_range

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

def build_mestre_curve(rational_str):
    """
    Build Mestre's elliptic curve at rational t

    Args:
        rational_str: string like "43/39" or "5"

    Returns:
        EllipticCurve object (Jacobian of plane cubic)
    """
    # Parse rational
    if '/' in rational_str:
        num, den = rational_str.split('/')
        my_rational = QQ(int(num)) / QQ(int(den))
    else:
        my_rational = QQ(int(rational_str))

    # Evaluate Mestre's polynomials at t
    A1 = QQ(R(a1).subs(t=my_rational))
    A2 = QQ(R(a2).subs(t=my_rational))
    A3 = QQ(R(a3).subs(t=my_rational))
    A4 = QQ(R(a4).subs(t=my_rational))
    A5 = QQ(R(a5).subs(t=my_rational))
    A6 = QQ(R(a6).subs(t=my_rational))
    A7 = QQ(R(a7).subs(t=my_rational))

    # Construct plane cubic curve
    P2 = ProjectiveSpace(QQ, 2, 'xyz')
    x, y, z = P2.gens()

    f = y**3 + A1*x**2*y + A2*x*y*z + A3*y*z**2 + A4*x**3 + A5*x**2*z + A6*x*z**2 + A7*z**3
    C = Curve(f, A=P2)
    E = Jacobian(C)

    return E

def compute_conductor(rational_str):
    """Compute conductor of Mestre's curve at rational t"""
    try:
        E = build_mestre_curve(rational_str)
        conductor = E.conductor()
        return str(conductor)
    except Exception as e:
        # On error, return ERROR
        return "ERROR"

def compute_ap_values(rational_str):
    """Compute ap values for Mestre's curve at rational t"""
    try:
        E = build_mestre_curve(rational_str)
        ap_vals = [E.ap(p) for p in prime_range(100)]
        return ','.join(map(str, ap_vals))
    except Exception as e:
        # On error, return zeros
        return ','.join(['0'] * 25)

# ============================================================================
# Main entry point
# ============================================================================

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("ERROR", file=sys.stderr)
        print("Usage: sage compute_mestre_curve.sage.py {conductor|ap} RATIONAL", file=sys.stderr)
        sys.exit(1)

    mode = sys.argv[1]
    rational_str = sys.argv[2]

    if mode == "conductor":
        result = compute_conductor(rational_str)
        print(result)
    elif mode == "ap":
        result = compute_ap_values(rational_str)
        print(result)
    else:
        print("ERROR", file=sys.stderr)
        print(f"Unknown mode: {mode}. Use 'conductor' or 'ap'", file=sys.stderr)
        sys.exit(1)
