#!/usr/bin/env sage
"""
Verify that ap values computed from:
1. Plane cubic (Mestre's original form)
2. Weierstrass form (converted via Jacobian)
are the same
"""

from sage.all import QQ, PolynomialRing, ProjectiveSpace, Curve, Jacobian, EllipticCurve, prime_range

# Test rational
test_t = QQ(43)/QQ(39)

print("="*80)
print(f"Testing at t = {test_t}")
print("="*80)
print()

# ============================================================================
# METHOD 1: Plane Cubic (the correct way from user's Python code)
# ============================================================================

R = PolynomialRing(QQ, 't')
t = R.gen()

# Mestre's coefficients from the paper
a1 = [-26940, 51220, -26940]
a2 = [-1320, 17280, 17280, -1320]
a3 = [-18876, -153828, 301221, -153828, -18776]
a4 = [-1489600, 1489600, 1489600, -1489600]
a5 = [5816880, 8043880, -27463500, 8043880, 5816880]
a6 = [3416160, -24166320, 19202040, 19202040, -24166320, 3416160]
a7 = [-745360, -15468024, 18853764, -138394, 18853764, -15468024, -745360]

# Evaluate at test_t
A1 = QQ(R(a1).subs(t=test_t))
A2 = QQ(R(a2).subs(t=test_t))
A3 = QQ(R(a3).subs(t=test_t))
A4 = QQ(R(a4).subs(t=test_t))
A5 = QQ(R(a5).subs(t=test_t))
A6 = QQ(R(a6).subs(t=test_t))
A7 = QQ(R(a7).subs(t=test_t))

P2 = ProjectiveSpace(QQ, 2, 'xyz')
x, y, z = P2.gens()

f = y**3 + A1*x**2*y + A2*x*y*z + A3*y*z**2 + A4*x**3 + A5*x**2*z + A6*x*z**2 + A7*z**3
C = Curve(f, A=P2)
E_plane = Jacobian(C)

print("METHOD 1: From plane cubic")
print(f"E = {E_plane}")
print()

ap_plane = [E_plane.ap(p) for p in prime_range(100)]
print(f"ap values (first 10): {ap_plane[:10]}")
print(f"Full ap values: {ap_plane}")
print()

# ============================================================================
# METHOD 2: Weierstrass form (using A4_NUMER/A6_NUMER coefficients)
# ============================================================================

# These are the polynomial coefficients from problem_mestre_rank12.jl
A4_NUMER_COEFFS = [
    -5013558100494028800,
    -31223157822576691200,
    36850819289195923200,
    375767910105823324800,
    -909015201933727021200,
    627735360031439126400,
    -190333976336927546400,
    627968231066603606400,
    -909272532073107781200,
    376003880650827484800,
    36714196278367683200,
    -31220058312737011200,
    -4996097133777388800,
]
A4_DENOM = 3

A6_NUMER_COEFFS = [
    22668186571459155780894720000,
    203387292236079580258713600000,
    134498583657781883343636480000,
    -3633649404971677060116049920000,
    -1672385445663582687556688640000,
    40212991310412954466845627840000,
    -66623930317405025667601698960000,
    -46203257144166699494281956480000,
    292205940790151219565054099600000,
    -429332179402122264136129783680000,
    292343874284328539800854200400000,
    -46390902011588484059487857280000,
    -66508882717375442192538796560000,
    40200621730127022808322677440000,
    -1693906103722165170601285440000,
    -3626101530671701430935032320000,
    135762415484271030280730880000,
    202896505624301923064140800000,
    22566782308794025025387520000,
]
A6_DENOM = 27

# Evaluate a4(t) and a6(t)
def eval_poly(coeffs, t_val):
    """Evaluate polynomial with coeffs in descending degree order"""
    result = QQ(0)
    for coeff in coeffs:
        result = result * t_val + QQ(coeff)
    return result

a4_val = eval_poly(A4_NUMER_COEFFS, test_t) / A4_DENOM
a6_val = eval_poly(A6_NUMER_COEFFS, test_t) / A6_DENOM

print("METHOD 2: From Weierstrass form (using Julia coefficients)")
print(f"a4({test_t}) = {a4_val}")
print(f"a6({test_t}) = {a6_val}")
print()

# Create elliptic curve in Weierstrass form: y^2 = x^3 + a4*x + a6
E_weierstrass = EllipticCurve([0, 0, 0, a4_val, a6_val])
print(f"E = {E_weierstrass}")
print()

ap_weierstrass = [E_weierstrass.ap(p) for p in prime_range(100)]
print(f"ap values (first 10): {ap_weierstrass[:10]}")
print(f"Full ap values: {ap_weierstrass}")
print()

# ============================================================================
# COMPARISON
# ============================================================================

print("="*80)
print("COMPARISON")
print("="*80)

if ap_plane == ap_weierstrass:
    print("✓ ap values MATCH! Both methods give the same result.")
    print()
    print("CONCLUSION: The Weierstrass coefficients in Julia ARE correct!")
    print("The bug must be elsewhere (perhaps in how ap values are computed at runtime)")
else:
    print("✗ ap values DO NOT MATCH!")
    print()
    print("Differences:")
    for i, p in enumerate(prime_range(100)):
        if ap_plane[i] != ap_weierstrass[i]:
            print(f"  p={p}: plane={ap_plane[i]}, weierstrass={ap_weierstrass[i]}")
    print()
    print("CONCLUSION: The Weierstrass coefficients in Julia are WRONG!")

print("="*80)

# Also check conductors
print()
print("Checking conductors:")
cond_plane = E_plane.conductor()
cond_weierstrass = E_weierstrass.conductor()
print(f"Conductor (plane cubic): {cond_plane}")
print(f"Conductor (Weierstrass): {cond_weierstrass}")
if cond_plane == cond_weierstrass:
    print("✓ Conductors match")
else:
    print("✗ Conductors DO NOT match!")
