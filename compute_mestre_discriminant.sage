"""compute_mestre_discriminant.sage

Computes the discriminant D(t) and its derivative D'(t) as polynomials
for Mestre's family of rank >= 12 elliptic curves.

Run with: sage compute_mestre_discriminant.sage
"""

from sage.all import QQ, PolynomialRing, ProjectiveSpace, Curve, Jacobian

# Polynomial ring in t
R = PolynomialRing(QQ, 't')
t = R.gen()

# Mestre's coefficients (from the paper)
a1_coeffs = [-26940, 51220, -26940]
a2_coeffs = [-1320, 17280, 17280, -1320]
a3_coeffs = [-18876, -153828, 301221, -153828, -18776]
a4_coeffs = [-1489600, 1489600, 1489600, -1489600]
a5_coeffs = [5816880, 8043880, -27463500, 8043880, 5816880]
a6_coeffs = [3416160, -24166320, 19202040, 19202040, -24166320, 3416160]
a7_coeffs = [-745360, -15468024, 18853764, -138394, 18853764, -15468024, -745360]

# Build polynomials in t
A1 = R(a1_coeffs)
A2 = R(a2_coeffs)
A3 = R(a3_coeffs)
A4 = R(a4_coeffs)
A5 = R(a5_coeffs)
A6 = R(a6_coeffs)
A7 = R(a7_coeffs)

print("Mestre's polynomial coefficients A1(t) through A7(t):")
print(f"A1(t) = {A1}")
print(f"A2(t) = {A2}")
print(f"A3(t) = {A3}")
print(f"A4(t) = {A4}")
print(f"A5(t) = {A5}")
print(f"A6(t) = {A6}")
print(f"A7(t) = {A7}")
print()

# Define the projective space and curve
# We'll work over the function field QQ(t)
K = R.fraction_field()  # QQ(t)
P2 = ProjectiveSpace(K, 2, 'xyz')
x, y, z = P2.gens()

# The plane cubic curve
f = y^3 + A1*x^2*y + A2*x*y*z + A3*y*z^2 + A4*x^3 + A5*x^2*z + A6*x*z^2 + A7*z^3

print("Defining the plane cubic curve...")
print(f"f = {f}")
print()

try:
    C = Curve(f, A=P2)
    print("Converting to Weierstrass form via Jacobian...")
    E = Jacobian(C)

    print(f"Elliptic curve E: {E}")
    print()

    # Get the Weierstrass coefficients
    ainvs = E.ainvs()
    print("Weierstrass coefficients [a1, a2, a3, a4, a6]:")
    for i, coeff in enumerate(['a1', 'a2', 'a3', 'a4', 'a6']):
        print(f"{coeff} = {ainvs[i]}")
    print()

    # Compute the discriminant
    disc = E.discriminant()
    print("="*80)
    print("DISCRIMINANT D(t):")
    print("="*80)
    print(disc)
    print()

    # Compute the derivative
    disc_derivative = disc.derivative(t)
    print("="*80)
    print("DERIVATIVE D'(t):")
    print("="*80)
    print(disc_derivative)
    print()

    # Print in a format easier to copy
    print("="*80)
    print("FOR JULIA/GP-PARI (numerator and denominator):")
    print("="*80)
    disc_numer = disc.numerator()
    disc_denom = disc.denominator()
    print(f"D_numerator(t) = {disc_numer}")
    print(f"D_denominator(t) = {disc_denom}")
    print()

    deriv_numer = disc_derivative.numerator()
    deriv_denom = disc_derivative.denominator()
    print(f"D'_numerator(t) = {deriv_numer}")
    print(f"D'_denominator(t) = {deriv_denom}")

except Exception as e:
    print(f"Error: {e}")
    print()
    print("Alternative: Computing at a specific value...")
    # Fallback: compute at t=1 to verify
    X_test = QQ(1)
    A1_val = QQ(R(a1_coeffs).subs(t=X_test))
    A2_val = QQ(R(a2_coeffs).subs(t=X_test))
    A3_val = QQ(R(a3_coeffs).subs(t=X_test))
    A4_val = QQ(R(a4_coeffs).subs(t=X_test))
    A5_val = QQ(R(a5_coeffs).subs(t=X_test))
    A6_val = QQ(R(a6_coeffs).subs(t=X_test))
    A7_val = QQ(R(a7_coeffs).subs(t=X_test))

    P2_val = ProjectiveSpace(QQ, 2, 'xyz')
    x, y, z = P2_val.gens()
    f_val = y^3 + A1_val*x^2*y + A2_val*x*y*z + A3_val*y*z^2 + A4_val*x^3 + A5_val*x^2*z + A6_val*x*z^2 + A7_val*z^3
    C_val = Curve(f_val, A=P2_val)
    E_val = Jacobian(C_val)

    print(f"At t={X_test}: E = {E_val}")
    print(f"Discriminant at t={X_test}: {E_val.discriminant()}")
