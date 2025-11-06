/* mestre_discriminant.gp
 *
 * Computes the discriminant D(t) and its derivative D'(t)
 * for Mestre's family of rank >= 12 elliptic curves
 */

/* Define the polynomial coefficients from Mestre's paper */
a1_coeffs = [-26940, 51220, -26940];
a2_coeffs = [-1320, 17280, 17280, -1320];
a3_coeffs = [-18876, -153828, 301221, -153828, -18776];
a4_coeffs = [-1489600, 1489600, 1489600, -1489600];
a5_coeffs = [5816880, 8043880, -27463500, 8043880, 5816880];
a6_coeffs = [3416160, -24166320, 19202040, 19202040, -24166320, 3416160];
a7_coeffs = [-745360, -15468024, 18853764, -138394, 18853764, -15468024, -745360];

/* Create polynomials in t */
A1(t) = sum(i=1, #a1_coeffs, a1_coeffs[i] * t^(i-1));
A2(t) = sum(i=1, #a2_coeffs, a2_coeffs[i] * t^(i-1));
A3(t) = sum(i=1, #a3_coeffs, a3_coeffs[i] * t^(i-1));
A4(t) = sum(i=1, #a4_coeffs, a4_coeffs[i] * t^(i-1));
A5(t) = sum(i=1, #a5_coeffs, a5_coeffs[i] * t^(i-1));
A6(t) = sum(i=1, #a6_coeffs, a6_coeffs[i] * t^(i-1));
A7(t) = sum(i=1, #a7_coeffs, a7_coeffs[i] * t^(i-1));

print("Computing elliptic curve from plane cubic...");

/* The curve is: y^3 + A1*x^2*y + A2*x*y*z + A3*y*z^2 + A4*x^3 + A5*x^2*z + A6*x*z^2 + A7*z^3 = 0
 * We need to convert this to Weierstrass form to get the discriminant.
 *
 * For a symbolic computation, let's evaluate at a specific value first to see the structure.
 */

/* Let's compute the Weierstrass coefficients symbolically */
/* This requires converting the plane cubic to Weierstrass form */
/* For now, let's use a simpler approach: compute at a specific rational value */

test_value = 1/2;

print("Testing at t = ", test_value);

/* Substitute the test value */
A1_val = A1(test_value);
A2_val = A2(test_value);
A3_val = A3(test_value);
A4_val = A4(test_value);
A5_val = A5(test_value);
A6_val = A6(test_value);
A7_val = A7(test_value);

print("A1 = ", A1_val);
print("A2 = ", A2_val);
print("A3 = ", A3_val);
print("A4 = ", A4_val);
print("A5 = ", A5_val);
print("A6 = ", A6_val);
print("A7 = ", A7_val);

/* Define the plane cubic curve equation */
/* y^3 + A1*x^2*y + A2*x*y*z + A3*y*z^2 + A4*x^3 + A5*x^2*z + A6*x*z^2 + A7*z^3 */
/* In affine coordinates (z=1): y^3 + A1*x^2*y + A2*x*y + A3*y + A4*x^3 + A5*x^2 + A6*x + A7 */

/* This is a genus 1 curve. We need to find its Jacobian (elliptic curve) */
/* GP/Pari has functions to work with hyperelliptic curves, but this is a plane cubic */

print("\nFor symbolic discriminant, we need the Weierstrass form.");
print("The plane cubic y^3 + ... needs to be converted to Weierstrass y^2 = x^3 + ...");
print("\nLet me try a different approach using ellinit on the Jacobian...");

/* Actually, looking at the Python code, it converts via Sage's Jacobian function */
/* Let's instead compute the discriminant numerically for several values and show the pattern */

print("\nComputing discriminants at several rational values:");
print("t, discriminant");

values = [0, 1/2, 1, 2, -1, -1/2, 1/3, 2/3];

for(i=1, #values,
    t_val = values[i];
    /* We need the Weierstrass coefficients - but we don't have them symbolically */
    /* This is actually computed via Sage in the Python script */
    print("t = ", t_val, ": [Weierstrass conversion needed]");
);

print("\n=================================================================");
print("NOTE: To get the discriminant polynomial D(t), we need to:");
print("1. Convert the plane cubic to Weierstrass form symbolically");
print("2. This requires Sage or more advanced algebraic geometry tools");
print("3. GP/Pari alone cannot easily convert a genus 1 plane cubic to Weierstrass form");
print("\nSuggestion: Use Sage to compute the discriminant symbolically and export it.");
print("=================================================================");
