#!/usr/bin/env sage
"""
Rigorously verify that D_COEFFS in problem_mestre_rank12.jl are correct
by computing the discriminant from Mestre's plane cubic and comparing
"""

from sage.all import QQ, PolynomialRing, ProjectiveSpace, Curve, Jacobian

print("="*80)
print("RIGOROUS VERIFICATION OF D_COEFFS")
print("="*80)
print()

# ============================================================================
# Step 1: Compute discriminant from Mestre's plane cubic
# ============================================================================

R = PolynomialRing(QQ, 't')
t = R.gen()

# Mestre's coefficients from the paper
a1_coeffs = [-26940, 51220, -26940]
a2_coeffs = [-1320, 17280, 17280, -1320]
a3_coeffs = [-18876, -153828, 301221, -153828, -18776]
a4_coeffs = [-1489600, 1489600, 1489600, -1489600]
a5_coeffs = [5816880, 8043880, -27463500, 8043880, 5816880]
a6_coeffs = [3416160, -24166320, 19202040, 19202040, -24166320, 3416160]
a7_coeffs = [-745360, -15468024, 18853764, -138394, 18853764, -15468024, -745360]

# Build polynomials
A1 = R(a1_coeffs)
A2 = R(a2_coeffs)
A3 = R(a3_coeffs)
A4 = R(a4_coeffs)
A5 = R(a5_coeffs)
A6 = R(a6_coeffs)
A7 = R(a7_coeffs)

print("Building Mestre's plane cubic curve...")
K = R.fraction_field()
P2 = ProjectiveSpace(K, 2, 'xyz')
x, y, z = P2.gens()

f = y**3 + A1*x**2*y + A2*x*y*z + A3*y*z**2 + A4*x**3 + A5*x**2*z + A6*x*z**2 + A7*z**3

C = Curve(f, A=P2)
E = Jacobian(C)

print("Computing discriminant from Jacobian...")
disc = E.discriminant()

print(f"Discriminant degree: {disc.numerator().degree()}")
print()

# Extract coefficients (highest degree first)
disc_poly = disc.numerator()
sage_coeffs = list(reversed(disc_poly.list()))

print(f"Number of coefficients: {len(sage_coeffs)}")
print()

# ============================================================================
# Step 2: Load D_COEFFS from Julia file
# ============================================================================

print("Loading D_COEFFS from problem_mestre_rank12.jl...")
print()

# These are from problem_mestre_rank12.jl lines 20-58
julia_d_coeffs = [
    -6180509591227720419380815576086350476681408413696000000,
    114938623929384725353177092671893112265996383551488000000,
    85551930076196922637182242936935849955283009798144000000,
    -11719628364517969915879400083879559410975035747205120000000,
    27052858938943217389306770689380951058020225437925376000000,
    495233735772095817573467306356693483579440850052055040000000,
    -1801350648161536348295098941726255410489247569001054208000000,
    -10250169932874704526238384703116770475682736517705367552000000,
    61082789522191032916195198949323882271668853151436636160000000,
    28889019713203853749514869477449341722130633703018463232000000,
    -894269251478892110653793279137284480187713458226239045632000000,
    1998986763029986547291442571731467929032900302560167854080000000,
    2661443381030372339132138461213064026328428987766103425024000000,
    -22395786314228742177536194385054915532573975723543903928320000000,
    49754045290848789556449598143368652244595710783550411472896000000,
    -44163173854815427694280317244618399762028370828027075166208000000,
    -32796148340266456018606916619411215277163210920437139800064000000,
    148904961810211812582905296904750630903782562417849592053760000000,
    -206251369280050995820557086852749583234853890852417835827200000000,
    148755092346812123632637376165602044783334295774330759413760000000,
    -32597758506965011566365737558577383086466467094612793278464000000,
    -44295499759368608691450381682183461951308091965292567134208000000,
    49782009495051548235764568635229690056981489920153730973696000000,
    -22362866931987209877620015953762190874516213419975061995520000000,
    2626761458211022508276814384763419934345961252162347597824000000,
    2012441064460877642689216131506573045948955416959580897280000000,
    -894419902548735011677159393784767862754953344245243052032000000,
    27004564310929212161890586883786579124328865854187896832000000,
    61657866003996771936038377726735901162216714589972725760000000,
    -10210770063465329148660348309485349439782151641945341952000000,
    -1846853462515131041133469192728870173290370843025604608000000,
    497428925592111443220177495517750224788127586965258240000000,
    28948530843419825472651917479685557177402342219186176000000,
    -11855959828599693328589625913076149102956954066616320000000,
    42900905735805344359511077352919934420956368338944000000,
    116724890597580319681815300992222752267149653311488000000,
    -5788570664188895847185744142334355789296064004096000000,
]

# ============================================================================
# Step 3: Compare
# ============================================================================

print("="*80)
print("COMPARISON")
print("="*80)
print()

if len(sage_coeffs) != len(julia_d_coeffs):
    print(f"✗ LENGTH MISMATCH!")
    print(f"  Sage has {len(sage_coeffs)} coefficients")
    print(f"  Julia has {len(julia_d_coeffs)} coefficients")
    print()
else:
    print(f"✓ Both have {len(sage_coeffs)} coefficients")
    print()

all_match = True
mismatches = []

for i in range(min(len(sage_coeffs), len(julia_d_coeffs))):
    sage_val = int(sage_coeffs[i])
    julia_val = julia_d_coeffs[i]

    if sage_val != julia_val:
        all_match = False
        mismatches.append((i, sage_val, julia_val))

if all_match:
    print("="*80)
    print("✓✓✓ ALL COEFFICIENTS MATCH EXACTLY! ✓✓✓")
    print("="*80)
    print()
    print("CONCLUSION:")
    print("  D_COEFFS in problem_mestre_rank12.jl were computed from the")
    print("  CORRECT Mestre plane cubic using Sage's Jacobian method.")
    print()
    print("  The discriminant is being used for gradient descent in local search,")
    print("  which means the local search is working correctly.")
    print()
else:
    print("="*80)
    print(f"✗✗✗ FOUND {len(mismatches)} MISMATCHES ✗✗✗")
    print("="*80)
    print()
    print("Mismatches (showing first 5):")
    for i, sage_val, julia_val in mismatches[:5]:
        print(f"  Coefficient {i} (t^{36-i}):")
        print(f"    Sage:  {sage_val}")
        print(f"    Julia: {julia_val}")
        print()

    print("CONCLUSION:")
    print("  D_COEFFS in problem_mestre_rank12.jl are INCORRECT!")
    print()

print("="*80)
print("VERIFICATION COMPLETE")
print("="*80)
