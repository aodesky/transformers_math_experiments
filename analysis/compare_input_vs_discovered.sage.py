"""Compare original input conductors vs PatternBoost discoveries using Sage's GP interface"""

from sage.all import QQ, PolynomialRing, gp

R = PolynomialRing(QQ, 't')
t = R.gen()

a1 = [-26940, 51220, -26940]
a2 = [-1320, 17280, 17280, -1320]
a3 = [-18876, -153828, 301221, -153828, -18776]
a4 = [-1489600, 1489600, 1489600, -1489600]
a5 = [5816880, 8043880, -27463500, 8043880, 5816880]
a6 = [3416160, -24166320, 19202040, 19202040, -24166320, 3416160]
a7 = [-745360, -15468024, 18853764, -138394, 18853764, -15468024, -745360]

def compute_conductor(X):
    try:
        A1 = QQ(R(a1).subs(t=X))
        A2 = QQ(R(a2).subs(t=X))
        A3 = QQ(R(a3).subs(t=X))
        A4 = QQ(R(a4).subs(t=X))
        A5 = QQ(R(a5).subs(t=X))
        A6 = QQ(R(a6).subs(t=X))
        A7 = QQ(R(a7).subs(t=X))

        gp.eval("x='x; y='y")
        f_str = f"y^3 + ({A1})*x^2*y + ({A2})*x*y + ({A3})*y + ({A4})*x^3 + ({A5})*x^2 + ({A6})*x + ({A7})"
        f_gp = gp(f_str)

        E_gp = gp.ellinit(gp.ellfromeqn(f_gp))

        if E_gp == 0:
            return None

        globalred = E_gp.ellglobalred()
        conductor = int(globalred[0])

        if conductor <= 0:
            return None

        return conductor
    except:
        return None

def parse_rational(s):
    s = s.strip()
    if '/' in s:
        parts = s.split('/')
        return QQ(int(parts[0])) / QQ(int(parts[1]))
    else:
        return QQ(int(s))

# Read original input
with open('elliptic_input_small.txt') as f:
    input_rationals = [line.split(',')[0].strip() for line in f if line.strip()]

print("="*80)
print("Computing conductors for ORIGINAL INPUT dataset (20 rationals)")
print("="*80)

input_results = []
for i, rat_str in enumerate(input_rationals, 1):
    X = parse_rational(rat_str)
    cond = compute_conductor(X)
    if cond:
        input_results.append((rat_str, cond))
        print(f"[{i}/20] {rat_str:15s}: {cond}")
    else:
        print(f"[{i}/20] {rat_str:15s}: FAILED")

input_results.sort(key=lambda x: x[1])

print("\n" + "="*80)
print("TOP 10 FROM ORIGINAL INPUT:")
print("="*80)
for i, (rat, cond) in enumerate(input_results[:10], 1):
    print(f"{i:2d}. t = {rat:15s}  conductor ≈ {float(cond):.2e}")

print("\n" + "="*80)
print("TOP 10 FROM PATTERNBOOST:")
print("="*80)

patternboost_top10 = [
    ('-8/3', 87654633297206983217828903738001686524560421993783479273052752),
    ('-8/13', 199011678383820637349691153572208063247253133346181204589574546365136),
    ('22/19', 28873664962283626861455048557077576859860055239061828030470253297285974),
    ('-6/5', 27800305059667362679889223064726232359183717932870998802322221026855782050),
    ('-18/13', 41282644642135529662285433707605061209837390189115500464176495473722416094),
    ('-8/7', 369528603095978369520466225357857900989540273066404261228054052481262535696),
    ('-3/5', 1169372431006312281631802334811043808099912354039071205706602727335792163200),
    ('-9/5', 25761127938525800403051978904777468842086739804051939340219755632708802089600),
    ('-11/4', 55320557873630627619556262821967929416185771671954371538127540545858261828400),
    ('-8/5', 57012211695948370834414875811979238022169903627680283178334608286990675745200),
]

for i, (rat, cond) in enumerate(patternboost_top10, 1):
    print(f"{i:2d}. t = {rat:15s}  conductor ≈ {float(cond):.2e}")

print("\n" + "="*80)
print("CRITICAL COMPARISON:")
print("="*80)

if input_results:
    best_input_rat, best_input_cond = input_results[0]
    best_pb_rat, best_pb_cond = patternboost_top10[0]

    print(f"\nBest from ORIGINAL INPUT:")
    print(f"  t = {best_input_rat}")
    print(f"  Conductor = {best_input_cond}")
    print(f"  Log10(conductor) ≈ {float(best_input_cond).bit_length() * 0.301:.1f}")

    print(f"\nBest from PATTERNBOOST:")
    print(f"  t = {best_pb_rat}")
    print(f"  Conductor = {best_pb_cond}")
    print(f"  Log10(conductor) ≈ {float(best_pb_cond).bit_length() * 0.301:.1f}")

    print("\n" + "="*80)
    if best_pb_cond < best_input_cond:
        improvement = best_input_cond / best_pb_cond
        print(f"✓✓✓ SUCCESS! PatternBoost found a {improvement:.2e}x BETTER conductor!")
        print(f"✓✓✓ The ML approach WORKED - we discovered better elliptic curves!")
    else:
        ratio = best_pb_cond / best_input_cond
        print(f"✗✗✗ FAILURE: PatternBoost result is {ratio:.2e}x WORSE")
        print(f"✗✗✗ The ML did NOT improve on the original dataset")
    print("="*80)
