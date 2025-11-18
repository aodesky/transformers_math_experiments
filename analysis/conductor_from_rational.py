"""curve_from_rational.py

    Generate the elliptic curve from a given rational number using Mestre's family.

"""

from sage.all import QQ, PolynomialRing, ProjectiveSpace, Curve, Jacobian, log, RR, ZZ
from sage.libs.pari.all import pari
from cysignals.signals import AlarmInterrupt

class pari_alarm:
    def __init__(self, seconds: int): self.seconds = int(seconds)
    def __enter__(self): pari(f"alarm({self.seconds})")
    def __exit__(self, exc_type, exc, tb): pari("alarm(0)")

R = PolynomialRing(QQ, 't')
t = R.gen()

CONDUCTOR_TIMEOUT = 3  # seconds  
BIG_LOG_CONDUCTOR = 1000 # large value to represent timeout

# The following are from Mestre's paper
a1 = [-26940, 51220, -26940]
a2 = [-1320, 17280, 17280, -1320]
a3 = [-18876, -153828, 301221, -153828, -18776]
a4 = [-1489600, 1489600, 1489600, -1489600]
a5 = [5816880, 8043880, -27463500, 8043880, 5816880]
a6 = [3416160, -24166320, 19202040, 19202040, -24166320, 3416160]
a7 = [-745360, -15468024, 18853764, -138394, 18853764, -15468024, -745360]

def ainvs_from_rational(X):
    # X is the rational

    A1 = QQ(R(a1).subs(t=X))
    A2 = QQ(R(a2).subs(t=X))
    A3 = QQ(R(a3).subs(t=X))
    A4 = QQ(R(a4).subs(t=X))
    A5 = QQ(R(a5).subs(t=X))
    A6 = QQ(R(a6).subs(t=X))
    A7 = QQ(R(a7).subs(t=X))

    P2 = ProjectiveSpace(QQ, 2, 'xyz')
    x, y, z = P2.gens()

    f = y**3 + A1*x**2*y + A2*x*y*z + A3*y*z**2 + A4*x**3 + A5*x**2*z + A6*x*z**2 + A7*z**3
    C = Curve(f, A=P2)
    E = Jacobian(C)

    ainvs = list(E.ainvs())  # [a1,a2,a3,a4,a6]
    return ainvs

def conductor_from_rational(X):
    # X is the rational

    ainvs = ainvs_from_rational(X)

    try:
        with pari_alarm(CONDUCTOR_TIMEOUT):          # 10 second timeout
            Epari = pari.ellinit(ainvs)  # cypari2 call, no gp/pexpect involved
            red = pari.ellglobalred(Epari)
            conductor = ZZ(red[0])       # first entry is the conductor
            return RR(log(conductor))
    except AlarmInterrupt as e:
        print(f"Timeout on X={X}, skipping.")
        return RR(BIG_LOG_CONDUCTOR)


if __name__ == "__main__":
    # Example usage
    X = QQ('3/4')  # Example rational
    result = conductor_from_rational(X)
    print(f"Logarithm of the conductor for X={X} is approximately {result}")