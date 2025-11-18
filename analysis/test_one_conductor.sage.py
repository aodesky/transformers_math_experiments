"""Test conductor computation on a single rational"""

from conductor_from_rational import conductor_from_rational
from sage.all import QQ
import time

# Test with first rational from elliptic_input_small.txt
rational_str = "43/39"
X = QQ(rational_str)

print(f"Testing conductor computation for t = {rational_str}")
print(f"Timeout is set to 60 seconds in conductor_from_rational.py")
print(f"Starting computation...")
start_time = time.time()

try:
    log_conductor = conductor_from_rational(X)
    elapsed = time.time() - start_time
    print(f"\nSUCCESS!")
    print(f"  Rational: t = {rational_str}")
    print(f"  Log(Conductor): {log_conductor}")
    print(f"  Time taken: {elapsed:.2f} seconds")
except Exception as e:
    elapsed = time.time() - start_time
    print(f"\nERROR!")
    print(f"  Rational: t = {rational_str}")
    print(f"  Error: {e}")
    print(f"  Time taken: {elapsed:.2f} seconds")
