"""
Compute log(conductors) in parallel using Sage with timeout handling

Reads data/ainvs_of_big_starter_dataset.txt and computes
the log of the conductor for each elliptic curve using parallel processing.

Input format: rational,a1,a2,a3,a4,a6
Output format: rational,log_conductor

Uses timeout mechanism to handle curves with very large conductors.
"""

from sage.all import QQ, ZZ, RR, log
from sage.libs.pari.all import pari
from cysignals.signals import AlarmInterrupt
import multiprocessing as mp
from functools import partial
import sys

# Global configuration parameters
CONDUCTOR_TIMEOUT = 1  # Timeout in seconds - adjust as needed
BIG_LOG_CONDUCTOR = 10000  # Default value for timeout cases
MAX_CURVES = 0  # Max curves to process (0 = process all)
NUM_PROCESSES = 0  # Number of parallel processes (0 = use all cores)

class pari_alarm:
    def __init__(self, seconds: int):
        self.seconds = int(seconds)
    def __enter__(self):
        pari(f"alarm({self.seconds})")
    def __exit__(self, exc_type, exc, tb):
        pari("alarm(0)")


def compute_conductor_single(ainvs_str, timeout_seconds=CONDUCTOR_TIMEOUT, big_value=BIG_LOG_CONDUCTOR):
    """
    Compute log(conductor) for a single curve from ainvs string.

    Args:
        ainvs_str: String in format "rational,a1,a2,a3,a4,a6"
        timeout_seconds: Timeout in seconds
        big_value: Default value to return on timeout

    Returns:
        Tuple: (rational_str, log_conductor)
    """
    parts = ainvs_str.strip().split(',')
    rational_str = parts[0]
    ainvs = [QQ(parts[i]) for i in range(1, 6)]

    try:
        with pari_alarm(timeout_seconds):
            Epari = pari.ellinit(ainvs)
            red = pari.ellglobalred(Epari)
            conductor = ZZ(red[0])  # First entry is the conductor
            log_conductor = RR(log(conductor))
    except AlarmInterrupt:
        print(f"Timeout on {rational_str}, using default value {big_value}")
        log_conductor = RR(big_value)
    except Exception as e:
        print(f"Error on {rational_str}: {e}, using default value {big_value}")
        log_conductor = RR(big_value)

    return (rational_str, log_conductor)


def main():
    print("Reading ../data/ainvs_of_big_starter_dataset.txt...")

    # Read the entire file
    with open("../data/ainvs_of_big_starter_dataset.txt", 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    total_curves = len(lines)

    # Limit number of curves if MAX_CURVES > 0
    if MAX_CURVES > 0 and MAX_CURVES < total_curves:
        n = MAX_CURVES
        lines = lines[:n]
        print(f"Found {total_curves} curves, processing first {n} for testing")
    else:
        n = total_curves
        print(f"Found {n} curves, processing all")

    print()

    # Determine number of processes
    if NUM_PROCESSES == 0:
        num_processes = mp.cpu_count()
    else:
        num_processes = NUM_PROCESSES

    print(f"Parallelization enabled: using {num_processes} processes")
    print(f"Timeout per curve: {CONDUCTOR_TIMEOUT} seconds")
    print()

    print("Computing log(conductors) in parallel...")
    print()

    # Create a pool of workers and process in parallel
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(compute_conductor_single, lines)

    print()
    print("Writing results to ../data/one_sec_log_conductors_of_big_starter_dataset.txt...")

    # Write results to file
    with open("../data/one_sec_log_conductors_of_big_starter_dataset.txt", 'w') as f:
        for rational_str, log_conductor in results:
            f.write(f"{rational_str},{log_conductor}\n")

    print()
    print(f"Done! Processed {n} curves.")

    # Optional: show some statistics
    log_conductors = [result[1] for result in results]
    print()
    print("Statistics:")
    print(f"  Minimum log(conductor): {min(log_conductors)}")
    print(f"  Maximum log(conductor): {max(log_conductors)}")


if __name__ == "__main__":
    main()
