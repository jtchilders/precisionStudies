# gen_ddsqrt.py
import random
import struct
import numpy as np
from decimal import Decimal, getcontext

# Set high precision for calculations
getcontext().prec = 50

def generate_double_double(base_value):
    """
    Generate a double-double representation (hi, lo) of a high-precision value.
    """
    hi = float(base_value)
    lo = float(base_value - Decimal(hi))  # Compute the remainder
    return hi, lo

def verify_double_double(original, hi, lo, tolerance=1e-30):
    """
    Verify that hi and lo accurately represent the original value.
    """
    reconstructed = Decimal(hi) + Decimal(lo)
    error = abs(original - reconstructed)
    return error < Decimal(tolerance), error

def generate_ddsqrt_test_case():
    """
    Generate a test case for ddsqrt.
    """
    # Generate a random factor for variety
    scale = Decimal(random.uniform(0.1, 100))

    # Derive the high-precision value
    a_value = scale**2  # Ensure a positive number to sqrt
    result_value = a_value.sqrt()

    # Convert to double-double format
    a_hi, a_lo = generate_double_double(a_value)
    result_hi, result_lo = generate_double_double(result_value)

    # Verify correctness of double-double representation
    valid_a, error_a = verify_double_double(a_value, a_hi, a_lo)
    valid_result, error_result = verify_double_double(result_value, result_hi, result_lo)

    if not (valid_a and valid_result):
        raise ValueError(
            f"Verification failed:\n"
            f"  a: {a_value} (hi={a_hi}, lo={a_lo}, error={error_a})\n"
            f"  result: {result_value} (hi={result_hi}, lo={result_lo}, error={error_result})"
        )

    return {
        "dda": (a_hi, a_lo),
        "expected": (result_hi, result_lo),
    }


def generate_test_cases(num_cases):
    return [generate_ddsqrt_test_case() for _ in range(num_cases)]

def write_test_cases_to_text_file(filename, test_cases):
    """
    Save test cases to a text file.
    Each case: [hi_a, lo_a, expected_hi, expected_lo].
    """
    with open(filename, "w") as f:
        for case in test_cases:
            f.write(
                f"{case['dda'][0]:.16e} {case['dda'][1]:.16e} "
                f"{case['expected'][0]:.16e} {case['expected'][1]:.16e}\n"
            )
    print(f"Test cases successfully written to {filename}")

def write_test_cases_to_binary(filename, test_cases):
    """
    Save test cases to a binary file.
    Each case: [hi_a, lo_a, expected_hi, expected_lo].
    """
    with open(filename, "wb") as f:
        for case in test_cases:
            f.write(struct.pack(
                "4d",  # 4 double-precision floats
                case["dda"][0], case["dda"][1],
                case["expected"][0], case["expected"][1]
            ))
    print(f"Test cases successfully written to {filename}")


if __name__ == "__main__":
    test_cases = generate_test_cases(10)
    write_test_cases_to_text_file("data/ddsqrt_test_cases.txt", test_cases)
    write_test_cases_to_binary("data/ddsqrt_test_cases.bin", test_cases)