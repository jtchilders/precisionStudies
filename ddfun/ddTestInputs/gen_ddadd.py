import random
import struct
import numpy as np
from decimal import Decimal, getcontext

# Set high precision for calculations
getcontext().prec = 50
PI = Decimal("3.14159265358979323846264338327950288419716939937510")  # High-precision PI


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


def generate_ddadd_test_case():
    """
    Generate a test case for ddadd using PI as the base value.
    """
    # Generate two random scaling factors
    scale_a = Decimal(random.uniform(0.1, 10))  # Scale factor for first number
    scale_b = Decimal(random.uniform(0.1, 10))  # Scale factor for second number

    # Derive the high-precision values from PI
    a = PI * scale_a
    b = PI / scale_b  # Ensure variety in operations

    # Compute expected result
    result = a + b

    # Convert to double-double format
    a_hi, a_lo = generate_double_double(a)
    b_hi, b_lo = generate_double_double(b)
    result_hi, result_lo = generate_double_double(result)

    # Verify correctness of double-double representation
    valid_a, error_a = verify_double_double(a, a_hi, a_lo)
    valid_b, error_b = verify_double_double(b, b_hi, b_lo)
    valid_result, error_result = verify_double_double(result, result_hi, result_lo)

    if not (valid_a and valid_b and valid_result):
        raise ValueError(
            f"Verification failed:\n"
            f"  a: {a} (hi={a_hi}, lo={a_lo}, error={error_a})\n"
            f"  b: {b} (hi={b_hi}, lo={b_lo}, error={error_b})\n"
            f"  result: {result} (hi={result_hi}, lo={result_lo}, error={error_result})"
        )

    return {
        "dda": (a_hi, a_lo),
        "ddb": (b_hi, b_lo),
        "expected": (result_hi, result_lo),
    }


def generate_test_cases(num_cases):
    return [generate_ddadd_test_case() for _ in range(num_cases)]

def write_test_cases_to_text_file(filename, test_cases):
    """
    Generate test cases for ddadd and write them to a file.
    Each line contains:
      hi_a lo_a hi_b lo_b expected_hi expected_lo
    """
    with open(filename, "w") as f:
        for case in test_cases:
            f.write(
                f"{case['dda'][0]:.16e} {case['dda'][1]:.16e} "
                f"{case['ddb'][0]:.16e} {case['ddb'][1]:.16e} "
                f"{case['expected'][0]:.16e} {case['expected'][1]:.16e}\n"
            )
    print(f"Test cases successfully written to {filename}")

def write_test_cases_to_binary(filename, test_cases):
    """
    Save test cases to a binary file.
    Each case: [hi_a, lo_a, hi_b, lo_b, expected_hi, expected_lo].
    """
    with open(filename, "wb") as f:
        for case in test_cases:
            f.write(struct.pack(
                "6d",  # 6 double-precision floats
                case["dda"][0], case["dda"][1],
                case["ddb"][0], case["ddb"][1],
                case["expected"][0], case["expected"][1]
            ))
    print(f"Test cases successfully written to {filename}")

if __name__ == "__main__":
    test_cases = generate_test_cases(10)
    write_test_cases_to_text_file("data/ddadd_test_cases.txt", test_cases)
    write_test_cases_to_binary("data/ddadd_test_cases.bin", test_cases)
