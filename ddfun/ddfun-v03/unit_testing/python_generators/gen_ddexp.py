# gen_ddexp.py
import random
import struct
import numpy as np
from decimal import Decimal, getcontext

# Set high precision for calculations
getcontext().prec = 50

def generate_double_double(value):
    """
    Generate a double-double representation (hi, lo) of a high-precision value.
    """
    hi = float(value)
    lo = float(value - Decimal(hi))  # Compute the remainder
    return hi, lo

def verify_double_double(original, hi, lo, tolerance=1e-30):
    """
    Verify that hi and lo accurately represent the original value.
    """
    reconstructed = Decimal(hi) + Decimal(lo)
    error = abs(original - reconstructed)
    return error < Decimal(tolerance), error

def generate_ddexp_test_case():
    """
    Generate a test case for ddexp.
    """
    # Generate a random value in range [-1, 1] for exponentiation
    scale = Decimal(random.uniform(-1.0, 1.0))
    a = scale

    # Compute expected result
    result = a.exp()  # Decimal exponential

    # Convert to double-double format
    a_hi, a_lo = generate_double_double(a)
    result_hi, result_lo = generate_double_double(result)

    # Verify correctness of double-double representation
    valid_a, error_a = verify_double_double(a, a_hi, a_lo)
    valid_result, error_result = verify_double_double(result, result_hi, result_lo)

    if not (valid_a and valid_result):
        raise ValueError(
            f"Verification failed:\n"
            f"  a: {a} (hi={a_hi}, lo={a_lo}, error={error_a})\n"
            f"  result: {result} (hi={result_hi}, lo={result_lo}, error={error_result})"
        )

    return {
        "a": (a_hi, a_lo),
        "expected": (result_hi, result_lo),
    }

def write_test_cases_to_file(filename, num_cases=10):
    """
    Generate test cases for ddexp and write them to a file.
    Each line contains:
      hi_a lo_a expected_hi expected_lo
    """
    with open(filename, "w") as f:
        for _ in range(num_cases):
            test_case = generate_ddexp_test_case()
            f.write(
                f"{test_case['a'][0]:.16e} {test_case['a'][1]:.16e} "
                f"{test_case['expected'][0]:.16e} {test_case['expected'][1]:.16e}\n"
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
                case["a"][0], case["a"][1],
                case["expected"][0], case["expected"][1]
            ))

def generate_binary_test_file(filename, num_cases=10):
    """
    Generate and save test cases for ddexp to a binary file.
    """
    test_cases = [generate_ddexp_test_case() for _ in range(num_cases)]
    write_test_cases_to_binary(filename, test_cases)
    print(f"Test cases successfully written to {filename}")

if __name__ == "__main__":
    # write_test_cases_to_file("ddexp_test_cases.txt", num_cases=10)
    generate_binary_test_file("ddexp_test_cases.bin", num_cases=10)