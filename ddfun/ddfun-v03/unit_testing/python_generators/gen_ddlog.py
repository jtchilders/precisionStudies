# gen_ddlog.py
import random
import struct
from decimal import Decimal, getcontext
import numpy as np

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


def generate_ddlog_test_case():
    """
    Generate a test case for ddlog.
    """
    # Generate a random positive high-precision value
    base_value = Decimal(random.uniform(0.1, 10))

    # Compute the natural logarithm value
    log_value = base_value.ln()

    # Convert to double-double format
    a_hi, a_lo = generate_double_double(base_value)
    result_hi, result_lo = generate_double_double(log_value)

    # Verify correctness of double-double representation
    valid_a, error_a = verify_double_double(base_value, a_hi, a_lo)
    valid_result, error_result = verify_double_double(log_value, result_hi, result_lo)

    if not (valid_a and valid_result):
        raise ValueError(
            f"Verification failed:\n"
            f"  a: {base_value} (hi={a_hi}, lo={a_lo}, error={error_a})\n"
            f"  result: {log_value} (hi={result_hi}, lo={result_lo}, error={error_result})"
        )

    return {
        "a": (a_hi, a_lo),
        "expected": (result_hi, result_lo),
    }


def write_test_cases_to_file(filename, num_cases=10):
    """
    Generate test cases for ddlog and write them to a file.
    Each line contains:
      hi_a lo_a expected_hi expected_lo
    """
    with open(filename, "w") as f:
        for _ in range(num_cases):
            test_case = generate_ddlog_test_case()
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
    Generate and save test cases for ddlog to a binary file.
    """
    test_cases = [generate_ddlog_test_case() for _ in range(num_cases)]
    write_test_cases_to_binary(filename, test_cases)
    print(f"Test cases successfully written to {filename}")


if __name__ == "__main__":
    # write_test_cases_to_file("ddlog_test_cases.txt", num_cases=10)
    generate_binary_test_file("ddlog_test_cases.bin", num_cases=10)