import random
import struct
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


def generate_dddivd_test_case():
    """
    Generate a test case for dddivd.
    """
    # Generate a random double-double value and a divisor
    scale = Decimal(random.uniform(1.0, 10.0))
    dd_value = PI * scale  # Generate a double-double value using PI as base
    divisor = Decimal(random.uniform(0.1, 10.0))  # Random divisor

    # Compute expected result
    result = dd_value / divisor

    # Convert to double-double format
    dd_hi, dd_lo = generate_double_double(dd_value)
    result_hi, result_lo = generate_double_double(result)

    # Verify correctness of double-double representation
    valid_dd, error_dd = verify_double_double(dd_value, dd_hi, dd_lo)
    valid_result, error_result = verify_double_double(result, result_hi, result_lo)

    if not (valid_dd and valid_result):
        raise ValueError(
            f"Verification failed:\n"
            f"  dd_value: {dd_value} (hi={dd_hi}, lo={dd_lo}, error={error_dd})\n"
            f"  divisor: {divisor}\n"
            f"  result: {result} (hi={result_hi}, lo={result_lo}, error={error_result})"
        )

    return {
        "dd": (dd_hi, dd_lo),
        "divisor": float(divisor),
        "expected": (result_hi, result_lo),
    }

def generate_test_cases(num_cases):
    return [generate_dddivd_test_case() for _ in range(num_cases)]


def write_test_cases_to_text_file(filename, test_cases):
    """
    Save test cases for dddivd to a text file.
    Each case: [dd_hi, dd_lo, divisor, expected_hi, expected_lo].
    Each line contains:
      hi_a lo_a hi_b lo_b expected_hi expected_lo
    """
    with open(filename, "w") as f:
        for case in test_cases:
            f.write(
                f"{case['dd'][0]:.16e} {case['dd'][1]:.16e} "
                f"{case['divisor']:.16e} "
                f"{case['expected'][0]:.16e} {case['expected'][1]:.16e}\n"
            )

    print(f"Test cases successfully written to {filename}")

def write_test_cases_to_binary(filename, test_cases):
    """
    Save test cases for dddivd to a binary file.
    Each case: [dd_hi, dd_lo, divisor, expected_hi, expected_lo].
    """
    with open(filename, "wb") as f:
        for case in test_cases:
            f.write(struct.pack(
                "5d",  # 5 double-precision floats
                case["dd"][0], case["dd"][1],
                case["divisor"],
                case["expected"][0], case["expected"][1]
            ))

    print(f"Test cases successfully written to {filename}")


if __name__ == "__main__":
    test_cases = generate_test_cases(10)
    write_test_cases_to_text_file("dddivd_test_cases.txt", test_cases)
    write_test_cases_to_binary("dddivd_test_cases.bin", test_cases)
