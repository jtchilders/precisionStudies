# gen_ddneg.py
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

def generate_ddneg_test_case():
    """
    Generate a test case for ddneg.
    """
    # Generate a random scaling factor
    scale = Decimal(random.uniform(0.1, 10))  # Scale factor for the number

    # Derive a high-precision value
    value = Decimal(random.random()) * scale

    # Compute expected result
    negative_value = -value

    # Convert to double-double format
    value_hi, value_lo = generate_double_double(value)
    negative_value_hi, negative_value_lo = generate_double_double(negative_value)

    # Verify correctness of double-double representation
    valid_value, error_value = verify_double_double(value, value_hi, value_lo)
    valid_negative_value, error_negative_value = verify_double_double(negative_value, negative_value_hi, negative_value_lo)

    if not (valid_value and valid_negative_value):
        print(
            f"Verification failed:\n"
            f"  value: {value} (hi={value_hi}, lo={value_lo}, error={error_value})\n"
            f"  negative_value: {negative_value} (hi={negative_value_hi}, lo={negative_value_lo}, error={error_negative_value})"
        )

    return {
        "a": (value_hi, value_lo),
        "expected": (negative_value_hi, negative_value_lo),
    }

def generate_test_cases(num_cases):
    return [generate_ddneg_test_case() for _ in range(num_cases)]

def write_test_cases_to_text_file(filename, test_cases):
    """
    Generate test cases for ddneg and write them to a file.
    Each line contains:
      hi_a lo_a expected_hi expected_lo
    """
    with open(filename, "w") as f:
        for case in test_cases:
            f.write(
                f"{case['a'][0]:.16e} {case['a'][1]:.16e} "
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
                case["a"][0], case["a"][1],
                case["expected"][0], case["expected"][1]
            ))
    print(f"Test cases successfully written to {filename}")

if __name__ == "__main__":
    test_cases = generate_test_cases(10)
    write_test_cases_to_text_file("data/ddneg_test_cases.txt", test_cases)
    write_test_cases_to_binary("data/ddneg_test_cases.bin", test_cases)