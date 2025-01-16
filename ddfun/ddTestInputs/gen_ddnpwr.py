# gen_ddnpwr.py
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

def generate_ddnpwr_test_case():
    """
    Generate a test case for ddnpwr.
    """
    # Generate random base and exponent
    base = Decimal(random.uniform(0.1, 10))
    exponent = random.randint(-5, 5)
    
    # Compute expected result
    if base == 0 and exponent < 0:
        raise ValueError('Cannot compute negative power of zero')
    result = base ** exponent if base != 0 else Decimal(0)

    # Convert to double-double format
    base_hi, base_lo = generate_double_double(base)
    result_hi, result_lo = generate_double_double(result)

    # Verify correctness of double-double representation
    valid_base, error_base = verify_double_double(base, base_hi, base_lo)
    valid_result, error_result = verify_double_double(result, result_hi, result_lo)

    if not (valid_base and valid_result):
        print(
            f"Verification failed:\n"
            f"  base: {base} (hi={base_hi}, lo={base_lo}, error={error_base})\n"
            f"  result: {result} (hi={result_hi}, lo={result_lo}, error={error_result})"
        )

    return {
        "a": (base_hi, base_lo),
        "n": exponent,
        "expected": (result_hi, result_lo),
    }

def generate_test_cases(num_cases):
    return [generate_ddnpwr_test_case() for _ in range(num_cases)]

def write_test_cases_to_text_file(filename, test_cases):
    """
    Generate test cases for ddnpwr and write them to a file.
    Each line contains:
      hi_a lo_a n expected_hi expected_lo
    """
    with open(filename, "w") as f:
        for case in test_cases:
            f.write(
                f"{case['a'][0]:.16e} {case['a'][1]:.16e} {case['n']} "
                f"{case['expected'][0]:.16e} {case['expected'][1]:.16e}\n"
            )
    print(f"Test cases successfully written to {filename}")


def write_test_cases_to_binary(filename, test_cases):
    """
    Save test cases to a binary file.
    Each case: [hi_a, lo_a, n, expected_hi, expected_lo].
    """
    with open(filename, "wb") as f:
        for case in test_cases:
            f.write(struct.pack(
                "2di2d",  # 2 double-precision floats, 1 int, 2 double-precision floats
                case["a"][0], case["a"][1], case["n"], case["expected"][0], case["expected"][1]
            ))
    print(f"Test cases successfully written to {filename}")



if __name__ == "__main__":
    test_cases = generate_test_cases(10)
    write_test_cases_to_text_file("data/ddnpwr_test_cases.txt", test_cases)
    write_test_cases_to_binary("data/ddnpwr_test_cases.bin", test_cases)