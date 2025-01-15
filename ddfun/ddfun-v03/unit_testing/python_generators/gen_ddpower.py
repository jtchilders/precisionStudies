# gen_ddpower.py
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

def compute_power(a, b):
    """
    Compute a^b using high-precision decimal arithmetic.
    """
    if a <= 0:
        raise ValueError("DDPOWER: A <= 0")
    log_a = Decimal(a).ln()
    power_result = (log_a * Decimal(b)).exp()
    return power_result

def generate_ddpower_test_case():
    """
    Generate a test case for ddpower using random positive values.
    """
    base_a = Decimal(random.uniform(0.1, 10))
    exp_b = Decimal(random.uniform(0.1, 5))

    try:
        result = compute_power(base_a, exp_b)
    except ValueError as e:
        raise ValueError(f"Invalid test case: {str(e)}")

    a_hi, a_lo = generate_double_double(base_a)
    b_hi, b_lo = generate_double_double(exp_b)
    result_hi, result_lo = generate_double_double(result)

    valid_a, error_a = verify_double_double(base_a, a_hi, a_lo)
    valid_b, error_b = verify_double_double(exp_b, b_hi, b_lo)
    valid_result, error_result = verify_double_double(result, result_hi, result_lo)

    if not (valid_a and valid_b and valid_result):
        print(
            f"Verification failed:\n"
            f"  a: {base_a} (hi={a_hi}, lo={a_lo}, error={error_a})\n"
            f"  b: {exp_b} (hi={b_hi}, lo={b_lo}, error={error_b})\n"
            f"  result: {result} (hi={result_hi}, lo={result_lo}, error={error_result})"
        )

    return {
        "a": (a_hi, a_lo),
        "b": (b_hi, b_lo),
        "expected": (result_hi, result_lo),
    }

def generate_test_cases(num_cases):
    cases = []
    while len(cases) < num_cases:
        try:
            cases.append(generate_ddpower_test_case())
        except ValueError:
            continue
    return cases

def write_test_cases_to_text_file(filename, test_case):
    """
    Generate test cases for ddpower and write them to a file.
    Each line contains:
      hi_a lo_a hi_b lo_b expected_hi expected_lo
    """
    with open(filename, "w") as f:
        for case in test_cases:
            f.write(
                f"{case['a'][0]:.16e} {case['a'][1]:.16e} "
                f"{case['b'][0]:.16e} {case['b'][1]:.16e} "
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
                case["a"][0], case["a"][1], 
                case["b"][0], case["b"][1],
                case["expected"][0], case["expected"][1]
            ))
    print(f"Test cases successfully written to {filename}")

if __name__ == "__main__":
    test_cases = generate_test_cases(10)
    write_test_cases_to_text_file("ddpower_test_cases.txt", test_cases)
    write_test_cases_to_binary("ddpower_test_cases.bin", test_cases)