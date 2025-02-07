import random
import struct
import math
import sys
from mpmath import mp

# Set high precision for calculations
mp.dps = 50

def generate_double_double(base_value):
    """
    Generate a double-double representation (hi, lo) of a high-precision value.
    """
    hi = float(base_value)
    lo = float(base_value - mp.mpf(hi))  # Compute the remainder
    return hi, lo

def calculate_exponent_difference(a, b):
    """
    Calculate the difference of scale between two double values
    """
    a_abs = math.fabs(a)
    b_abs = math.fabs(b)
    try:
        a_exp = int(math.log10(a_abs))
    except ValueError:
        a_exp = 0
    try:
        b_exp = int(math.log10(b_abs))
    except ValueError:
        b_exp = 0
    out = math.fabs(a_exp - b_exp)
    return out

def verify_double_double(original, hi, lo, tolerance):
    """
    Verify that hi and lo accurately represent the original value.
    """
    reconstructed = mp.mpf(hi) + mp.mpf(lo)
    error = abs(original - reconstructed)
    if error == 0.0:
        return True, error
    scale_diff = calculate_exponent_difference(original, error)
    return scale_diff > tolerance, error

def generate_ddcsshr_test_case():
    """
    Generate a test case for ddcsshr using a random value.
    """
    # Generate a random value with a range suitable for hyperbolic functions
    a_value = mp.rand() * 10  # Adjust range as needed
    # Calculate hyperbolic sine and cosine using high precision
    sinh_value = mp.sinh(a_value)
    cosh_value = mp.cosh(a_value)

    # Convert to inputs and outputs to double-double format
    a_hi, a_lo = generate_double_double(a_value)
    x_hi, x_lo = generate_double_double(cosh_value)
    y_hi, y_lo = generate_double_double(sinh_value)

    # Verify correctness of double-double representation
    tolerance = 30
    valid_a, error_a = verify_double_double(a_value, a_hi, a_lo, tolerance)
    valid_x, error_x = verify_double_double(cosh_value, x_hi, x_lo, tolerance)
    valid_y, error_y = verify_double_double(sinh_value, y_hi, y_lo, tolerance)

    if not (valid_a and valid_x and valid_y):
        print(
            f"Verification of ddcsshr failed for tolerance {tolerance}:\n"
            f"  a: {a_value} (hi={a_hi}, lo={a_lo}, error={error_a})\n"
            f"  cosh: {cosh_value} (hi={x_hi}, lo={x_lo}, error={error_x})\n"
            f"  sinh: {sinh_value} (hi={y_hi}, lo={y_lo}, error={error_y})"
        )

    return {
        "a": (a_hi, a_lo),
        "cosh_expected": (x_hi, x_lo),
        "sinh_expected": (y_hi, y_lo),
    }

def generate_test_cases(num_cases):
    return [generate_ddcsshr_test_case() for _ in range(num_cases)]

def write_test_cases_to_text_file(filename, test_cases):
    """
    Generate test cases for ddcsshr and write them to a file.
    Each line contains:
    hi_a lo_a cosh_expected_hi cosh_expected_lo sinh_expected_hi sinh_expected_lo
    """
    with open(filename, "w") as f:
        for case in test_cases:
            f.write(
                f"{case['a'][0]:.16e} {case['a'][1]:.16e} "
                f"{case['cosh_expected'][0]:.16e} {case['cosh_expected'][1]:.16e} "
                f"{case['sinh_expected'][0]:.16e} {case['sinh_expected'][1]:.16e}\n"
            )
    print(f"Test cases successfully written to {filename}")

def write_test_cases_to_binary(filename, test_cases):
    """
    Save test cases to a binary file.
    Each case: [hi_a, lo_a, cosh_expected_hi, cosh_expected_lo, sinh_expected_hi, sinh_expected_lo].
    """
    with open(filename, "wb") as f:
        for case in test_cases:
            f.write(struct.pack(
                "6d",  # 6 double-precision floats
                case["a"][0], case["a"][1],
                case["cosh_expected"][0], case["cosh_expected"][1],
                case["sinh_expected"][0], case["sinh_expected"][1]
            ))
    print(f"Test cases successfully written to {filename}")

if __name__ == "__main__":
    print("Generating test cases for ddcsshr...")
    ntests = 10
    if len(sys.argv) > 1:
        ntests = int(sys.argv[1])
    test_cases = generate_test_cases(ntests)
    write_test_cases_to_text_file("data/ddcsshr_test_cases.txt", test_cases)
    write_test_cases_to_binary("data/ddcsshr_test_cases.bin", test_cases)