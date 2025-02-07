# gen_ddcssnr.py
import random
import struct
import math
import sys
from mpmath import mp

# Set high precision for calculations
mp.dps = 50

def generate_double_double(base_value):
    hi = float(base_value)
    lo = float(base_value - mp.mpf(hi))
    return hi, lo

def calculate_exponent_difference(a, b):
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
    reconstructed = mp.mpf(hi) + mp.mpf(lo)
    error = abs(original - reconstructed)
    if error == 0.0:
        return True, error
    scale_diff = calculate_exponent_difference(original, error)
    return scale_diff > tolerance, error

def generate_ddcssnr_test_case():
    # Generate a random angle within [-π, π)
    a = mp.rand() * 2 * mp.pi - mp.pi

    # Compute expected cos and sin results using mpmath
    cos_result = mp.cos(a)
    sin_result = mp.sin(a)

    # Convert inputs and outputs to double-double format
    a_hi, a_lo = generate_double_double(a)
    cos_hi, cos_lo = generate_double_double(cos_result)
    sin_hi, sin_lo = generate_double_double(sin_result)

    # Verify correctness of double-double representation
    tolerance = 30
    valid_a, error_a = verify_double_double(a, a_hi, a_lo, tolerance)
    valid_cos, error_cos = verify_double_double(cos_result, cos_hi, cos_lo, tolerance)
    valid_sin, error_sin = verify_double_double(sin_result, sin_hi, sin_lo, tolerance)

    if not (valid_a and valid_cos and valid_sin):
        print(
            f"Verification of ddcssnr failed for tolerance {tolerance}:\n"
            f"  a: {a} (hi={a_hi}, lo={a_lo}, error={error_a})\n"
            f"  cos_result: {cos_result} (hi={cos_hi}, lo={cos_lo}, error={error_cos})\n"
            f"  sin_result: {sin_result} (hi={sin_hi}, lo={sin_lo}, error={error_sin})"
        )

    return {
        "a": (a_hi, a_lo),
        "expected_cos": (cos_hi, cos_lo),
        "expected_sin": (sin_hi, sin_lo),
    }

def generate_test_cases(num_cases):
    return [generate_ddcssnr_test_case() for _ in range(num_cases)]

def write_test_cases_to_text_file(filename, test_cases):
    with open(filename, "w") as f:
        for case in test_cases:
            f.write(
                f"{case['a'][0]:.16e} {case['a'][1]:.16e} "
                f"{case['expected_cos'][0]:.16e} {case['expected_cos'][1]:.16e} "
                f"{case['expected_sin'][0]:.16e} {case['expected_sin'][1]:.16e}\n"
            )
    print(f"Test cases successfully written to {filename}")

def write_test_cases_to_binary(filename, test_cases):
    with open(filename, "wb") as f:
        for case in test_cases:
            f.write(struct.pack(
                "6d",  # 6 double-precision floats
                case["a"][0], case["a"][1],
                case["expected_cos"][0], case["expected_cos"][1],
                case["expected_sin"][0], case["expected_sin"][1]
            ))
    print(f"Test cases successfully written to {filename}")

if __name__ == "__main__":
    print("Generating test cases for ddcssnr...")
    ntests = 10
    if len(sys.argv) > 1:
        ntests = int(sys.argv[1])
    test_cases = generate_test_cases(ntests)
    write_test_cases_to_text_file("data/ddcssnr_test_cases.txt", test_cases)
    write_test_cases_to_binary("data/ddcssnr_test_cases.bin", test_cases)