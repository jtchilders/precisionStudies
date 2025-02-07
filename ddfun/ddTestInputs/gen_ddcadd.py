# gen_ddcadd.py
import random
import struct
import math
import sys
from mpmath import mp, mpc

# Set high precision for calculations
mp.dps = 50

def generate_double_double_complex(base_real, base_imag):
    """
    Generate a double-double representation (hi, lo) of a high-precision complex value.
    """
    hi_real = float(base_real)
    lo_real = float(base_real - mp.mpf(hi_real))
    hi_imag = float(base_imag)
    lo_imag = float(base_imag - mp.mpf(hi_imag))
    return (hi_real, lo_real), (hi_imag, lo_imag)

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

def generate_ddcadd_test_case():
    """
    Generate a test case for ddcadd using random complex numbers as the base value.
    """
    # Generate random complex numbers in high precision
    a_real, a_imag = mp.rand(), mp.rand()
    b_real, b_imag = mp.rand(), mp.rand()

    # Compute expected result
    a = mpc(real=a_real, imag=a_imag)
    b = mpc(real=b_real, imag=b_imag)
    result = a + b

    # Convert to inputs and outputs in double-double format
    (a_real_hi, a_real_lo), (a_imag_hi, a_imag_lo) = generate_double_double_complex(a_real, a_imag)
    (b_real_hi, b_real_lo), (b_imag_hi, b_imag_lo) = generate_double_double_complex(b_real, b_imag)
    (result_real_hi, result_real_lo), (result_imag_hi, result_imag_lo) = generate_double_double_complex(result.real, result.imag)

    # Verify correctness of double-double representation
    tolerance = 30
    valid_a_real, error_a_real = verify_double_double(a_real, a_real_hi, a_real_lo, tolerance)
    valid_a_imag, error_a_imag = verify_double_double(a_imag, a_imag_hi, a_imag_lo, tolerance)
    valid_b_real, error_b_real = verify_double_double(b_real, b_real_hi, b_real_lo, tolerance)
    valid_b_imag, error_b_imag = verify_double_double(b_imag, b_imag_hi, b_imag_lo, tolerance)
    valid_result_real, error_result_real = verify_double_double(result.real, result_real_hi, result_real_lo, tolerance)
    valid_result_imag, error_result_imag = verify_double_double(result.imag, result_imag_hi, result_imag_lo, tolerance)

    if not (valid_a_real and valid_a_imag and valid_b_real and valid_b_imag and valid_result_real and valid_result_imag):
        print(
            f"Verification of ddcadd failed for tolerance {tolerance}:\n"
            f"  a_real: {a_real} (hi={a_real_hi}, lo={a_real_lo}, error={error_a_real})\n"
            f"  a_imag: {a_imag} (hi={a_imag_hi}, lo={a_imag_lo}, error={error_a_imag})\n"
            f"  b_real: {b_real} (hi={b_real_hi}, lo={b_real_lo}, error={error_b_real})\n"
            f"  b_imag: {b_imag} (hi={b_imag_hi}, lo={b_imag_lo}, error={error_b_imag})\n"
            f"  result_real: {result.real} (hi={result_real_hi}, lo={result_real_lo}, error={error_result_real})\n"
            f"  result_imag: {result.imag} (hi={result_imag_hi}, lo={result_imag_lo}, error={error_result_imag})"
        )

    return {
        "dda": (a_real_hi, a_real_lo, a_imag_hi, a_imag_lo),
        "ddb": (b_real_hi, b_real_lo, b_imag_hi, b_imag_lo),
        "expected": (result_real_hi, result_real_lo, result_imag_hi, result_imag_lo)
    }

def generate_test_cases(num_cases):
    return [generate_ddcadd_test_case() for _ in range(num_cases)]

def write_test_cases_to_text_file(filename, test_cases):
    """
    Generate test cases for ddcadd and write them to a file.
    Each line contains:
    hi_a_real lo_a_real hi_a_imag lo_a_imag hi_b_real lo_b_real hi_b_imag lo_b_imag expected_hi_real expected_lo_real expected_hi_imag expected_lo_imag
    """
    with open(filename, "w") as f:
        for case in test_cases:
            f.write(
                f"{case['dda'][0]:.16e} {case['dda'][1]:.16e} "
                f"{case['dda'][2]:.16e} {case['dda'][3]:.16e} "
                f"{case['ddb'][0]:.16e} {case['ddb'][1]:.16e} "
                f"{case['ddb'][2]:.16e} {case['ddb'][3]:.16e} "
                f"{case['expected'][0]:.16e} {case['expected'][1]:.16e} "
                f"{case['expected'][2]:.16e} {case['expected'][3]:.16e}\n"
            )
    print(f"Test cases successfully written to {filename}")

def write_test_cases_to_binary(filename, test_cases):
    """
    Save test cases to a binary file.
    Each case: [hi_a_real, lo_a_real, hi_a_imag, lo_a_imag, hi_b_real, lo_b_real, hi_b_imag, lo_b_imag, expected_hi_real, expected_lo_real, expected_hi_imag, expected_lo_imag].
    """
    with open(filename, "wb") as f:
        for case in test_cases:
            f.write(struct.pack(
                "12d",  # 12 double-precision floats
                case["dda"][0], case["dda"][1], case["dda"][2], case["dda"][3],
                case["ddb"][0], case["ddb"][1], case["ddb"][2], case["ddb"][3],
                case["expected"][0], case["expected"][1], case["expected"][2], case["expected"][3]
            ))
    print(f"Test cases successfully written to {filename}")

if __name__ == "__main__":
    print("Generating test cases for ddcadd...")
    ntests = 10
    if len(sys.argv) > 1:
        ntests = int(sys.argv[1])
    test_cases = generate_test_cases(ntests)
    write_test_cases_to_text_file("data/ddcadd_test_cases.txt", test_cases)
    write_test_cases_to_binary("data/ddcadd_test_cases.bin", test_cases)