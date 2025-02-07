# gen_ddcmul.py
import random
import struct
import math
import sys
from mpmath import mp, mpc

# Set high precision for calculations
mp.dps = 100

def generate_double_double_complex(base_value):
    """
    Generate a double-double representation (hi, lo) for complex numbers.
    """
    real_hi = float(base_value.real)
    real_lo = float(base_value.real - mp.mpf(real_hi))
    imag_hi = float(base_value.imag)
    imag_lo = float(base_value.imag - mp.mpf(imag_hi))
    return (real_hi, real_lo, imag_hi, imag_lo)

def calculate_exponent_difference(a,b):
    """
    Calculate the difference of scale between two double values
    """
    a_abs = abs(a)
    b_abs = abs(b)
    try:
        a_exp = int(math.log10(float(a_abs)))
    except ValueError:
        a_exp = 0
    try:
        b_exp = int(math.log10(float(b_abs)))
    except ValueError:
        b_exp = 0
    out = abs(a_exp - b_exp)
    return out

def verify_double_double_complex(original, hi_real, lo_real, hi_imag, lo_imag, tolerance):
    """
    Verify that hi and lo accurately represent the original complex value.
    """
    reconstructed = mp.mpc(hi_real, hi_imag) + mp.mpc(lo_real, lo_imag)
    error = abs(original - reconstructed)
    if error == 0.0:
        return True, error
    scale_diff = calculate_exponent_difference(original, error)
    return scale_diff > tolerance, error

def generate_ddcmul_test_case():
    """
    Generate a test case for ddcmul using random complex numbers as base values.
    """
    # Generate two random complex numbers
    a = mpc(mp.rand(), mp.rand())
    b = mpc(mp.rand(), mp.rand())

    # Compute expected result
    result = a * b

    # Convert inputs and outputs to double-double complex format
    a_hi_real, a_lo_real, a_hi_imag, a_lo_imag = generate_double_double_complex(a)
    b_hi_real, b_lo_real, b_hi_imag, b_lo_imag = generate_double_double_complex(b)
    result_hi_real, result_lo_real, result_hi_imag, result_lo_imag = generate_double_double_complex(result)

    # Verify correctness of double-double representation
    tolerance = 30
    valid_a, error_a = verify_double_double_complex(a, a_hi_real, a_lo_real, a_hi_imag, a_lo_imag, tolerance)
    valid_b, error_b = verify_double_double_complex(b, b_hi_real, b_lo_real, b_hi_imag, b_lo_imag, tolerance)
    valid_result, error_result = verify_double_double_complex(result, result_hi_real, result_lo_real, result_hi_imag, result_lo_imag, tolerance)

    if not (valid_a and valid_b and valid_result):
        print(
            f"Verification of ddcmul failed for tolerance {tolerance}:\n"
            f"  a: {a} (hi=({a_hi_real}, {a_hi_imag}), lo=({a_lo_real}, {a_lo_imag}), error={error_a})\n"
            f"  b: {b} (hi=({b_hi_real}, {b_hi_imag}), lo=({b_lo_real}, {b_lo_imag}), error={error_b})\n"
            f"  result: {result} (hi=({result_hi_real}, {result_hi_imag}), lo=({result_lo_real}, {result_lo_imag}), error={error_result})"
        )

    return {
        "a": (a_hi_real, a_lo_real, a_hi_imag, a_lo_imag),
        "b": (b_hi_real, b_lo_real, b_hi_imag, b_lo_imag),
        "expected": (result_hi_real, result_lo_real, result_hi_imag, result_lo_imag),
    }

def generate_test_cases(num_cases):
    return [generate_ddcmul_test_case() for _ in range(num_cases)]

def write_test_cases_to_text_file(filename, test_cases):
    """
    Generate test cases for ddcmul and write them to a file.
    Each line contains:
    real_hi_a, real_lo_a, imag_hi_a, imag_lo_a, real_hi_b, real_lo_b, imag_hi_b, imag_lo_b, expected_real_hi, expected_real_lo, expected_imag_hi, expected_imag_lo
    """
    with open(filename, "w") as f:
        for case in test_cases:
            f.write(
                f"{case['a'][0]:.16e} {case['a'][1]:.16e} {case['a'][2]:.16e} {case['a'][3]:.16e} "
                f"{case['b'][0]:.16e} {case['b'][1]:.16e} {case['b'][2]:.16e} {case['b'][3]:.16e} "
                f"{case['expected'][0]:.16e} {case['expected'][1]:.16e} {case['expected'][2]:.16e} {case['expected'][3]:.16e}\n"
            )
    print(f"Test cases successfully written to {filename}")

def write_test_cases_to_binary(filename, test_cases):
    """
    Save test cases to a binary file.
    Each case: [real_hi_a, real_lo_a, imag_hi_a, imag_lo_a, real_hi_b, real_lo_b, imag_hi_b, imag_lo_b, expected_real_hi, expected_real_lo, expected_imag_hi, expected_imag_lo].
    """
    with open(filename, "wb") as f:
        for case in test_cases:
            f.write(struct.pack(
                "12d",  # 12 double-precision floats
                case["a"][0], case["a"][1], case["a"][2], case["a"][3],
                case["b"][0], case["b"][1], case["b"][2], case["b"][3],
                case["expected"][0], case["expected"][1], case["expected"][2], case["expected"][3]
            ))
    print(f"Test cases successfully written to {filename}")

if __name__ == "__main__":
    print("Generating test cases for ddcmul...")
    ntests = 10
    if len(sys.argv) > 1:
        ntests = int(sys.argv[1])
    test_cases = generate_test_cases(ntests)
    write_test_cases_to_text_file("data/ddcmul_test_cases.txt", test_cases)
    write_test_cases_to_binary("data/ddcmul_test_cases.bin", test_cases)