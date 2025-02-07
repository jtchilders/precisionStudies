# gen_ddcdiv.py
import random
import struct
import math
import sys
from mpmath import mp

# Set high precision for calculations
mp.dps = 50

def generate_double_double_complex(real_part, imag_part):
    """
    Generate a double-double complex representation from high-precision real and imaginary parts.
    """
    real_hi = float(real_part)
    real_lo = float(real_part - mp.mpf(real_hi))
    imag_hi = float(imag_part)
    imag_lo = float(imag_part - mp.mpf(imag_hi))
    return (real_hi, real_lo, imag_hi, imag_lo)

def verify_double_double_complex(original_real, original_imag, hi_real, lo_real, hi_imag, lo_imag, tolerance):
    """
    Verify that hi and lo accurately represent the original value in double-double complex form.
    """
    reconstructed_real = mp.mpf(hi_real) + mp.mpf(lo_real)
    reconstructed_imag = mp.mpf(hi_imag) + mp.mpf(lo_imag)
    error_real = abs(original_real - reconstructed_real)
    error_imag = abs(original_imag - reconstructed_imag)

    scale_diff_real = calculate_exponent_difference(original_real, error_real)
    scale_diff_imag = calculate_exponent_difference(original_imag, error_imag)

    return scale_diff_real > tolerance and scale_diff_imag > tolerance, (error_real, error_imag)

def calculate_exponent_difference(a, b):
    """
    Calculate the difference of scale between two values
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

def generate_ddcdiv_test_case():
    """
    Generate a test case for ddcdiv using random complex numbers.
    """
    # Generate random complex numbers for input
    a_real = mp.rand()
    a_imag = mp.rand()
    b_real = mp.rand()
    b_imag = mp.rand()

    while b_real == 0.0 and b_imag == 0.0:
        b_real = mp.rand()
        b_imag = mp.rand()

    a = mp.mpc(a_real, a_imag)
    b = mp.mpc(b_real, b_imag)

    # Compute expected result using high precision division
    result = a / b

    # Convert inputs and outputs to double-double complex format
    a_components = generate_double_double_complex(a.real, a.imag)
    b_components = generate_double_double_complex(b.real, b.imag)
    result_components = generate_double_double_complex(result.real, result.imag)

    # Verify correctness of double-double complex representation
    tolerance = 30
    valid_a, error_a = verify_double_double_complex(a.real, a.imag, *a_components, tolerance)
    valid_b, error_b = verify_double_double_complex(b.real, b.imag, *b_components, tolerance)
    valid_result, error_result = verify_double_double_complex(result.real, result.imag, *result_components, tolerance)

    if not (valid_a and valid_b and valid_result):
        print(
            f"Verification of ddcdiv failed for tolerance {tolerance}:\n"
            f"  a: {a} (components={a_components}, error={error_a})\n"
            f"  b: {b} (components={b_components}, error={error_b})\n"
            f"  result: {result} (components={result_components}, error={error_result})"
        )

    return {
        "a": a_components,
        "b": b_components,
        "expected": result_components,
    }

def generate_test_cases(num_cases):
    return [generate_ddcdiv_test_case() for _ in range(num_cases)]

def write_test_cases_to_text_file(filename, test_cases):
    """
    Generate test cases for ddcdiv and write them to a file.
    Each line contains:
    [hi_real_a lo_real_a hi_imag_a lo_imag_a hi_real_b lo_real_b hi_imag_b lo_imag_b expected_hi_real expected_lo_real expected_hi_imag expected_lo_imag]
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
    Each case: [hi_real_a, lo_real_a, hi_imag_a, lo_imag_a, hi_real_b, lo_real_b, hi_imag_b, lo_imag_b, expected_hi_real, expected_lo_real, expected_hi_imag, expected_lo_imag].
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
    print("Generating test cases for ddcdiv...")
    ntests = 10
    if len(sys.argv) > 1:
        ntests = int(sys.argv[1])
    test_cases = generate_test_cases(ntests)
    write_test_cases_to_text_file("data/ddcdiv_test_cases.txt", test_cases)
    write_test_cases_to_binary("data/ddcdiv_test_cases.bin", test_cases)