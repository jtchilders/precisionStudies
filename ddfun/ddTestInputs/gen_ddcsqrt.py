# gen_ddcsqrt.py
import random
import struct
import math
import sys
from mpmath import mp, sqrt

# Set high precision for calculations
mp.dps = 50

def generate_double_double_complex(real_part, imag_part):
    """
    Generate a double-double representation of a complex number with high-precision values for real and imaginary parts.
    """
    real_hi = float(real_part)
    real_lo = float(real_part - mp.mpf(real_hi))
    imag_hi = float(imag_part)
    imag_lo = float(imag_part - mp.mpf(imag_hi))
    return (real_hi, real_lo), (imag_hi, imag_lo)

def verify_double_double_complex(original, hi_lo_real, hi_lo_imag, tolerance):
    """
    Verify that hi_lo_real and hi_lo_imag accurately represent the original complex value.
    """
    reconstructed_real = mp.mpf(hi_lo_real[0]) + mp.mpf(hi_lo_real[1])
    reconstructed_imag = mp.mpf(hi_lo_imag[0]) + mp.mpf(hi_lo_imag[1])
    error_real = abs(original.real - reconstructed_real)
    error_imag = abs(original.imag - reconstructed_imag)
    scale_diff_real = calculate_exponent_difference(original.real, error_real)
    scale_diff_imag = calculate_exponent_difference(original.imag, error_imag)
    return scale_diff_real > tolerance and scale_diff_imag > tolerance, (error_real, error_imag)

def calculate_exponent_difference(a, b):
    """
    Calculate the difference of scale between two double values.
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

def generate_ddcsqrt_test_case():
    """
    Generate a test case for ddcsqrt using random complex numbers.
    """
    # Generate two random components for the complex number A
    real_a = mp.rand() - 0.5
    imag_a = mp.rand() - 0.5
    a = mp.mpc(real=real_a, imag=imag_a)
    
    # Compute expected result
    magnitude = sqrt(real_a**2 + imag_a**2)
    sqrt_real = sqrt((magnitude + real_a) / 2)
    sqrt_imag = sqrt((magnitude - real_a) / 2) * (1 if imag_a >= 0 else -1)
    result = mp.mpc(real=sqrt_real, imag=sqrt_imag)

    # Convert to inputs and outputs to double-double format
    a_hi_lo_real, a_hi_lo_imag = generate_double_double_complex(real_a, imag_a)
    result_hi_lo_real, result_hi_lo_imag = generate_double_double_complex(sqrt_real, sqrt_imag)

    # Verify correctness of double-double representation
    tolerance = 30
    valid_a, error_a = verify_double_double_complex(a, a_hi_lo_real, a_hi_lo_imag, tolerance)
    valid_result, error_result = verify_double_double_complex(result, result_hi_lo_real, result_hi_lo_imag, tolerance)

    if not (valid_a and valid_result):
        print(
            f"Verification of ddcsqrt failed for tolerance {tolerance}:\n"
            f"  a: {a} (real_hi={a_hi_lo_real[0]}, real_lo={a_hi_lo_real[1]}, imag_hi={a_hi_lo_imag[0]}, imag_lo={a_hi_lo_imag[1]}, error_real={error_a[0]}, error_imag={error_a[1]})\n"
            f"  result: {result} (real_hi={result_hi_lo_real[0]}, real_lo={result_hi_lo_real[1]}, imag_hi={result_hi_lo_imag[0]}, imag_lo={result_hi_lo_imag[1]}, error_real={error_result[0]}, error_imag={error_result[1]})"
        )

    return {
        "a_real": (a_hi_lo_real[0], a_hi_lo_real[1]),
        "a_imag": (a_hi_lo_imag[0], a_hi_lo_imag[1]),
        "expected_real": (result_hi_lo_real[0], result_hi_lo_real[1]),
        "expected_imag": (result_hi_lo_imag[0], result_hi_lo_imag[1]),
    }

def generate_test_cases(num_cases):
    return [generate_ddcsqrt_test_case() for _ in range(num_cases)]

def write_test_cases_to_text_file(filename, test_cases):
    """
    Generate test cases for ddcsqrt and write them to a file.
    Each line contains:
    hi_a_real lo_a_real hi_a_imag lo_a_imag expected_hi_real expected_lo_real expected_hi_imag expected_lo_imag
    """
    with open(filename, "w") as f:
        for case in test_cases:
            f.write(
                f"{case['a_real'][0]:.16e} {case['a_real'][1]:.16e} "
                f"{case['a_imag'][0]:.16e} {case['a_imag'][1]:.16e} "
                f"{case['expected_real'][0]:.16e} {case['expected_real'][1]:.16e} "
                f"{case['expected_imag'][0]:.16e} {case['expected_imag'][1]:.16e}\n"
            )
    print(f"Test cases successfully written to {filename}")

def write_test_cases_to_binary(filename, test_cases):
    """
    Save test cases to a binary file.
    Each case: [hi_a_real, lo_a_real, hi_a_imag, lo_a_imag, expected_hi_real, expected_lo_real, expected_hi_imag, expected_lo_imag].
    """
    with open(filename, "wb") as f:
        for case in test_cases:
            f.write(struct.pack(
                "8d",  # 8 double-precision floats
                case["a_real"][0], case["a_real"][1],
                case["a_imag"][0], case["a_imag"][1],
                case["expected_real"][0], case["expected_real"][1],
                case["expected_imag"][0], case["expected_imag"][1]
            ))
    print(f"Test cases successfully written to {filename}")

if __name__ == "__main__":
    print("Generating test cases for ddcsqrt...")
    ntests = 10
    if len(sys.argv) > 1:
        ntests = int(sys.argv[1])
    test_cases = generate_test_cases(ntests)
    write_test_cases_to_text_file("data/ddcsqrt_test_cases.txt", test_cases)
    write_test_cases_to_binary("data/ddcsqrt_test_cases.bin", test_cases)