import struct
import math
import random
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

def generate_ddmuldd_test_case():
    """
    Generate a test case for ddmuldd using random values as base values.
    """
    # Generate two random scaling factors [0,1) in high precision
    a = random.random()  # Adjust range as needed
    b = random.random()  # Adjust range as needed

    # Compute expected result
    result = mp.fmul(mp.mpf(a), mp.mpf(b))

    # Convert to inputs and outputs to double-double format
    result_hi, result_lo = generate_double_double(result)

    # Verify correctness of double-double representation
    tolerance = 30
    valid_result, error_result = verify_double_double(result, result_hi, result_lo, tolerance)

    if not valid_result:
        print(
            f"Verification of ddmuldd failed for tolerance {tolerance}:\n"
            f"  a: {a} b: {b}\n"
            f"  result: {result} (hi={result_hi}, lo={result_lo}, error={error_result})"
        )

    return {
        "da": float(a),
        "db": float(b),
        "expected": (result_hi, result_lo),
    }

def generate_test_cases(num_cases):
    return [generate_ddmuldd_test_case() for _ in range(num_cases)]

def write_test_cases_to_text_file(filename, test_cases):
    """
    Generate test cases for ddmuldd and write them to a file.
    Each line contains:
    da db expected_hi expected_lo
    """
    with open(filename, "w") as f:
        for case in test_cases:
            f.write(
                f"{case['da']:.16e} {case['db']:.16e} "
                f"{case['expected'][0]:.16e} {case['expected'][1]:.16e}\n"
            )
    print(f"Test cases successfully written to {filename}")

def write_test_cases_to_binary(filename, test_cases):
    """
    Save test cases to a binary file.
    Each case: [da, db, expected_hi, expected_lo].
    """
    with open(filename, "wb") as f:
        for case in test_cases:
            f.write(struct.pack(
                "4d",  # 4 double-precision floats
                case["da"], case["db"],
                case["expected"][0], case["expected"][1]
            ))
    print(f"Test cases successfully written to {filename}")

if __name__ == "__main__":
   print("Generating test cases for ddmuldd...")
   ntests = 10
   if len(sys.argv) > 1:
      ntests = int(sys.argv[1])
   test_cases = generate_test_cases(ntests)
   write_test_cases_to_text_file("data/ddmuldd_test_cases.txt", test_cases)
   write_test_cases_to_binary("data/ddmuldd_test_cases.bin", test_cases)