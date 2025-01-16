import random
import struct
from decimal import Decimal, getcontext

# Set high precision for calculations
getcontext().prec = 50
PI = Decimal("3.14159265358979323846264338327950288419716939937510")  # High-precision PI


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


def generate_ddmuldd_test_case():
    """
    Generate a test case for ddmuldd.
    """
    # Generate two random double-precision values
    a = Decimal(random.uniform(0.1, 10))
    b = Decimal(random.uniform(0.1, 10))

    # Compute expected result
    result = a * b

    # Convert result to double-double format
    result_hi, result_lo = generate_double_double(result)

    # Verify correctness of double-double representation
    valid_result, error_result = verify_double_double(result, result_hi, result_lo)

    if not valid_result:
        raise ValueError(
            f"Verification failed:\n"
            f"  a: {a}\n"
            f"  b: {b}\n"
            f"  result: {result} (hi={result_hi}, lo={result_lo}, error={error_result})"
        )

    return {
        "a": float(a),
        "b": float(b),
        "expected": (result_hi, result_lo),
    }

def generate_test_cases(num_cases):
    return [generate_ddmuldd_test_case() for _ in range(num_cases)]

def write_test_cases_to_text_file(filename, test_cases):
    """
    Generate test cases for ddmuldd and write them to a file.
    Each line contains:
      a b expected_hi expected_lo
    """
    with open(filename, "w") as f:
        for case in test_cases:
            f.write(
                f"{case['a']} {case['b']} "
                f"{case['expected'][0]} {case['expected'][1]}\\n"
            )
    print(f"Test cases successfully written to {filename}")


def write_test_cases_to_binary(filename, test_cases):
    """
    Save test cases for ddmuldd to a binary file.
    Each case: [a, b, expected_hi, expected_lo].
    """
    with open(filename, "wb") as f:
        for case in test_cases:
            f.write(struct.pack(
                "4d",  # 4 double-precision floats
                case["a"], case["b"],
                case["expected"][0], case["expected"][1]
            ))

    print(f"Test cases successfully written to {filename}")


if __name__ == "__main__":
    test_cases = generate_test_cases(10)
    write_test_cases_to_text_file("data/ddmuldd_test_cases.txt", test_cases)
    write_test_cases_to_binary("data/ddmuldd_test_cases.bin", test_cases)
