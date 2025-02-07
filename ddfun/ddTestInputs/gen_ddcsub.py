# gen_ddcsub.py
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

def calculate_exponent_difference(a,b):
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

def generate_ddcsub_test_case():
   """
   Generate a test case for ddcsub using random complex numbers.
   """
   # Generate two random complex numbers
   a_real = mp.rand() * 2 - 1
   a_imag = mp.rand() * 2 - 1
   b_real = mp.rand() * 2 - 1
   b_imag = mp.rand() * 2 - 1

   a = a_real + a_imag*1j
   b = b_real + b_imag*1j

   # Compute expected result
   result = a - b

   # Convert to inputs and outputs to double-double complex format
   a_real_hi, a_real_lo = generate_double_double(a_real)
   a_imag_hi, a_imag_lo = generate_double_double(a_imag)
   b_real_hi, b_real_lo = generate_double_double(b_real)
   b_imag_hi, b_imag_lo = generate_double_double(b_imag)
   result_real_hi, result_real_lo = generate_double_double(result.real)
   result_imag_hi, result_imag_lo = generate_double_double(result.imag)

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
         f"Verification of ddcsub failed for tolerance {tolerance}:\n"
         f"  a: ({a_real}, {a_imag}) (real_hi={a_real_hi}, real_lo={a_real_lo}, imag_hi={a_imag_hi}, imag_lo={a_imag_lo}, error_real={error_a_real}, error_imag={error_a_imag})\n"
         f"  b: ({b_real}, {b_imag}) (real_hi={b_real_hi}, real_lo={b_real_lo}, imag_hi={b_imag_hi}, imag_lo={b_imag_lo}, error_real={error_b_real}, error_imag={error_b_imag})\n"
         f"  result: ({result.real}, {result.imag}) (real_hi={result_real_hi}, real_lo={result_real_lo}, imag_hi={result_imag_hi}, imag_lo={result_imag_lo}, error_real={error_result_real}, error_imag={error_result_imag})"
      )

   return {
      "a": ((a_real_hi, a_real_lo), (a_imag_hi, a_imag_lo)),
      "b": ((b_real_hi, b_real_lo), (b_imag_hi, b_imag_lo)),
      "expected": ((result_real_hi, result_real_lo), (result_imag_hi, result_imag_lo)),
   }


def generate_test_cases(num_cases):
   return [generate_ddcsub_test_case() for _ in range(num_cases)]

def write_test_cases_to_text_file(filename, test_cases):
   """
   Generate test cases for ddcsub and write them to a file.
   Each line contains:
   real_hi_a real_lo_a imag_hi_a imag_lo_a real_hi_b real_lo_b imag_hi_b imag_lo_b expected_real_hi expected_real_lo expected_imag_hi expected_imag_lo
   """
   with open(filename, "w") as f:
      for case in test_cases:
         f.write(
               f"{case['a'][0][0]:.16e} {case['a'][0][1]:.16e} "
               f"{case['a'][1][0]:.16e} {case['a'][1][1]:.16e} "
               f"{case['b'][0][0]:.16e} {case['b'][0][1]:.16e} "
               f"{case['b'][1][0]:.16e} {case['b'][1][1]:.16e} "
               f"{case['expected'][0][0]:.16e} {case['expected'][0][1]:.16e} "
               f"{case['expected'][1][0]:.16e} {case['expected'][1][1]:.16e}\n"
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
               case["a"][0][0], case["a"][0][1],
               case["a"][1][0], case["a"][1][1],
               case["b"][0][0], case["b"][0][1],
               case["b"][1][0], case["b"][1][1],
               case["expected"][0][0], case["expected"][0][1],
               case["expected"][1][0], case["expected"][1][1]
         ))
   print(f"Test cases successfully written to {filename}")

if __name__ == "__main__":
   print("Generating test cases for ddcsub...")
   ntests = 10
   if len(sys.argv) > 1:
      ntests = int(sys.argv[1])
   test_cases = generate_test_cases(ntests)
   write_test_cases_to_text_file("data/ddcsub_test_cases.txt", test_cases)
   write_test_cases_to_binary("data/ddcsub_test_cases.bin", test_cases)