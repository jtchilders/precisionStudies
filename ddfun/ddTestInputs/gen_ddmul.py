# gen_ddmul.py
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

def generate_ddmul_test_case():
   """
   Generate a test case for ddmul using PI and a random value as the base values.
   """
   # Generate two random scaling factors [-1,1) in high precision
   a = mp.rand() * 2 - 1  # Adjusted for a wider range
   b = mp.rand() * 2 - 1  # Adjusted for a wider range

   # Compute expected result
   result = mp.fmul(a, b)

   # Convert to inputs and outputs to double-double format
   a_hi, a_lo = generate_double_double(a)
   b_hi, b_lo = generate_double_double(b)
   result_hi, result_lo = generate_double_double(result)

   # Verify correctness of double-double representation
   tolerance = 30
   valid_a, error_a = verify_double_double(a, a_hi, a_lo, tolerance)
   valid_b, error_b = verify_double_double(b, b_hi, b_lo, tolerance)
   valid_result, error_result = verify_double_double(result, result_hi, result_lo, tolerance)

   if not (valid_a and valid_b and valid_result):
      print(
         f"Verification of ddmul failed for tolerance {tolerance}:\n"
         f"  a: {a} (hi={a_hi}, lo={a_lo}, error={error_a})\n"
         f"  b: {b} (hi={b_hi}, lo={b_lo}, error={error_b})\n"
         f"  result: {result} (hi={result_hi}, lo={result_lo}, error={error_result})"
      )

   return {
      "dda": (a_hi, a_lo),
      "ddb": (b_hi, b_lo),
      "expected": (result_hi, result_lo),
   }

def generate_test_cases(num_cases):
   return [generate_ddmul_test_case() for _ in range(num_cases)]

def write_test_cases_to_text_file(filename, test_cases):
   """
   Generate test cases for ddmul and write them to a file.
   Each line contains:
   hi_a lo_a hi_b lo_b expected_hi expected_lo
   """
   with open(filename, "w") as f:
      for case in test_cases:
         f.write(
               f"{case['dda'][0]:.16e} {case['dda'][1]:.16e} "
               f"{case['ddb'][0]:.16e} {case['ddb'][1]:.16e} "
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
               case["dda"][0], case["dda"][1],
               case["ddb"][0], case["ddb"][1],
               case["expected"][0], case["expected"][1]
         ))
   print(f"Test cases successfully written to {filename}")

if __name__ == "__main__":
   print("Generating test cases for ddmul...")
   ntests = 10
   if len(sys.argv) > 1:
      ntests = int(sys.argv[1])
   test_cases = generate_test_cases(ntests)
   write_test_cases_to_text_file("data/ddmul_test_cases.txt", test_cases)
   write_test_cases_to_binary("data/ddmul_test_cases.bin", test_cases)