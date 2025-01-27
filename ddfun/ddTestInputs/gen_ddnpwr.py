# gen_ddnpwr.py
import random
import struct
import math
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

def generate_ddnpwr_test_case():
   """
   Generate a test case for ddnpwr.
   """
   # Generate random base in the range and an integer exponent
   base = mp.rand() * 100 - 50  # range [-50, 50)
   exponent = random.randint(-10, 10)

   # Compute expected result
   if base == 0 and exponent < 0:
      return None  # Skip invalid test case
   if exponent == 0:
      result = mp.mpf(1)
   elif exponent > 0:
      result = mp.power(base, exponent)
   else:
      result = mp.fdiv(mp.mpf(1), mp.power(base, abs(exponent)))

   # Convert to inputs and outputs to double-double format
   base_hi, base_lo = generate_double_double(base)
   result_hi, result_lo = generate_double_double(result)

   # Verify correctness of double-double representation
   tolerance = 30
   valid_base, error_base = verify_double_double(base, base_hi, base_lo, tolerance)
   valid_result, error_result = verify_double_double(result, result_hi, result_lo, tolerance)

   if not (valid_base and valid_result):
      print(
         f"Verification of ddnpwr failed for tolerance {tolerance}:\n"
         f"  base: {base} (hi={base_hi}, lo={base_lo}, error={error_base})\n"
         f"  result: {result} (hi={result_hi}, lo={result_lo}, error={error_result})"
      )
      return None

   return {
      "base": (base_hi, base_lo),
      "exponent": exponent,
      "expected": (result_hi, result_lo),
   }

def generate_test_cases(num_cases):
   test_cases = []
   while len(test_cases) < num_cases:
      case = generate_ddnpwr_test_case()
      if case is not None:
         test_cases.append(case)
   return test_cases

def write_test_cases_to_text_file(filename, test_cases):
   """
   Generate test cases for ddnpwr and write them to a file.
   Each line contains:
   hi_base lo_base exponent expected_hi expected_lo
   """
   with open(filename, "w") as f:
      for case in test_cases:
         f.write(
               f"{case['base'][0]:.16e} {case['base'][1]:.16e} "
               f"{case['exponent']} "
               f"{case['expected'][0]:.16e} {case['expected'][1]:.16e}\n"
         )
   print(f"Test cases successfully written to {filename}")

def write_test_cases_to_binary(filename, test_cases):
   """
   Save test cases to a binary file.
   Each case: [hi_base, lo_base, exponent, expected_hi, expected_lo].
   """
   with open(filename, "wb") as f:
      for case in test_cases:
         f.write(struct.pack(
               # write 2 double-precision floats and an integer and 2 double-precision floats
               "2di2d",
               case["base"][0], case["base"][1],
               case["exponent"],
               case["expected"][0], case["expected"][1]
         ))
   print(f"Test cases successfully written to {filename}")

if __name__ == "__main__":
   print("Generating test cases for ddnpwr...")
   test_cases = generate_test_cases(10)
   write_test_cases_to_text_file("data/ddnpwr_test_cases.txt", test_cases)
   write_test_cases_to_binary("data/ddnpwr_test_cases.bin", test_cases)
