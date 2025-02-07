# gen_ddcpwr.py
import random
import struct
import math
import sys
from mpmath import mp, mpc

# Set high precision for calculations
mp.dps = 50

def generate_double_double_complex(base_value):
   """
   Generate a double-double representation for a complex number (hi, lo) of a high-precision complex value.
   """
   real_hi = float(base_value.real)
   real_lo = float(base_value.real - mp.mpf(real_hi))
   imag_hi = float(base_value.imag)
   imag_lo = float(base_value.imag - mp.mpf(imag_hi))
   return (real_hi, real_lo, imag_hi, imag_lo)

def calculate_exponent_difference(a, b):
   """
   Calculate the difference of scale between two double values or complex numbers.
   """
   if isinstance(a, complex) or isinstance(b, complex):
       a_abs = math.fabs(a.real) + math.fabs(a.imag)
       b_abs = math.fabs(b.real) + math.fabs(b.imag)
   else:
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

def generate_ddcpwr_test_case():
   """
   Generate a test case for ddcpwr using a random complex number as the base value and a random integer exponent.
   """
   # Generate a random base complex number and random integer exponent
   real_part = mp.rand() - 0.5
   imag_part = mp.rand() - 0.5
   a = mpc(real_part, imag_part)
   n = random.randint(-5, 5)

   # Compute expected result
   result = a**n

   # Convert to double-double complex format
   a_dd = generate_double_double_complex(a)
   result_dd = generate_double_double_complex(result)
   
   # Verify correctness of double-double complex representation
   tolerance = 30
   valid_a_real, error_a_real = verify_double_double(a.real, a_dd[0], a_dd[1], tolerance)
   valid_a_imag, error_a_imag = verify_double_double(a.imag, a_dd[2], a_dd[3], tolerance)
   valid_result_real, error_result_real = verify_double_double(result.real, result_dd[0], result_dd[1], tolerance)
   valid_result_imag, error_result_imag = verify_double_double(result.imag, result_dd[2], result_dd[3], tolerance)
   valid_a = valid_a_real and valid_a_imag
   valid_result = valid_result_real and valid_result_imag

   if not (valid_a and valid_result):
      print(
         f"Verification of ddcpwr failed for tolerance {tolerance}:\n"
         f"  a: {a} (real_hi={a_dd[0]}, real_lo={a_dd[1]}, imag_hi={a_dd[2]}, imag_lo={a_dd[3]}, error={error_a})\n"
         f"  result: {result} (real_hi={result_dd[0]}, real_lo={result_dd[1]}, imag_hi={result_dd[2]}, imag_lo={result_dd[3]}, error={error_result})"
      )

   return {
      "a": a_dd,
      "n": n,
      "expected": result_dd,
   }

def generate_test_cases(num_cases):
   return [generate_ddcpwr_test_case() for _ in range(num_cases)]

def write_test_cases_to_text_file(filename, test_cases):
   """
   Generate test cases for ddcpwr and write them to a file.
   Each line contains:
   real_hi_a real_lo_a imag_hi_a imag_lo_a n expected_real_hi expected_real_lo expected_imag_hi expected_imag_lo
   """
   with open(filename, "w") as f:
      for case in test_cases:
         f.write(
               f"{case['a'][0]:.16e} {case['a'][1]:.16e} {case['a'][2]:.16e} {case['a'][3]:.16e} {case['n']} "
               f"{case['expected'][0]:.16e} {case['expected'][1]:.16e} {case['expected'][2]:.16e} {case['expected'][3]:.16e}\n"
         )
   print(f"Test cases successfully written to {filename}")

def write_test_cases_to_binary(filename, test_cases):
   """
   Save test cases to a binary file.
   Each case: [real_hi_a, real_lo_a, imag_hi_a, imag_lo_a, n, expected_real_hi, expected_real_lo, expected_imag_hi, expected_imag_lo].
   """
   with open(filename, "wb") as f:
      for case in test_cases:
         f.write(struct.pack(
               "4d i 4d",  # 4 double-precision floats, 1 integer, 4 double-precision floats
               case["a"][0], case["a"][1], case["a"][2], case["a"][3], case["n"],
               case["expected"][0], case["expected"][1], case["expected"][2], case["expected"][3]
         ))
   print(f"Test cases successfully written to {filename}")

if __name__ == "__main__":
   print("Generating test cases for ddcpwr...")
   ntests = 10
   if len(sys.argv) > 1:
      ntests = int(sys.argv[1])
   test_cases = generate_test_cases(ntests)
   write_test_cases_to_text_file("data/ddcpwr_test_cases.txt", test_cases)
   write_test_cases_to_binary("data/ddcpwr_test_cases.bin", test_cases)