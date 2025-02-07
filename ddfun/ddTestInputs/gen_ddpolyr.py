# gen_ddpolyr.py
import random
import struct
import math
import sys
from mpmath import mp


# Set high precision for calculations
mp.dps = 50

def ddpolyr(n, a, x0, max_iterations=20, eps=mp.mpf('1e-29')):
   """
   Finds the root of a polynomial near x0 using Newton's method with high precision.

   Parameters:
      n (int): Degree of the polynomial.
      a (list of mpf): Coefficients of the polynomial in ascending order (a[0] + a[1]*x + ...).
      x0 (mpf): Initial guess for the root.
      max_iterations (int): Maximum number of iterations for convergence.
      eps (mpf): Convergence threshold.

   Returns:
      mpf: The computed root of the polynomial.
   """

   # Compute derivative coefficients
   ad = [a[i] * i for i in range(1, n + 1)]

   # Initialize root estimate
   x = mp.mpf(x0)

   # Newton's method
   for _ in range(max_iterations):
      # Compute f(x) and f'(x)
      fx = sum(a[i] * (x ** i) for i in range(n + 1))
      f_prime_x = sum(ad[i - 1] * (x ** (i - 1)) for i in range(1, n + 1))

      # Avoid division by zero
      if f_prime_x == 0:
         print("Warning: Derivative is zero. Newton's method fails.")
         return 0

      # Update x using Newton's formula
      delta_x = fx / f_prime_x
      x -= delta_x

      # Check for convergence
      if abs(delta_x) <= eps:
         return x

   # If no convergence, raise an error
   print("Newton's method did not converge within the maximum number of iterations.")
   return 0


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

def generate_ddpolyr_test_case(degree):
   """
   Generate a test case for ddpolyr function.
   """
   # Generate random coefficients for the polynomial
   coefficients = [mp.rand() for _ in range(degree + 1)]
   # Generate a random initial guess for the root
   x0 = mp.rand()

   
   # Compute expected result (root) using a known root or solver
   # Note: Finding exact roots might require solving the polynomial using numerical methods
   # For now, assume the root is near the x0 (initial guess)
   expected_root = ddpolyr(degree, coefficients, x0)

   # Convert inputs and expected output to double-double format
   coefficients_double_double = [generate_double_double(coef) for coef in coefficients]
   x0_hi, x0_lo = generate_double_double(x0)
   expected_hi, expected_lo = generate_double_double(expected_root)

   # Verify correctness of double-double representation
   tolerance = 30
   for i, coef in enumerate(coefficients):
      valid, error = verify_double_double(coef, coefficients_double_double[i][0], coefficients_double_double[i][1], tolerance)
      if not valid:
         print(f"Verification of coefficient {i} failed for tolerance {tolerance}.")

   valid_x0, error_x0 = verify_double_double(x0, x0_hi, x0_lo, tolerance)
   valid_expected, error_expected = verify_double_double(expected_root, expected_hi, expected_lo, tolerance)

   if not (valid_x0 and valid_expected):
      print(
         f"Verification of ddpolyr failed for tolerance {tolerance}:\n"
         f"  x0: {x0} (hi={x0_hi}, lo={x0_lo}, error={error_x0})\n"
         f"  expected_root: {expected_root} (hi={expected_hi}, lo={expected_lo}, error={error_expected})"
      )

   return {
      "coefficients": coefficients_double_double,
      "x0": (x0_hi, x0_lo),
      "expected": (expected_hi, expected_lo),
   }

def generate_test_cases(num_cases, degree):
   return [generate_ddpolyr_test_case(degree) for _ in range(num_cases)]

def write_test_cases_to_text_file(filename, test_cases):
   """
   Generate test cases for ddpolyr and write them to a file.
   Each line contains:
   degree coeffs_hi coeffs_lo .... x0_hi x0_lo expected_hi expected_lo
   """
   with open(filename, "w") as f:
      for case in test_cases:
         coeffs_str = ' '.join(
            f"{coef[0]:.16e} {coef[1]:.16e}" for coef in case['coefficients']
         )
         f.write(
            f"{len(case['coefficients']) - 1} {coeffs_str} "
            f"{case['x0'][0]:.16e} {case['x0'][1]:.16e} "
            f"{case['expected'][0]:.16e} {case['expected'][1]:.16e}\n"
         )
   print(f"Test cases successfully written to {filename}")

def write_test_cases_to_binary(filename, test_cases):
   """
   Save test cases to a binary file.
   Each case: degree, coeffs_hi, coeffs_lo, x0_hi, x0_lo, expected_hi, expected_lo.
   """
   with open(filename, "wb") as f:
      for case in test_cases:
         f.write(struct.pack(
            "i" + "2d" * len(case['coefficients']) + "4d",  # degree + 2*len(coeffs) doubles + 4 doubles
            len(case['coefficients']) - 1,
            *(elem for coef in case['coefficients'] for elem in coef),
            case['x0'][0], case['x0'][1],
            case['expected'][0], case['expected'][1]
         ))
   print(f"Test cases successfully written to {filename}")

if __name__ == "__main__":
   print("Generating test cases for ddpolyr...")
   degree = 3  # Choose the polynomial degree
   ntests = 10
   if len(sys.argv) > 1:
      ntests = int(sys.argv[1])
   test_cases = generate_test_cases(ntests, degree)
   write_test_cases_to_text_file("data/ddpolyr_test_cases.txt", test_cases)
   write_test_cases_to_binary("data/ddpolyr_test_cases.bin", test_cases)