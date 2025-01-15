# python file that uses GPT to generate C++ code
# the generated functions will be math functions that
# operate on a ddouble struct containing hi and lo

import re
import os
from openai import OpenAI


# Initialize the OpenAI client
client = OpenAI()

# OpenAI API key setup via environment variable OPENAI_API_KEY

# list of functions to extract from the fortran file
functions_to_extract = [
#    "ddadd",
#    "ddsub",
#    "ddmul",
#    "dddiv",
#    "ddmuld",
#    "ddmuldd",
#    "ddsqrt",
   # "ddexp",
   # "ddlog",
   # "ddpower",
   "ddnpwr",
#    "ddlog10",
#    "ddlog1p",
#    "ddlog2",
#    "ddlogb",
#    "ddexpm1",
#    "ddlog1p",
]

gen_code = """# gen_ddadd.py
import random
import struct
import numpy as np
from decimal import Decimal, getcontext

# Set high precision for calculations
getcontext().prec = 50
PI = Decimal("3.14159265358979323846264338327950288419716939937510")  # High-precision PI


def generate_double_double(base_value):
    \"""
    Generate a double-double representation (hi, lo) of a high-precision value.
    \"""
    hi = float(base_value)
    lo = float(base_value - Decimal(hi))  # Compute the remainder
    return hi, lo


def verify_double_double(original, hi, lo, tolerance=1e-30):
    \"""
    Verify that hi and lo accurately represent the original value.
    \"""
    reconstructed = Decimal(hi) + Decimal(lo)
    error = abs(original - reconstructed)
    return error < Decimal(tolerance), error


def generate_ddadd_test_case():
    \"""
    Generate a test case for ddadd using PI as the base value.
    \"""
    # Generate two random scaling factors
    scale_a = Decimal(random.uniform(0.1, 10))  # Scale factor for first number
    scale_b = Decimal(random.uniform(0.1, 10))  # Scale factor for second number

    # Derive the high-precision values from PI
    a = PI * scale_a
    b = PI / scale_b  # Ensure variety in operations

    # Compute expected result
    result = a + b

    # Convert to double-double format
    a_hi, a_lo = generate_double_double(a)
    b_hi, b_lo = generate_double_double(b)
    result_hi, result_lo = generate_double_double(result)

    # Verify correctness of double-double representation
    valid_a, error_a = verify_double_double(a, a_hi, a_lo)
    valid_b, error_b = verify_double_double(b, b_hi, b_lo)
    valid_result, error_result = verify_double_double(result, result_hi, result_lo)

    if not (valid_a and valid_b and valid_result):
        print(
            f"Verification failed:\\n"
            f"  a: {{a}} (hi={{a_hi}}, lo={{a_lo}}, error={{error_a}})\\n"
            f"  b: {{b}} (hi={{b_hi}}, lo={{b_lo}}, error={{error_b}})\\n"
            f"  result: {{result}} (hi={{result_hi}}, lo={{result_lo}}, error={{error_result}})"
        )

    return {
        "dda": (a_hi, a_lo),
        "ddb": (b_hi, b_lo),
        "expected": (result_hi, result_lo),
    }


def generate_test_cases(num_cases):
    return [generate_ddadd_test_case() for _ in range(num_cases)]

def write_test_cases_to_text_file(filename, test_cases):
    \"""
    Generate test cases for ddadd and write them to a file.
    Each line contains:
      hi_a lo_a hi_b lo_b expected_hi expected_lo
    \"""
    with open(filename, "w") as f:
        for case in test_cases:
            f.write(
                f"{{case['dda'][0]:.16e}} {{case['dda'][1]:.16e}} "
                f"{{case['ddb'][0]:.16e}} {{case['ddb'][1]:.16e}} "
                f"{{case['expected'][0]:.16e}} {{case['expected'][1]:.16e}}\\n"
            )
    print(f"Test cases successfully written to {{filename}}")

def write_test_cases_to_binary(filename, test_cases):
    \"""
    Save test cases to a binary file.
    Each case: [hi_a, lo_a, hi_b, lo_b, expected_hi, expected_lo].
    \"""
    with open(filename, "wb") as f:
        for case in test_cases:
            f.write(struct.pack(
                "6d",  # 6 double-precision floats
                case["dda"][0], case["dda"][1],
                case["ddb"][0], case["ddb"][1],
                case["expected"][0], case["expected"][1]
            ))
    print(f"Test cases successfully written to {{filename}}")

if __name__ == "__main__":
    test_cases = generate_test_cases(10)
    write_test_cases_to_text_file("ddadd_test_cases.txt", test_cases)
    write_test_cases_to_binary("ddadd_test_cases.bin", test_cases)
"""


test_code = """! test_ddadd.f90
module test_ddadd
  use ddfuna  ! Import the module containing the ddadd subroutine
  implicit none

contains

  subroutine unittest_ddadd(filename)
    implicit none
    character(len=*), intent(in) :: filename
    real(8) :: hi_a, lo_a, hi_b, lo_b, expected_hi, expected_lo
    real(8) :: dda(2), ddb(2), ddc(2), expected_ddc(2)
    real(8) :: tolerance
    logical :: test_passed
    integer :: iunit, ios, total_tests, passed_tests

    tolerance = 1.0e-30  ! Double-double precision tolerance
    total_tests = 0
    passed_tests = 0

    ! Open the binary file
    open(unit=10, file=filename, form="unformatted", access="stream", status="old", action="read")
    write(*,*) "Running tests from:", filename

    do
      ! Read six double-precision values (one test case)
      read(10, iostat=ios) hi_a, lo_a, hi_b, lo_b, expected_hi, expected_lo
      if (ios /= 0) exit  ! Exit loop at end of file

      total_tests = total_tests + 1

      ! Pack inputs
      dda(1) = hi_a
      dda(2) = lo_a
      ddb(1) = hi_b
      ddb(2) = lo_b

      ! Expected result
      expected_ddc(1) = expected_hi
      expected_ddc(2) = expected_lo

      ! Call the ddadd subroutine
      call ddadd(dda, ddb, ddc)

      ! Compare results with expected values
      test_passed = abs(ddc(1) - expected_ddc(1)) < tolerance .and. &
                    abs(ddc(2) - expected_ddc(2)) < tolerance

      ! Print results
      if (test_passed) then
        passed_tests = passed_tests + 1
      else
        write(*,*) "Test Failed:", &
                   "inputs: [", dda(1), ", ", dda(2), "] + [", ddb(1), ", ", ddb(2), "]", &
                   "result: [", ddc(1), ", ", ddc(2), "]", &
                   "expected: [", expected_ddc(1), ", ", expected_ddc(2), "]", &
                   "error: [", abs(ddc(1) - expected_ddc(1)), ", ", abs(ddc(2) - expected_ddc(2)), "]"
      end if
    end do

    close(10)

    ! Summary
    write(*,*) "Tests completed:", passed_tests, "passed out of", total_tests

  end subroutine unittest_ddadd

end module test_ddadd
"""


ddadd_func = """subroutine ddadd (dda, ddb, ddc)

!   This subroutine computes ddc = dda + ddb, where dda, ddb and ddc are type DDR.

implicit none
real (ddknd), intent(in):: dda(2), ddb(2)
real (ddknd), intent(out):: ddc(2)
real (ddknd) e, t1, t2

!   Compute dda + ddb using Knuth's trick.

t1 = dda(1) + ddb(1)
e = t1 - dda(1)
t2 = ((ddb(1) - e) + (dda(1) - (t1 - e))) + dda(2) + ddb(2)

!   The result is t1 + t2, after normalization.

ddc(1) = t1 + t2
ddc(2) = t2 - (ddc(1) - t1)
return
end subroutine ddadd
"""

def extract_functions_with_bodies(filename):
   with open(filename, 'r') as file:
      code = file.read()

   # Compile the regex pattern
   function_pattern = re.compile(
      r'\b(subroutine|function)\s+(\w+)\s*\((.*?)\)\s*'  # Function signature
      r'(.*?)'                                          # Function body (lazy match)
      r'\b(end\s+\1)',                                  # Match the corresponding end
      re.DOTALL | re.IGNORECASE                        # Dot matches newlines; case-insensitive
   )

   # Use the compiled pattern to find matches
   matches = function_pattern.findall(code)

   functions = []
   for match in matches:
      func_type, name, args, body, _ = match
      functions.append({
         "type": func_type.strip(),
         "name": name.strip(),
         "args": [arg.strip() for arg in args.split(',') if arg.strip()],
         "body": f"{func_type} {name}({', '.join(args.split())})\n{body}\nend {func_type}"
      })

   return functions


def generate_gen_with_gpt(function):
   # Create a prompt for GPT
   prompt = f"""

I am porting the ddfun library from Fortran to C++.
The ddfun library contains subroutines for performing arithmetic operations on double-double numbers.

Given this example ddadd function:
---------------
{ddadd_func}
---------------
I have the following python file that generates test data for the ddadd function:
----------------
{gen_code}
----------------
Then I have the fortran unit test for the ddadd function:
----------------
{test_code}
----------------

I need a similar python file that generates test data for the {function} function, which is defined as follows:

{function['body']}

Can you please generate a generator following the same format as in the ddadd case? Please avoid any extra text or markdown formatting because your reply will go directly into a file to be compiled.
Do not even include the markdown formatting around the code.
"""   

   # Call the GPT API
   completion = client.chat.completions.create(
      model="gpt-4o",  # Use the desired model
      store=True,
      messages=[
         {"role": "user", "content": prompt}
      ]
   )

   # Extract and return the generated code
   return completion.choices[0].message.content


def generate_test_with_gpt(function):
   # Create a prompt for GPT
   prompt = f"""

I am porting the ddfun library from Fortran to C++.
The ddfun library contains subroutines for performing arithmetic operations on double-double numbers.

Given this example ddadd function:
---------------
{ddadd_func}
---------------
I have the following python file that generates test data for the ddadd function:
----------------
{gen_code}
----------------
Then I have the fortran unit test for the ddadd function:
----------------
{test_code}
----------------
I need a similar fortran file that reads the generated test data and checks the results for the {function} function, which is defined as follows:

{function['body']}

Can you please generate a test fortran code following the same format as in the ddadd case? Please avoid any extra text or markdown formatting because your reply will go directly into a file to be compiled.
Do not even include the markdown formatting around the code.
"""   

   # Call the GPT API
   completion = client.chat.completions.create(
      model="gpt-4o",  # Use the desired model
      store=True,
      messages=[
         {"role": "user", "content": prompt}
      ]
   )

   # Extract and return the generated test code
   return completion.choices[0].message.content


# Main script
if __name__ == "__main__":
   # input fortran file contains the original version of David Bailey's ddfun code
   # we want to extract the functions from this file
   input_file = "/home/jchilders/git/precisionStudies/ddfun/ddfun-v03/fortran/ddfuna.f90"
   functions = extract_functions_with_bodies(input_file)

   # path for python generator file
   gen_path = "python_generators/gen_%s.py"
   test_path = "tests/test_%s.f90"

   # loop over each function and generate C++ code
   for func in functions:
      if func["name"] in functions_to_extract:
         print(f"Generating C++ code for {func['name']}")
         # generate C++ code using GPT
         python_code = generate_gen_with_gpt(func)

         # write the generated C++ code to the output file
         with open(gen_path % func["name"], "w") as file:
            file.write(python_code)
         
         # generate C++ code using GPT
         test_code = generate_test_with_gpt(func)

         # write the generated C++ code to the output file
         with open(test_path % func["name"], "w") as file:
            file.write(test_code)
        

        
         

         