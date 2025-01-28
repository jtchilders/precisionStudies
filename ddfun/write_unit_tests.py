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
   # "ddabs",
   # "ddacosh",
   # "ddadd",
   # "ddasinh",
   # "ddatanh",
   # "dddiv",
   # "dddivd",
   # "ddexp",
   # "ddlog",
   # "ddmul",
   # "ddmuld",
   # "ddmuldd",
   # "ddneg",
   # "ddnint",
   # "ddnpwr",
   "ddpolyr"
   # "ddpower",
   # "ddsqrt",
   # "ddsub",
]

py_gen_ddadd_inputs = """# gen_ddadd.py
import random
import struct
import math
from mpmath import mp


# Set high precision for calculations
mp.dps = 50

def generate_double_double(base_value):
   \"""
   Generate a double-double representation (hi, lo) of a high-precision value.
   \"""
   hi = float(base_value)
   lo = float(base_value - mp.mpf(hi))  # Compute the remainder
   return hi, lo


def calculate_exponent_difference(a,b):
   \"""
   Calculate the difference of scale between two double values
   \"""
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
   \"""
   Verify that hi and lo accurately represent the original value.
   \"""
   reconstructed = mp.mpf(hi) + mp.mpf(lo)
   error = abs(original - reconstructed)
   if error == 0.0:
      return True, error
   scale_diff = calculate_exponent_difference(original, error)
   return scale_diff > tolerance, error

def generate_ddadd_test_case():
   \"""
   Generate a test case for ddadd using PI as the base value.
   \"""
   # Generate two random scaling factors [0,1) in high precision
   # May need to adjust random number inputs for different test functions
   a = mp.rand()
   b = mp.rand()

   # Compute expected result
   result = mp.fadd(a,b)

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
         f"Verification of ddadd failed for tolerance {tolerance}:\\n"
         f"  a: {a} (hi={a_hi}, lo={a_lo}, error={error_a})\\n"
         f"  b: {b} (hi={b_hi}, lo={b_lo}, error={error_b})\\n"
         f"  result: {result} (hi={result_hi}, lo={result_lo}, error={error_result})"
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
               f"{case['dda'][0]:.16e} {case['dda'][1]:.16e} "
               f"{case['ddb'][0]:.16e} {case['ddb'][1]:.16e} "
               f"{case['expected'][0]:.16e} {case['expected'][1]:.16e}\n"
         )
   print(f"Test cases successfully written to {filename}")

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
   print(f"Test cases successfully written to {filename}")

if __name__ == "__main__":
   print("Generating test cases for ddadd...")
   test_cases = generate_test_cases(10)
   write_test_cases_to_text_file("data/ddadd_test_cases.txt", test_cases)
   write_test_cases_to_binary("data/ddadd_test_cases.bin", test_cases)

"""


fortran_ddadd_unittest = """! test_ddadd.f90
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

cpp_ddadd_unittest = """#ifndef TEST_DDADD_H
#define TEST_DDADD_H

#include "ddmath.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

void unittest_ddadd(const std::string &filename) {
   std::ifstream infile(filename, std::ios::binary);
   if (!infile) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      return;
   }

   std::cout << "Running tests for ddadd from: " << filename << std::endl;

   double hi_a, lo_a, hi_b, lo_b, expected_hi, expected_lo;
   int total_tests = 0;
   int passed_tests = 0;
   const double tolerance = 1.0e-30;

   while (infile.read(reinterpret_cast<char *>(&hi_a), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&lo_a), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&hi_b), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&lo_b), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_lo), sizeof(double))) {

      total_tests++;

      // Perform ddadd operation
      ddouble a{hi_a, lo_a}, b{hi_b, lo_b};
      ddouble expected{expected_hi, expected_lo};
      ddouble result = ddadd(a, b);

      // Compare results
      bool test_passed = (std::abs(result.hi - expected.hi) < tolerance) &&
                         (std::abs(result.lo - expected.lo) < tolerance);

      if (test_passed) {
         passed_tests++;
      } else {
         std::cout << "Test Failed: " << std::setprecision(16)
                   << "inputs: [" << a.hi << ", " << a.lo << "] + [" << b.hi << ", " << b.lo << "] "
                   << "result: [" << result.hi << ", " << result.lo << "] "
                   << "expected: [" << expected.hi << ", " << expected.lo << "] " 
                   << "error: [" << std::abs(result.hi - expected.hi) << ", " << std::abs(result.lo - expected.lo) << "]" << std::endl;
      }
   }

   infile.close();

   std::cout << "Tests completed: " << passed_tests << " passed out of " << total_tests << std::endl;
}

#endif // TEST_DDADD_H
"""

fortran_ddadd_func = """subroutine ddadd (dda, ddb, ddc)

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

cpp_ddadd_func = """ddouble ddadd(const ddouble& dda, const ddouble& ddb) {
    ddouble ddc;
    double e, t1, t2;

    // Compute dda + ddb using Knuth's trick.
    t1 = dda.hi + ddb.hi;
    e = t1 - dda.hi;
    t2 = ((ddb.hi - e) + (dda.hi - (t1 - e))) + dda.lo + ddb.lo;

    // The result is t1 + t2, after normalization.
    ddc.hi = t1 + t2;
    ddc.lo = t2 - (ddc.hi - t1);
    return ddc;
}
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


def generate_cpp_func(function, outfilename):
   # Create a prompt for GPT
   prompt = f"""Based on the following Fortran {function['type']}:
---
{function['body']}
---

Generate a corresponding C++ function that performs the same operation as the Fortran function.
Assume that a ddouble struct containing a double hi and a double lo have already been defined.
The function should return a ddouble struct containing the result of the operation.
Cases where the Fortran calls the ddabrt due to errors, the C++ function should print an error and return an empty ddouble struct.
The function should be named {function['name']} and can return a ddouble struct containing the result of the operation.

Please only reply with the function code using no formatting because your reply will go directly into a file to be compiled. Do not even include the markdown 
formatting around the code. No includes are needed, just the function code. Do not redefine the ddouble struct.
Any functions that appear in the Fortran code, you can assume that they are already defined in C++ and follow a similar function signature except that results are always returned, not passed as arguments.
"""   

   # Call the GPT API
   completion = client.chat.completions.create(
      model="gpt-4o",  # Use the desired model
      store=True,
      messages=[
         {"role": "user", "content": prompt}
      ]
   )

   # Extract the generated test code
   code = completion.choices[0].message.content

   # Write the generated C++ code to the output file
   with open(outfilename, "w") as file:
      file.write(code)
   
   return code


def gen_python_generator_with_gpt(function, outfilename):
   # Create a prompt for GPT
   prompt = f"""I am porting the ddfun library from Fortran to C++.
The ddfun library contains subroutines for performing arithmetic operations on double-double numbers.

I need to generate input test data for each Fortran and C++ function.

Given this example ddadd function:
---------------
{fortran_ddadd_func}
---------------
I have the following python file that generates test data that is used as input for the ddadd function in both C++ and Fortran:
----------------
{py_gen_ddadd_inputs}
----------------

Please create a similar python file for the {function['name']} function, which is defined as follows:

{function['body']}

In the python file, you may need to adjust the range of input random numbers generated to be appropriate for the {function['name']} function.

Only reply with the python code using no formatting because your reply will go directly into a file to be compiled. Do not even include the markdown formatting around the code.
"""

   # Call the GPT API
   completion = client.chat.completions.create(
      model="gpt-4o",  # Use the desired model
      store=True,
      messages=[
         {"role": "user", "content": prompt}
      ]
   )

   # Extract the generated test code
   code = completion.choices[0].message.content

   # Write the generated C++ code to the output file
   with open(outfilename, "w") as file:
      file.write(code)

def generate_fortran_unittest(function, outfilename):
   # Create a prompt for GPT
   prompt = f"""Given this unit test in fortran for the ddadd function:
---------------
{fortran_ddadd_unittest}
---------------

Please create a similar unit test in fortran for the {function['name']} function, which is defined as follows:

{function['body']}

Only reply with the fortran code using no formatting because your reply will go directly into a file to be compiled. Do not even include the markdown formatting around the code.
"""

   # Call the GPT API
   completion = client.chat.completions.create(
      model="gpt-4o",  # Use the desired model
      store=True,
      messages=[
         {"role": "user", "content": prompt}
      ]
   )

   # Extract the generated test code
   code = completion.choices[0].message.content

   # Write the generated C++ code to the output file
   with open(outfilename, "w") as file:
      file.write(code)

def generate_cpp_unittest(function, cpp_func, outfilename):
   # Create a prompt for GPT
   prompt = f"""Given this unit test in C++ for the ddadd function:
---------------
{cpp_ddadd_unittest}
---------------

Please create a similar unit test in C++ for this function:

{cpp_func}

Only reply with the c++ code using no formatting because your reply will go directly into a file to be compiled. Do not even include the markdown formatting around the code.
"""

   # Call the GPT API
   completion = client.chat.completions.create(
      model="gpt-4o",  # Use the desired model
      store=True,
      messages=[
         {"role": "user", "content": prompt}
      ]
   )

   # Extract the generated test code
   code = completion.choices[0].message.content

   # Write the generated C++ code to the output file
   with open(outfilename, "w") as file:
      file.write(code)


def update_support_files():
   # C++ outputs
   # for C++ build I need to update the CMakeLists.txt and the main.cpp
   cpp_base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),"cppPort")
   cpp_src_dir = os.path.join(cpp_base_dir,"src")
   cpp_unit_test_dir = os.path.join(cpp_base_dir,"unit_testing")

   # get list of all test_*.h files in cpp_unit_test_dir
   test_files = [f for f in os.listdir(cpp_unit_test_dir) if f.startswith("test_")]
   test_files = sorted(test_files)

   # update CMakeLists.txt
   cmake_code = "set(CPP_UNIT_TEST_DEPS\n"
   for test_file in test_files:
      cmake_code += f"   {test_file}\n"
   cmake_code += """)\n

# compile main into a binary, depends on the unit tests
add_executable(unit_tests main.cpp ${CPP_UNIT_TEST_DEPS})
target_link_libraries(unit_tests ddmath)
"""

   with open(os.path.join(cpp_unit_test_dir,"CMakeLists.txt"), "w") as file:
      file.write(cmake_code)
   print("Updated CMakeLists.txt")

   # update main.cpp
   main_code = "#include <string>\n"
   # add unit test headers
   for test_file in test_files:
      main_code += f"#include \"{test_file}\"\n"
   main_code += "\n"
   
   # begin main
   main_code += """int main(int argc, char* argv[]) {

   // user passes path to input files as first parameter on command line
   // extract the path
   std::string path = std::string(argv[1]);

"""
   # add unit test code
   for test_file in test_files:
      # extract function name from "test_<function_name>.h"
      func_name = test_file[5:-2]
      main_code += f"   unittest_{func_name}(path + \"/{func_name}_test_cases.bin\");\n"

   # end main
   main_code += """
   return 0;
}
"""

   with open(os.path.join(cpp_unit_test_dir,"main.cpp"), "w") as file:
      file.write(main_code)
   print("Updated main.cpp")


   # fortran outputs
   ftn_base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),"ddfun-v03")
   ftn_src_dir = os.path.join(ftn_base_dir,"fortran")
   fnt_unit_test_dir = os.path.join(ftn_base_dir,"unit_testing")

   # get list of all test_*.f90 files in fnt_unit_test_dir
   test_files = [f for f in os.listdir(fnt_unit_test_dir) if f.startswith("test_")]
   test_files = sorted(test_files)

   # create cmake file
   cmake_code = """
# make list of test files
set(TEST_SOURCES\n"""

   for test_file in test_files:
      cmake_code += f"   {test_file}\n"
   cmake_code += """)\n

# link into an executable
add_executable(unit_tests main.f90 ${TEST_SOURCES})
# target_link_libraries(unit_tests unit_test_lib)
target_link_libraries(unit_tests ddfun)
target_include_directories(unit_tests PRIVATE ${CMAKE_BINARY_DIR}/fortran)
"""

   with open(os.path.join(fnt_unit_test_dir,"CMakeLists.txt"), "w") as file:
      file.write(cmake_code)
   print("updated CMakeLists.txt")
   
   # update main.f90
   main_code = "program main\n"
   # add "use test_<function_name>" lines
   for test_file in test_files:
      # remove ".f90" from test_file
      test_name = test_file.split(".")[0]
      main_code += f"   use {test_name}\n"

   main_code += """   implicit none

   ! variables for input file path
   character(len=100) :: input_files_path

   ! read path to input files from command line
   call get_command_argument(1, input_files_path)

   write(*,*) "Starting unit tests for DDFUN library..."

   ! Call unit tests with respective input files
"""

   for test_file in test_files:
      # remove ".f90" from test_file
      test_name = test_file.split(".")[0]
      func_name = test_name.replace("test_","")
      main_code += f"   call unittest_{func_name}(trim(input_files_path) // \"/{func_name}_test_cases.bin\")\n"

   main_code += """
   write(*,*) "All tests completed."

end program main
"""

   with open(os.path.join(fnt_unit_test_dir,"main.f90"), "w") as file:
      file.write(main_code)
   print("updated main.f90")



# Main script
if __name__ == "__main__":
   # input fortran file contains the original version of David Bailey's ddfun code
   # we want to extract the functions from this file
   input_file = "/home/jchilders/git/precisionStudies/ddfun/ddfun-v03/fortran/ddfuna.f90"
   functions = extract_functions_with_bodies(input_file)

   # loop over each function and generate C++ code
   for func in functions:
      if func["name"] in functions_to_extract:
         
         # generate C++ code using GPT
         print(f"Generating C++ code for {func['name']}")
         cpp_func = generate_cpp_func(func,f"cppPort/{func['name']}.h")

         # generate python code using GPT
         print(f"Generating python generator for {func['name']}")
         gen_python_generator_with_gpt(func,f"ddTestInputs/gen_{func['name']}.py")

         # generate fortran unit test using GPT
         print(f"Generating fortran unit test for {func['name']}")
         generate_fortran_unittest(func,f"ddfun-v03/unit_testing/test_{func['name']}.f90")

         # generate C++ unit test using GPT
         print(f"Generating C++ unit test for {func['name']}")
         generate_cpp_unittest(func,cpp_func,f"cppPort/unit_testing/test_{func['name']}.h")

   # print("Updating support files")
   update_support_files()
        

        
         

         