#ifndef TEST_DDCPWR_H
#define TEST_DDCPWR_H

#include "ddmath.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

void unittest_ddcpwr(const std::string &filename) {
   std::ifstream infile(filename, std::ios::binary);
   if (!infile) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      return;
   }

   std::cout << "Running tests for ddcpwr from: " << filename << std::endl;

   double hi_real, lo_real, hi_imag, lo_imag, expected_hi, expected_lo;
   int n;
   int total_tests = 0;
   int passed_tests = 0;
   const double tolerance = 30;

   while (infile.read(reinterpret_cast<char *>(&hi_real), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&lo_real), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&hi_imag), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&lo_imag), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&n), sizeof(int)) &&
          infile.read(reinterpret_cast<char *>(&expected_hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_lo), sizeof(double))) {

      total_tests++;

      // Prepare inputs
      ddouble real_part{hi_real, lo_real};
      ddouble imag_part{hi_imag, lo_imag};
      ddcomplex a{real_part, imag_part};
      ddouble expected{expected_hi, expected_lo};

      // Perform ddcpwr operation
      ddouble result = ddcpwr(a, n);

      // Use scale difference for comparison
      int scale_diff = calculate_scale_difference(result, expected);
      bool test_passed = (scale_diff >= tolerance or scale_diff == 0);

      if (test_passed) {
         passed_tests++;
      } else {
         std::cout << "Test Failed: " << std::setprecision(16)
                   << "input: [" << a.real.hi << ", " << a.real.lo << "] + [" << a.imag.hi << ", " << a.imag.lo << "], n: " << n << " "
                   << "result: [" << result.hi << ", " << result.lo << "] "
                   << "expected: [" << expected.hi << ", " << expected.lo << "] "
                   << "error: [" << std::abs(result.hi - expected.hi) << ", " << std::abs(result.lo - expected.lo) << "]"
                   << "scale difference: [" << scale_diff << "]" << std::endl;
      }
   }

   infile.close();

   std::cout << "Tests completed: " << passed_tests << " passed out of " << total_tests << std::endl;
}

#endif // TEST_DDCPWR_H