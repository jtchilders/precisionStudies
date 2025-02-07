#ifndef TEST_DDCSQRT_H
#define TEST_DDCSQRT_H

#include "ddmath.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

void unittest_ddcsqrt(const std::string &filename) {
   std::ifstream infile(filename, std::ios::binary);
   if (!infile) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      return;
   }

   std::cout << "Running tests for ddcsqrt from: " << filename << std::endl;

   double hi_real, lo_real, hi_imag, lo_imag;
   double expected_hi_real, expected_lo_real, expected_hi_imag, expected_lo_imag;
   int total_tests = 0;
   int passed_tests = 0;
   const double tolerance = 30;

   while (infile.read(reinterpret_cast<char *>(&hi_real), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&lo_real), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&hi_imag), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&lo_imag), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_hi_real), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_lo_real), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_hi_imag), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_lo_imag), sizeof(double))) {

      total_tests++;

      // Perform ddcsqrt operation
      ddouble real{hi_real, lo_real}, imag{hi_imag, lo_imag};
      ddcomplex a{real, imag};
      ddcomplex expected{{expected_hi_real, expected_lo_real}, {expected_hi_imag, expected_lo_imag}};
      ddcomplex result = ddcsqrt(a);

      // Use scale difference for comparison
      int scale_diff_real = calculate_scale_difference(result.real, expected.real);
      int scale_diff_imag = calculate_scale_difference(result.imag, expected.imag);
      bool test_passed = (scale_diff_real >= tolerance or scale_diff_real == 0) &&
                         (scale_diff_imag >= tolerance or scale_diff_imag == 0);

      if (test_passed) {
         passed_tests++;
      } else {
         std::cout << "Test Failed: "
                   << " inputs: [" << a.real.hi << ", " << a.real.lo << "] + [" << a.imag.hi << ", " << a.imag.lo << "]i "
                   << " result: [" << result.real.hi << ", " << result.real.lo << "] + [" << result.imag.hi << ", " << result.imag.lo << "]i "
                   << " expected: [" << expected.real.hi << ", " << expected.real.lo << "] + [" << expected.imag.hi << ", " << expected.imag.lo << "]i "
                   << " error: [" << std::abs(result.real.hi - expected.real.hi) << ", " << std::abs(result.real.lo - expected.real.lo) << "]real "
                   << " [" << std::abs(result.imag.hi - expected.imag.hi) << ", " << std::abs(result.imag.lo - expected.imag.lo) << "]imag "
                   << " scale difference: [" << scale_diff_real << "]real [" << scale_diff_imag << "]imag "
                   << std::endl;
      }
   }

   infile.close();

   std::cout << "Tests completed: " << passed_tests << " passed out of " << total_tests << std::endl;
}

#endif // TEST_DDCSQRT_H