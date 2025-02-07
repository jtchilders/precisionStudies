#ifndef TEST_DDEXP_H
#define TEST_DDEXP_H

#include "ddmath.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

void unittest_ddexp(const std::string &filename) {
   std::ifstream infile(filename, std::ios::binary);
   if (!infile) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      return;
   }

   std::cout << "Running tests for ddexp from: " << filename << std::endl;

   double hi_a, lo_a, expected_hi, expected_lo;
   int total_tests = 0;
   int passed_tests = 0;
   const int tolerance = 28;

   while (infile.read(reinterpret_cast<char *>(&hi_a), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&lo_a), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_lo), sizeof(double))) {

      total_tests++;

      // Perform ddexp operation
      ddouble a{hi_a, lo_a};
      ddouble expected{expected_hi, expected_lo};

      try {
         ddouble result = ddexp(a);

         // Use scale difference for comparison
         int scale_diff = calculate_scale_difference(result, expected);
         bool test_passed = (scale_diff >= tolerance or scale_diff == 0);

         if (test_passed) {
            passed_tests++;
         } else {
            std::cout << "Test Failed: " << std::setprecision(16)
                      << "input: [" << a.hi << ", " << a.lo << "] "
                      << "result: [" << result.hi << ", " << result.lo << "] "
                      << "expected: [" << expected.hi << ", " << expected.lo << "]"
                      << "error: [" << std::abs(result.hi - expected.hi) << ", " << std::abs(result.lo - expected.lo) << "]"
                      << "scale difference: [" << scale_diff << "]\n";
         }
      } catch (const std::overflow_error &e) {
         if (expected.hi == std::numeric_limits<double>::infinity() && expected.lo == 0.0) {
            passed_tests++;
         } else {
            std::cout << "Test Failed: Overflow with input: [" << a.hi << ", " << a.lo << "] "
                      << "expected: [Inf, 0.0]" << std::endl;
         }
      }
   }

   infile.close();

   std::cout << "Tests completed: " << passed_tests << " passed out of " << total_tests << std::endl;
}

#endif // TEST_DDEXP_H