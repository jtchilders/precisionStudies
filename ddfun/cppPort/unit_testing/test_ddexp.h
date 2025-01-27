#ifndef TEST_DDEXP_H
#define TEST_DDEXP_H

#include "ddmath.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

int calculate_scale_difference(const ddouble &result, const ddouble &expected){
   double error_hi_abs = std::abs(result.hi - expected.hi);

   if (error_hi_abs > 0.0){
      double error_hi_exponent = std::log10(error_hi_abs);
      double expected_hi_exponent = std::log10(std::abs(expected.hi));
      int scale_difference = static_cast<int>(std::abs(error_hi_exponent - expected_hi_exponent));
      return scale_difference;
   }

   double error_lo_abs = std::abs(result.lo - expected.lo);
   if(error_lo_abs > 0.0){
      double error_lo_exponent = std::log10(error_lo_abs);
      double expected_hi_exponent = std::log10(std::abs(expected.hi));
      int scale_difference = static_cast<int>(std::abs(error_lo_exponent - expected_hi_exponent));
      return scale_difference;
   }

   return 0;
}

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
   const int tolerance = 30;

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

         // Compare results
         bool test_passed = calculate_scale_difference(result, expected) >= tolerance;

         if (test_passed) {
            passed_tests++;
         } else {
            std::cout << "Test Failed: " << std::setprecision(16)
                      << "input: [" << a.hi << ", " << a.lo << "] "
                      << "result: [" << result.hi << ", " << result.lo << "] "
                      << "expected: [" << expected.hi << ", " << expected.lo << "]"
                      << "error: [" << std::abs(result.hi - expected.hi) << ", " << std::abs(result.lo - expected.lo) << "]" << std::endl;
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