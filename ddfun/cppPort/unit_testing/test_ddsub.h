#ifndef TEST_DDSUB_H
#define TEST_DDSUB_H

#include "ddmath.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

void unittest_ddsub(const std::string &filename) {
   std::ifstream infile(filename, std::ios::binary);
   if (!infile) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      return;
   }

   std::cout << "Running tests for ddsub from: " << filename << std::endl;

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

      // Perform ddsub operation
      ddouble a{hi_a, lo_a}, b{hi_b, lo_b};
      ddouble expected{expected_hi, expected_lo};
      ddouble result = a - b;

      // Compare results
      bool test_passed = (std::abs(result.hi - expected.hi) < tolerance) &&
                        (std::abs(result.lo - expected.lo) < tolerance);

      if (test_passed) {
         passed_tests++;
      } else {
         std::cout << "Test Failed: " << std::setprecision(16)
                     << "inputs: [" << a.hi << ", " << a.lo << "] - [" << b.hi << ", " << b.lo << "] "
                     << "result: [" << result.hi << ", " << result.lo << "] "
                     << "expected: [" << expected.hi << ", " << expected.lo << "]"
                     << "error: [" << std::abs(result.hi - expected.hi) << ", " << std::abs(result.lo - expected.lo) << "]" << std::endl;
      }
   }

   infile.close();

   std::cout << "Tests completed: " << passed_tests << " passed out of " << total_tests << std::endl;
}

#endif // TEST_DDSUB_H
