#ifndef TEST_DDLOG_H
#define TEST_DDLOG_H

#include "ddmath.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

void unittest_ddlog(const std::string &filename) {
   std::ifstream infile(filename, std::ios::binary);
   if (!infile) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      return;
   }

   std::cout << "Running tests for ddlog from: " << filename << std::endl;

   double hi_a, lo_a, expected_hi, expected_lo;
   int total_tests = 0;
   int passed_tests = 0;
   const double tolerance = 27;

   while (infile.read(reinterpret_cast<char *>(&hi_a), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&lo_a), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_lo), sizeof(double))) {

      total_tests++;

      // Perform ddlog operation
      ddouble a{hi_a, lo_a};
      ddouble expected{expected_hi, expected_lo};
      ddouble result = a.log();

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
   }

   infile.close();

   std::cout << "Tests completed: " << passed_tests << " passed out of " << total_tests << std::endl;
}

#endif // TEST_DDLOG_H