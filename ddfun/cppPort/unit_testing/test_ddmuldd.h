#ifndef TEST_DDMULDD_H
#define TEST_DDMULDD_H

#include "ddmath.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

void unittest_ddmuldd(const std::string &filename) {
   std::ifstream infile(filename, std::ios::binary);
   if (!infile) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      return;
   }

   std::cout << "Running tests for ddmuldd from: " << filename << std::endl;

   double da, db, expected_hi, expected_lo;
   int total_tests = 0;
   int passed_tests = 0;
   const double tolerance = 30;

   while (infile.read(reinterpret_cast<char*>(&da), sizeof(double)) &&
          infile.read(reinterpret_cast<char*>(&db), sizeof(double)) &&
          infile.read(reinterpret_cast<char*>(&expected_hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char*>(&expected_lo), sizeof(double))) {

      total_tests++;

      // Perform ddmuldd operation
      ddouble expected{expected_hi, expected_lo};
      ddouble result = ddmuldd(da, db);

      // Use scale difference for comparison
      int scale_diff = calculate_scale_difference(result, expected);
      bool test_passed = (scale_diff >= tolerance or scale_diff == 0);

      if (test_passed) {
         passed_tests++;
      } else {
         std::cout << "Test Failed: " << std::setprecision(16)
                   << "inputs: [" << da << ", " << db << "] "
                   << "result: [" << result.hi << ", " << result.lo << "] "
                   << "expected: [" << expected.hi << ", " << expected.lo << "]"
                   << "error: [" << std::abs(result.hi - expected.hi) << ", " << std::abs(result.lo - expected.lo) << "]"
                   << "scale difference: [" << scale_diff << "]\n";
      }
   }

   infile.close();

   std::cout << "Tests completed: " << passed_tests << " passed out of " << total_tests << std::endl;
}

#endif // TEST_DDMULDD_H