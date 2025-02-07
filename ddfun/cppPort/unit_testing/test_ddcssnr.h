#ifndef TEST_DDCSSNR_H
#define TEST_DDCSSNR_H

#include "ddmath.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

void unittest_ddcssnr(const std::string &filename) {
   std::ifstream infile(filename, std::ios::binary);
   if (!infile) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      return;
   }

   std::cout << "Running tests for ddcssnr from: " << filename << std::endl;

   double hi_a, lo_a, expected_x_hi, expected_x_lo, expected_y_hi, expected_y_lo;
   int total_tests = 0;
   int passed_tests = 0;
   const double tolerance = 28;

   while (infile.read(reinterpret_cast<char *>(&hi_a), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&lo_a), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_x_hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_x_lo), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_y_hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_y_lo), sizeof(double))) {

      total_tests++;

      // Perform ddcssnr operation
      ddouble a{hi_a, lo_a};
      ddouble expected_x{expected_x_hi, expected_x_lo}, expected_y{expected_y_hi, expected_y_lo};
      ddouble x,y;
      ddcssnr(a,x,y);

      // Use scale difference for comparison
      int scale_diff_cos = calculate_scale_difference(x, expected_x);
      int scale_diff_sin = calculate_scale_difference(y, expected_y);


      // Compare results with expected values within the tolerance
      bool test_passed = ((scale_diff_cos >= tolerance or scale_diff_cos == 0) &&
                          (scale_diff_sin >= tolerance or scale_diff_sin == 0));

      if (test_passed) {
         passed_tests++;
      } else {
         std::cout << "Test Failed: " << std::setprecision(16)
                   << "input: [" << a.hi << ", " << a.lo << "] "
                   << "result cos: [" << x.hi << ", " << x.lo << "] "
                   << "expected cos: [" << expected_x.hi << ", " << expected_x.lo << "] "
                   << "result sin: [" << y.hi << ", " << y.lo << "] "
                   << "expected sin: [" << expected_y.hi << ", " << expected_y.lo << "] "
                   << "scale difference cos: " << scale_diff_cos << " "
                   << "scale difference sin: " << scale_diff_sin
                   << std::endl;
      }
   }

   infile.close();

   std::cout << "Tests completed: " << passed_tests << " passed out of " << total_tests << std::endl;
}

#endif // TEST_DDCSSNR_H