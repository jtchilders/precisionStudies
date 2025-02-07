#ifndef TEST_DDCSSHR_H
#define TEST_DDCSSHR_H

#include "ddmath.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

void unittest_ddcsshr(const std::string &filename) {
   std::ifstream infile(filename, std::ios::binary);
   if (!infile) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      return;
   }

   std::cout << "Running tests for ddcsshr from: " << filename << std::endl;

   ddouble a,expected_x,expected_y;
   int total_tests = 0;
   int passed_tests = 0;
   const double tolerance = 29;

   while (infile.read(reinterpret_cast<char *>(&a.hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&a.lo), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_x.hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_x.lo), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_y.hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_y.lo), sizeof(double))) {
             
      total_tests++;

      // Perform ddcsshr operation
      ddouble x,y;
      ddcsshr(a,x,y);

      // Use scale difference for comparison
      int scale_diff_x = calculate_scale_difference(x, expected_x);
      int scale_diff_y = calculate_scale_difference(y, expected_y);
      bool test_passed = (scale_diff_x >= tolerance || scale_diff_x == 0) &&
                         (scale_diff_y >= tolerance || scale_diff_y == 0);

      if (test_passed) {
         passed_tests++;
      } else {
         std::cout << "Test Failed: " << std::setprecision(16)
                   << "a: [" << a << "] "
                   << "x: [" << x << "] "
                   << "y: [" << y << "] "
                   << "expected_x: [" << expected_x << "] "
                   << "expected_y: [" << expected_y << "] "
                   << "scale_diff: [" << scale_diff_x << ", " << scale_diff_y << "]"
                   << std::endl;
      }
   }

   infile.close();

   std::cout << "Tests completed: " << passed_tests << " passed out of " << total_tests << std::endl;
}

#endif // TEST_DDCSSHR_H