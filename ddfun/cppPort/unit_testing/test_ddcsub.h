#ifndef TEST_DDCSUB_H
#define TEST_DDCSUB_H

#include "ddmath.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

void unittest_ddcsub(const std::string &filename) {
   std::ifstream infile(filename, std::ios::binary);
   if (!infile) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      return;
   }

   std::cout << "Running tests for ddcsub from: " << filename << std::endl;

   ddcomplex a, b, expected;
   int total_tests = 0;
   int passed_tests = 0;
   const double tolerance = 30;

   while (infile.read(reinterpret_cast<char *>(&a.real.hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&a.real.lo), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&a.imag.hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&a.imag.lo), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&b.real.hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&b.real.lo), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&b.imag.hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&b.imag.lo), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected.real.hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected.real.lo), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected.imag.hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected.imag.lo), sizeof(double))) {


      total_tests++;

      // Perform ddcsub operation
      ddcomplex result = a - b;

      // Use scale difference for comparison
      int scale_diff_real = calculate_scale_difference(result.real, expected.real);
      int scale_diff_imag = calculate_scale_difference(result.imag, expected.imag);
      bool test_passed = (scale_diff_real >= tolerance or scale_diff_real == 0) and 
                         (scale_diff_imag >= tolerance or scale_diff_imag == 0);


      if (test_passed) {
         passed_tests++;
      } else {
         std::cout << "Test Failed: " << std::setprecision(16)
                   << "inputs: [" << a.real << ", " << a.imag << "] / [" << b.real << ", " << b.imag << "]"
                   << "result: [" << result.real << ", " << result.imag << "]"
                   << "expected: [" << expected.real << ", " << expected.imag << "]"
                   << "scale difference: [" << scale_diff_real << ", " << scale_diff_imag << "]"
                   << std::endl;
      }
   }

   infile.close();

   std::cout << "Tests completed: " << passed_tests << " passed out of " << total_tests << std::endl;
}

#endif // TEST_DDCSUB_H