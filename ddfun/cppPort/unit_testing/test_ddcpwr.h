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

   ddcomplex a,expected;
   int n,placeholder;
   int total_tests = 0;
   int passed_tests = 0;
   const double tolerance = 28;

   while (infile.read(reinterpret_cast<char *>(&a.real.hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&a.real.lo), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&a.imag.hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&a.imag.lo), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&n), sizeof(int)) &&
          infile.read(reinterpret_cast<char *>(&placeholder), sizeof(int)) &&
          infile.read(reinterpret_cast<char *>(&expected.real.hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected.real.lo), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected.imag.hi), sizeof(double)) &&  
          infile.read(reinterpret_cast<char *>(&expected.imag.lo), sizeof(double))) {
           

      total_tests++;

      // Perform ddcpwr operation
      ddcomplex result = ddcpwr(a, n);

      // Use scale difference for comparison
      int scale_diff_real = calculate_scale_difference(result.real, expected.real);
      int scale_diff_imag = calculate_scale_difference(result.imag, expected.imag);
      bool test_passed = (scale_diff_real >= tolerance || scale_diff_real == 0) && 
                         (scale_diff_imag >= tolerance || scale_diff_imag == 0);

      if (test_passed) {
         passed_tests++;
      } else {
         std::cout << "Test Failed: " << std::setprecision(16)
                   << "input: [" << a.real << "] + [" << a.imag << "], n: " << n << " "
                   << "result: [" << result.real << ", " << result.imag << "] "
                   << "expected: [" << expected.real << ", " << expected.imag << "] "
                   << "scale difference: [" << scale_diff_real << ", " << scale_diff_imag << "]" << std::endl;
      }
   }

   infile.close();

   std::cout << "Tests completed: " << passed_tests << " passed out of " << total_tests << std::endl;
}

#endif // TEST_DDCPWR_H