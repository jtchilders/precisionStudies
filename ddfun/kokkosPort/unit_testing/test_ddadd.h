#ifndef TEST_DDADD_H
#define TEST_DDADD_H

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <limits>

#include "Kokkos_Core.hpp"
#include "ddouble.h"
#include "ddmath.h"

// Kokkos-based unittest for ddadd
void unittest_ddadd(const std::string &filename) {
   std::ifstream infile(filename, std::ios::binary);
   if (!infile) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      return;
   }

   std::cout << "Running tests for ddadd from: " << filename << std::endl;

   double hi_a, lo_a, hi_b, lo_b, expected_hi, expected_lo;
   int total_tests = 0;
   int passed_tests = 0;
   const int tolerance = 30;

   std::vector<double> hi_a_vec, lo_a_vec, hi_b_vec, lo_b_vec;
   std::vector<double> exp_hi_vec, exp_lo_vec;

   while (infile.read(reinterpret_cast<char *>(&hi_a), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&lo_a), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&hi_b), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&lo_b), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_lo), sizeof(double))) 
   {
      total_tests++;
      hi_a_vec.push_back(hi_a);
      lo_a_vec.push_back(lo_a);
      hi_b_vec.push_back(hi_b);
      lo_b_vec.push_back(lo_b);
      exp_hi_vec.push_back(expected_hi);
      exp_lo_vec.push_back(expected_lo);
   }
   infile.close();

   // Create Kokkos Views for input (a, b) and expected
   Kokkos::View<ddouble *> a("a", total_tests);
   Kokkos::View<ddouble *> b("b", total_tests);
   auto ha = Kokkos::create_mirror_view(a);
   auto hb = Kokkos::create_mirror_view(b);

   Kokkos::View<ddouble *> expected("expected", total_tests);
   auto hexpected = Kokkos::create_mirror_view(expected);

   // Fill host mirrors
   for (int i = 0; i < total_tests; i++) {
      ha(i) = ddouble(hi_a_vec[i], lo_a_vec[i]);
      hb(i) = ddouble(hi_b_vec[i], lo_b_vec[i]);
      hexpected(i) = ddouble(exp_hi_vec[i], exp_lo_vec[i]);
   }

   // Copy to device
   Kokkos::deep_copy(a, ha);
   Kokkos::deep_copy(b, hb);
   Kokkos::deep_copy(expected, hexpected);

   // Prepare a result array
   Kokkos::View<ddouble *> result("result", total_tests);
   auto hresult = Kokkos::create_mirror_view(result);

   // Perform ddadd in parallel
   Kokkos::parallel_for("test_ddadd", total_tests, KOKKOS_LAMBDA(const int i) {
      result(i) = a(i) + b(i);
   });

   // Copy back to host
   Kokkos::deep_copy(hresult, result);

   // Check on host
   for (int i = 0; i < total_tests; i++) {
      ddouble val_a = ha(i);
      ddouble val_b = hb(i);
      ddouble val_expected = hexpected(i);
      ddouble val_result = hresult(i);

      // Use scale difference for comparison
      int scale_diff = calculate_scale_difference(val_result, val_expected);
      bool test_passed = (scale_diff >= tolerance or scale_diff == 0);

      if (test_passed) {
         passed_tests++;
      } else {
         std::cout << std::setprecision(16)
                   << "Test Failed: "
                   << "inputs: " << val_a << " + " << val_b << " "
                   << "result: " << val_result << " "
                   << "expected: " << val_expected << " "
                   << "error: ["
                   << std::abs(val_result.hi - val_expected.hi) << ", "
                   << std::abs(val_result.lo - val_expected.lo) << "] "
                   << "scale difference: [" << scale_diff << "]\n";
      }
   }

   std::cout << "Tests completed: " << passed_tests
             << " passed out of " << total_tests << std::endl;
}

#endif // TEST_DDADD_H
