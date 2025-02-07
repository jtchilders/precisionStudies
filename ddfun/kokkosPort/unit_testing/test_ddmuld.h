#ifndef TEST_DDMULD_H
#define TEST_DDMULD_H

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <limits>

#include "Kokkos_Core.hpp"
#include "ddouble.h"
#include "ddmath.h"

// Kokkos-based unittest for ddmuld
void unittest_ddmuld(const std::string &filename) {
   std::ifstream infile(filename, std::ios::binary);
   if (!infile) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      return;
   }

   std::cout << "Running tests for ddmuld from: " << filename << std::endl;

   double hi_a, lo_a, b, expected_hi, expected_lo;
   int total_tests = 0;
   int passed_tests = 0;
   const int tolerance = 30;

   std::vector<double> hi_a_vec, lo_a_vec, b_vec;
   std::vector<double> exp_hi_vec, exp_lo_vec;

   while (infile.read(reinterpret_cast<char *>(&hi_a), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&lo_a), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&b), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_lo), sizeof(double))) 
   {
      total_tests++;
      hi_a_vec.push_back(hi_a);
      lo_a_vec.push_back(lo_a);
      b_vec.push_back(b);
      exp_hi_vec.push_back(expected_hi);
      exp_lo_vec.push_back(expected_lo);
   }
   infile.close();

   // Create Kokkos Views for input (a, b) and expected
   Kokkos::View<ddouble *> da("a", total_tests);
   Kokkos::View<double *> db("b", total_tests);
   auto ha = Kokkos::create_mirror_view(da);
   auto hb = Kokkos::create_mirror_view(db);

   Kokkos::View<ddouble *> dexpected("expected", total_tests);
   auto hexpected = Kokkos::create_mirror_view(dexpected);

   // Fill host mirrors
   for (int i = 0; i < total_tests; i++) {
      ha(i) = ddouble(hi_a_vec[i], lo_a_vec[i]);
      hb(i) = b_vec[i];
      hexpected(i) = ddouble(exp_hi_vec[i], exp_lo_vec[i]);
   }

   // Copy to device
   Kokkos::deep_copy(da, ha);
   Kokkos::deep_copy(db, hb);
   Kokkos::deep_copy(dexpected, hexpected);

   // Prepare a result array
   Kokkos::View<ddouble *> dresult("result", total_tests);
   auto hresult = Kokkos::create_mirror_view(dresult);

   // Perform ddmuld in parallel
   Kokkos::parallel_for("test_ddmuld", total_tests, KOKKOS_LAMBDA(const int i) {
      dresult(i) = ddmuld(da(i),db(i));
   });

   // Copy back to host
   Kokkos::deep_copy(hresult, dresult);

   // Check on host
   for (int i = 0; i < total_tests; i++) {
      ddouble val_a = ha(i);
      double val_b = hb(i);
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
                   << "inputs: " << val_a << " * " << val_b << " "
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

#endif // TEST_DDMULD_H
