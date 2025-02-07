#ifndef TEST_DDASINH_H
#define TEST_DDASINH_H

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "Kokkos_Core.hpp"
#include "ddouble.h"
#include "ddmath.h"

void unittest_ddasinh(const std::string &filename) {
   std::ifstream infile(filename, std::ios::binary);
   if (!infile) {
      std::cerr << "Error: Cannot open file " << filename << std::endl;
      return;
   }

   std::cout << "Running tests for ddasinh from: " << filename << std::endl;

   // Each test: a.hi, a.lo, expected.hi, expected.lo
   double hi_a, lo_a, expected_hi, expected_lo;
   int total_tests = 0;
   int passed_tests = 0;
   const double tolerance = 29;

   // Host vectors to store the input data
   std::vector<double> hi_a_vec, lo_a_vec;
   std::vector<double> exp_hi_vec, exp_lo_vec;

   // Read until we cannot read 4 doubles per test
   while (infile.read(reinterpret_cast<char *>(&hi_a), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&lo_a), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_hi), sizeof(double)) &&
          infile.read(reinterpret_cast<char *>(&expected_lo), sizeof(double))) 
   {
      total_tests++;
      hi_a_vec.push_back(hi_a);
      lo_a_vec.push_back(lo_a);
      exp_hi_vec.push_back(expected_hi);
      exp_lo_vec.push_back(expected_lo);
   }
   infile.close();

   // Create Kokkos Views for input and expected output
   Kokkos::View<ddouble *> a("a", total_tests);
   auto ha = Kokkos::create_mirror_view(a);

   Kokkos::View<ddouble *> expected("expected", total_tests);
   auto hexp = Kokkos::create_mirror_view(expected);

   // Fill the host mirrors
   for (int i = 0; i < total_tests; i++) {
      ha(i)   = ddouble(hi_a_vec[i], lo_a_vec[i]);
      hexp(i) = ddouble(exp_hi_vec[i], exp_lo_vec[i]);
   }

   // Copy data to device
   Kokkos::deep_copy(a, ha);
   Kokkos::deep_copy(expected, hexp);

   // View for results
   Kokkos::View<ddouble *> result("result", total_tests);
   auto hresult = Kokkos::create_mirror_view(result);

   // Parallel kernel: result(i) = ddasinh(a(i))
   Kokkos::parallel_for("test_ddasinh", total_tests, KOKKOS_LAMBDA(const int i) {
      result(i) = ddasinh(a(i));
   });

   // Copy results back to host
   Kokkos::deep_copy(hresult, result);

   // Compare results on host
   for (int i = 0; i < total_tests; i++) {
      ddouble val_a       = ha(i);
      ddouble val_expected= hexp(i);
      ddouble val_result  = hresult(i);

      // Use scale difference for comparison
      int scale_diff = calculate_scale_difference(val_result, val_expected);
      bool test_passed = (scale_diff >= tolerance or scale_diff == 0);

      if (test_passed) {
         passed_tests++;
      } else {
         std::cout << std::setprecision(16)
                   << "Test Failed: "
                   << "input: " << val_a << " "
                   << "result: " << val_result << " "
                   << "expected: " << val_expected << " "
                   << "error: ["
                   << std::abs(val_result.hi - val_expected.hi) << ", "
                   << std::abs(val_result.lo - val_expected.lo)
                   << "]"
                   << "scale difference: [" << scale_diff << "]\n";
      }
   }

   std::cout << "Tests completed: " << passed_tests
             << " passed out of " << total_tests << std::endl;
}

#endif // TEST_DDASINH_H
