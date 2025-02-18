#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include "ddouble.h"
#include "ddmath.h"
#include "ddcomplex.h"
#include <fstream>
#include <cstring>
#include <cmath>
#include <iostream>

const std::string INPUT_FILES_DIR("/home/jchilders/git/precisionStudies/ddfun/ddTestInputs/");

// Fortran routines must be exported with ISO_C_BINDING.
// Uncomment and update the prototypes as needed.
extern "C" {
   // math operators
   void ddadd_fortran(const double a[2], const double b[2], double c[2]);
   void ddsub_fortran(const double a[2], const double b[2], double c[2]);
   void ddmul_fortran(const double a[2], const double b[2], double c[2]);
   void ddmuld_fortran(const double a[2], const double& b, double c[2]);
   void ddmuldd_fortran(const double& a, const double& b, double c[2]);
   void dddiv_fortran(const double a[2], const double b[2], double c[2]);
   void dddivd_fortran(const double a[2], const double& b, double c[2]);

   // For math functions:
   void ddexp_fortran(const double a[2], double b[2]);
   void ddlog_fortran(const double a[2], double b[2]);
   void ddnint_fortran(const double a[2], double b[2]);
   void ddabs_fortran(const double a[2], double b[2]);
   void ddsqrt_fortran(const double a[2], double b[2]);
   void ddnpwr_fortran(const double a[2], const int& b, double c[2]);
   void ddpower_fortran(const double a[2], const double b[2], double c[2]);

   // For trig functions
   void ddacosh_fortran(const double a[2], double b[2]);
   void ddasinh_fortran(const double a[2], double b[2]);
   void ddatanh_fortran(const double a[2], double b[2]);
   void ddcsshr_fortran(const double a[2], double b[2], double c[2]);
   void ddcssnr_fortran(const double a[2], double b[2], double c[2]);

   // complex math operators
   void ddcadd_fortran(const double a[4], const double b[4], double c[4]);
   void ddcsub_fortran(const double a[4], const double b[4], double c[4]);
   void ddcmul_fortran(const double a[4], const double b[4], double c[4]);
   void ddcdiv_fortran(const double a[4], const double b[4], double c[4]);

   // complex math functions
   void ddcpwr_fortran(const double a[4], const int& b, double c[4]);
   void ddcsqrt_fortran(const double a[4], double b[4]);
}

using namespace ddfun;

// Helper function that computes the scale difference between result and expected.
// We compute the number of orders of magnitude between the error and the expected value.
int calculate_scale_difference(const ddouble &result, const ddouble &expected) {
   double error_hi_abs = std::fabs(result.hi - expected.hi);
   if (error_hi_abs > 0.0) {
      double error_hi_exponent = std::log10(error_hi_abs);
      double expected_hi_exponent = std::log10(std::fabs(expected.hi));
      int scale_difference = static_cast<int>(std::fabs(error_hi_exponent - expected_hi_exponent));
      return scale_difference;
   }
   double error_lo_abs = std::fabs(result.lo - expected.lo);
   if (error_lo_abs > 0.0) {
      double error_lo_exponent = std::log10(error_lo_abs);
      double expected_hi_exponent = std::log10(std::fabs(expected.hi));
      int scale_difference = static_cast<int>(std::fabs(error_lo_exponent - expected_hi_exponent));
      return scale_difference;
   }
   return 0;
}

//
// ddadd Test Case
//
TEST_CASE("ddadd operator test", "[ddfun]") {
   // Binary file format for ddadd.bin:
   //   2 mp inputs (each: 2 doubles) followed by 1 mp output (2 doubles).
   struct ddaddRecord {
      double a_hi;
      double a_lo;
      double b_hi;
      double b_lo;
      double exp_hi;
      double exp_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddadd.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddaddRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddouble a(rec.a_hi, rec.a_lo);
      ddouble b(rec.b_hi, rec.b_lo);
      ddouble expected(rec.exp_hi, rec.exp_lo);
      
      // Compute the result using the C++ operator.
      ddouble cppResult = a + b;
      
      int scaleDiff = calculate_scale_difference(cppResult, expected);
      INFO("C++ ddadd scale difference: " << scaleDiff << ", expected vs result:\n hi " << expected.hi << " <> " << cppResult.hi << "; lo\n " << expected.lo << " <> " << cppResult.lo);
      // We expect at least 20 digits of precision.
      bool test_passed = (scaleDiff >= 20) || (scaleDiff == 0);
      REQUIRE(test_passed);
      
      // Optionally, call the Fortran version.
      double a_arr[2] = { rec.a_hi, rec.a_lo };
      double b_arr[2] = { rec.b_hi, rec.b_lo };
      double f_arr[2];
      // Uncomment if Fortran ddadd is available:
      ddadd_fortran(a_arr, b_arr, f_arr);
      ddouble fortranResult(f_arr[0], f_arr[1]);
      int scaleDiffFortran = calculate_scale_difference(fortranResult, expected);
      INFO("Fortran ddadd scale difference: " << scaleDiffFortran << ", expected vs result:\n hi " << expected.hi << " <> " << fortranResult.hi << "; lo\n " << expected.lo << " <> " << fortranResult.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffFortran >= 20) || (scaleDiffFortran == 0);
      REQUIRE(test_passed);
   }
}

//
// ddsub Test Case
//
TEST_CASE("ddsub operator test", "[ddfun]") {
   // Binary file format for ddsub.bin:
   //   2 mp inputs (each: 2 doubles) followed by 1 mp output (2 doubles).
   struct ddsubRecord {
      double a_hi;
      double a_lo;
      double b_hi;
      double b_lo;
      double exp_hi;
      double exp_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddsub.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddsubRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddouble a(rec.a_hi, rec.a_lo);
      ddouble b(rec.b_hi, rec.b_lo);
      ddouble expected(rec.exp_hi, rec.exp_lo);
      
      // Compute the result using the C++ operator.
      ddouble cppResult = a - b;
      
      int scaleDiff = calculate_scale_difference(cppResult, expected);
      INFO("C++ ddsub scale difference: " << scaleDiff << ", expected vs result:\n hi " << expected.hi << " <> " << cppResult.hi << "; lo\n " << expected.lo << " <> " << cppResult.lo);
      // We expect at least 20 digits of precision.
      bool test_passed = (scaleDiff >= 20) || (scaleDiff == 0);
      REQUIRE(test_passed);

      // call Fortran ddsub
      double a_arr[2] = { rec.a_hi, rec.a_lo };
      double b_arr[2] = { rec.b_hi, rec.b_lo };
      double f_arr[2];
      // Uncomment if Fortran ddsub is available:
      ddsub_fortran(a_arr, b_arr, f_arr);
      ddouble fortranResult(f_arr[0], f_arr[1]);
      int scaleDiffFortran = calculate_scale_difference(fortranResult, expected);
      INFO("Fortran ddsub scale difference: " << scaleDiffFortran << ", expected vs result:\n hi " << expected.hi << " <> " << fortranResult.hi << "; lo\n " << expected.lo << " <> " << fortranResult.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffFortran >= 20) || (scaleDiffFortran == 0);
      REQUIRE(test_passed);
   }
}

//
// ddmul Test Case
//
TEST_CASE("ddmul operator test", "[ddfun]") {
   // Binary file format for ddmul.bin:
   //   2 mp inputs (each: 2 doubles) followed by 1 mp output (2 doubles).
   struct ddmulRecord {
      double a_hi;
      double a_lo;
      double b_hi;
      double b_lo;
      double exp_hi;
      double exp_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddmul.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddmulRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddouble a(rec.a_hi, rec.a_lo);
      ddouble b(rec.b_hi, rec.b_lo);
      ddouble expected(rec.exp_hi, rec.exp_lo);
      
      // Compute the result using the C++ operator.
      ddouble cppResult = a * b;
      
      int scaleDiff = calculate_scale_difference(cppResult, expected);
      INFO("C++ ddmul scale difference: " << scaleDiff << ", expected vs result:\n hi " << expected.hi << " <> " << cppResult.hi << "; lo\n " << expected.lo << " <> " << cppResult.lo);
      // We expect at least 20 digits of precision.
      bool test_passed = (scaleDiff >= 20) || (scaleDiff == 0);
      REQUIRE(test_passed);

      // call Fortran ddmul
      double a_arr[2] = { rec.a_hi, rec.a_lo };
      double b_arr[2] = { rec.b_hi, rec.b_lo };
      double f_arr[2];
      // Uncomment if Fortran ddmul is available:
      ddmul_fortran(a_arr, b_arr, f_arr);
      ddouble fortranResult(f_arr[0], f_arr[1]);
      int scaleDiffFortran = calculate_scale_difference(fortranResult, expected);
      INFO("Fortran ddmul scale difference: " << scaleDiffFortran << ", expected vs result:\n hi " << expected.hi << " <> " << fortranResult.hi << "; lo\n " << expected.lo << " <> " << fortranResult.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffFortran >= 20) || (scaleDiffFortran == 0);
      REQUIRE(test_passed);
   }
}

//
// ddmuld Test Case
//
TEST_CASE("ddmuld operator test", "[ddfun]") {
   // Binary file format for ddmuld.bin:
   //   1 mp and 1 double input followed by 1 mp output (2 doubles).
   struct ddmuldRecord {
      double a_hi;
      double a_lo;
      double b;
      double exp_hi;
      double exp_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddmuld.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddmuldRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddouble a(rec.a_hi, rec.a_lo);
      double b(rec.b);
      ddouble expected(rec.exp_hi, rec.exp_lo);
      
      // Compute the result using the C++ operator.
      ddouble cppResult = ddmuld(a, b);
      
      int scaleDiff = calculate_scale_difference(cppResult, expected);
      INFO("C++ ddmuld scale difference: " << scaleDiff << ", expected vs result:\n hi " << expected.hi << " <> " << cppResult.hi << "; lo\n " << expected.lo << " <> " << cppResult.lo);
      // We expect at least 20 digits of precision.
      bool test_passed = (scaleDiff >= 20) || (scaleDiff == 0);
      REQUIRE(test_passed);

      // call Fortran ddmuld
      double a_arr[2] = { rec.a_hi, rec.a_lo };
      double f_arr[2];
      // Uncomment if Fortran ddmuld is available:
      ddmuld_fortran(a_arr, b, f_arr);
      ddouble fortranResult(f_arr[0], f_arr[1]);
      int scaleDiffFortran = calculate_scale_difference(fortranResult, expected);
      INFO("Fortran ddmuld scale difference: " << scaleDiffFortran << ", expected vs result:\n hi " << expected.hi << " <> " << fortranResult.hi << "; lo\n " << expected.lo << " <> " << fortranResult.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffFortran >= 20) || (scaleDiffFortran == 0);
      REQUIRE(test_passed);
   }
}

//
// ddmuldd Test Case
//
TEST_CASE("ddmuldd operator test", "[ddfun]") {
   // Binary file format for ddmuldd.bin:
   //   2 doubles followed by 1 mp output (2 doubles).
   struct ddmulddRecord {
      double a;
      double b;
      double exp_hi;
      double exp_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddmuldd.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddmulddRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      double a(rec.a);
      double b(rec.b);
      ddouble expected(rec.exp_hi, rec.exp_lo);
      
      // Compute the result using the C++ operator.
      ddouble cppResult = ddmuldd(a, b);
      
      int scaleDiff = calculate_scale_difference(cppResult, expected);
      INFO("C++ ddmuldd scale difference: " << scaleDiff << ", expected vs result:\n hi " << expected.hi << " <> " << cppResult.hi << "; lo\n " << expected.lo << " <> " << cppResult.lo);
      // We expect at least 20 digits of precision.
      bool test_passed = (scaleDiff >= 20) || (scaleDiff == 0);

      // call Fortran ddmuldd
      double f_arr[2];
      // Uncomment if Fortran ddmuldd is available:
      ddmuldd_fortran(a, b, f_arr);
      ddouble fortranResult(f_arr[0], f_arr[1]);
      int scaleDiffFortran = calculate_scale_difference(fortranResult, expected);
      INFO("Fortran ddmuldd scale difference: " << scaleDiffFortran << ", expected vs result:\n hi " << expected.hi << " <> " << fortranResult.hi << "; lo\n " << expected.lo << " <> " << fortranResult.lo);
      // We expect at least 20 digits of precision.
      bool fortran_test_passed = (scaleDiffFortran >= 20) || (scaleDiffFortran == 0);


      REQUIRE(test_passed);
      REQUIRE(fortran_test_passed);
   }
}

//
// dddiv Test Case
// 
TEST_CASE("dddiv operator test", "[ddfun]") {
   // Binary file format for dddiv.bin:
   //   2 mp inputs (each: 2 doubles) followed by 1 mp output (2 doubles).
   struct dddivRecord {
      double a_hi;
      double a_lo;
      double b_hi;
      double b_lo;
      double exp_hi;
      double exp_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/dddiv.bin", std::ios::binary);
   REQUIRE(infile.good());
   dddivRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddouble a(rec.a_hi, rec.a_lo);
      ddouble b(rec.b_hi, rec.b_lo);
      ddouble expected(rec.exp_hi, rec.exp_lo);
      
      // Compute the result using the C++ operator.
      ddouble cppResult = a / b;
      
      int scaleDiff = calculate_scale_difference(cppResult, expected);
      INFO("C++ dddiv scale difference: " << scaleDiff << ", expected vs result:\n hi " << expected.hi << " <> " << cppResult.hi << "; lo\n " << expected.lo << " <> " << cppResult.lo);
      // We expect at least 20 digits of precision.
      bool test_passed = (scaleDiff >= 20) || (scaleDiff == 0);
      REQUIRE(test_passed);

      // call Fortran dddiv
      double a_arr[2] = { rec.a_hi, rec.a_lo };
      double b_arr[2] = { rec.b_hi, rec.b_lo };
      double f_arr[2];
      // Uncomment if Fortran dddiv is available:
      dddiv_fortran(a_arr, b_arr, f_arr);
      ddouble fortranResult(f_arr[0], f_arr[1]);
      int scaleDiffFortran = calculate_scale_difference(fortranResult, expected);
      INFO("Fortran dddiv scale difference: " << scaleDiffFortran << ", expected vs result:\n hi " << expected.hi << " <> " << fortranResult.hi << "; lo\n " << expected.lo << " <> " << fortranResult.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffFortran >= 20) || (scaleDiffFortran == 0);
      REQUIRE(test_passed);
   }
}

//
// dddivd Test Case
//
TEST_CASE("dddivd operator test", "[ddfun]") {
   // Binary file format for dddivd.bin:
   //   1 mp input (2 doubles), 1 double input, followed by 1 mp output (2 doubles).
   struct dddivdRecord {
      double a_hi;
      double a_lo;
      double b;
      double exp_hi;
      double exp_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/dddivd.bin", std::ios::binary);
   REQUIRE(infile.good());
   dddivdRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddouble a(rec.a_hi, rec.a_lo);
      double b(rec.b);
      ddouble expected(rec.exp_hi, rec.exp_lo);
      
      // Compute the result using the C++ operator.
      ddouble cppResult = dddivd(a, b);
      
      int scaleDiff = calculate_scale_difference(cppResult, expected);
      INFO("C++ dddivd scale difference: " << scaleDiff << ", expected vs result:\n hi " << expected.hi << " <> " << cppResult.hi << "; lo\n " << expected.lo << " <> " << cppResult.lo);
      // We expect at least 20 digits of precision.
      bool test_passed = (scaleDiff >= 20) || (scaleDiff == 0);
      REQUIRE(test_passed);

      // call Fortran dddivd
      double a_arr[2] = { rec.a_hi, rec.a_lo };
      double f_arr[2];
      // Uncomment if Fortran dddivd is available:
      dddivd_fortran(a_arr, b, f_arr);
      ddouble fortranResult(f_arr[0], f_arr[1]);
      int scaleDiffFortran = calculate_scale_difference(fortranResult, expected);
      INFO("Fortran dddivd scale difference: " << scaleDiffFortran << ", expected vs result:\n hi " << expected.hi << " <> " << fortranResult.hi << "; lo\n " << expected.lo << " <> " << fortranResult.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffFortran >= 20) || (scaleDiffFortran == 0);
      REQUIRE(test_passed);
   }
}

//
// ddexp Test Case
//
TEST_CASE("ddexp function test", "[ddfun]") {
   // Binary file format for ddexp.bin:
   //   1 mp input (2 doubles) followed by 1 mp output (2 doubles).
   struct ddexpRecord {
      double a_hi;
      double a_lo;
      double exp_hi;
      double exp_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddexp.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddexpRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddouble a(rec.a_hi, rec.a_lo);
      ddouble expected(rec.exp_hi, rec.exp_lo);
      
      // Compute the result using the C++ function.
      ddouble cppResult = a.exp();
      
      int scaleDiff = calculate_scale_difference(cppResult, expected);
      INFO("C++ ddexp scale difference: " << scaleDiff << "; expected vs result:\n hi " << expected.hi << " <> " << cppResult.hi << "; lo\n " << expected.lo << " <> " << cppResult.lo);
      bool test_passed = (scaleDiff >= 20) || (scaleDiff == 0);
      REQUIRE(test_passed);
      
      // Call the Fortran version.
      double a_arr[2] = { rec.a_hi, rec.a_lo };
      double f_arr[2];
      ddexp_fortran(a_arr, f_arr);
      ddouble fortranResult(f_arr[0], f_arr[1]);
      int scaleDiffFortran = calculate_scale_difference(fortranResult, expected);
      INFO("Fortran ddexp scale difference: " << scaleDiffFortran << "; expected vs result:\n hi " << expected.hi << " <> " << fortranResult.hi << "; lo\n " << expected.lo << " <> " << fortranResult.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffFortran >= 20) || (scaleDiffFortran == 0);
      REQUIRE(test_passed);
   }
}

//
// ddlog Test Case
//
TEST_CASE("ddlog function test", "[ddfun]") {
   // Binary file format for ddlog.bin:
   //   1 mp input (2 doubles) followed by 1 mp output (2 doubles).
   struct ddlogRecord {
      double a_hi;
      double a_lo;
      double exp_hi;
      double exp_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddlog.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddlogRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddouble a(rec.a_hi, rec.a_lo);
      ddouble expected(rec.exp_hi, rec.exp_lo);
      
      // Compute the result using the C++ function.
      ddouble cppResult = a.log();
      
      int scaleDiff = calculate_scale_difference(cppResult, expected);
      INFO("C++ ddlog scale difference: " << scaleDiff << "; expected vs result:\n hi " << expected.hi << " <> " << cppResult.hi << "; lo\n " << expected.lo << " <> " << cppResult.lo);
      bool test_passed = (scaleDiff >= 20) || (scaleDiff == 0);
      REQUIRE(test_passed);
      
      // Call the Fortran version.
      double a_arr[2] = { rec.a_hi, rec.a_lo };
      double f_arr[2];
      ddlog_fortran(a_arr, f_arr);
      ddouble fortranResult(f_arr[0], f_arr[1]);
      int scaleDiffFortran = calculate_scale_difference(fortranResult, expected);
      INFO("Fortran ddlog scale difference: " << scaleDiffFortran << "; expected vs result:\n hi " << expected.hi << " <> " << fortranResult.hi << "; lo\n " << expected.lo << " <> " << fortranResult.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffFortran >= 20) || (scaleDiffFortran == 0);
      REQUIRE(test_passed);
   }
}

//
// ddnint Test Case
//
TEST_CASE("ddnint function test", "[ddfun]") {
   // Binary file format for ddnint.bin:
   //   1 mp input (2 doubles) followed by 1 mp output (2 doubles).
   struct ddnintRecord {
      double a_hi;
      double a_lo;
      double exp_hi;
      double exp_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddnint.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddnintRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddouble a(rec.a_hi, rec.a_lo);
      ddouble expected(rec.exp_hi, rec.exp_lo);
      
      // Compute the result using the C++ function.
      ddouble cppResult = a.nint();
      
      int scaleDiff = calculate_scale_difference(cppResult, expected);
      INFO("C++ ddnint scale difference: " << scaleDiff << "; expected vs result:\n hi " << expected.hi << " <> " << cppResult.hi << "; lo\n " << expected.lo << " <> " << cppResult.lo);
      bool test_passed = (scaleDiff >= 20) || (scaleDiff == 0);
      REQUIRE(test_passed);
      
      // Call the Fortran version.
      double a_arr[2] = { rec.a_hi, rec.a_lo };
      double f_arr[2];
      ddnint_fortran(a_arr, f_arr);
      ddouble fortranResult(f_arr[0], f_arr[1]);
      int scaleDiffFortran = calculate_scale_difference(fortranResult, expected);
      INFO("Fortran ddnint scale difference: " << scaleDiffFortran << "; expected vs result:\n hi " << expected.hi << " <> " << fortranResult.hi << "; lo\n " << expected.lo << " <> " << fortranResult.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffFortran >= 20) || (scaleDiffFortran == 0);
      REQUIRE(test_passed);
   }
}

//
// ddabs Test Case
//
TEST_CASE("ddabs function test", "[ddfun]") {
   // Binary file format for ddabs.bin:
   //   1 mp input (2 doubles) followed by 1 mp output (2 doubles).
   struct ddabsRecord {
      double a_hi;
      double a_lo;
      double exp_hi;
      double exp_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddabs.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddabsRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddouble a(rec.a_hi, rec.a_lo);
      ddouble expected(rec.exp_hi, rec.exp_lo);
      
      // Compute the result using the C++ function.
      ddouble cppResult = a.abs();
      
      int scaleDiff = calculate_scale_difference(cppResult, expected);
      INFO("C++ ddabs scale difference: " << scaleDiff << "; expected vs result:\n hi " << expected.hi << " <> " << cppResult.hi << "; lo\n " << expected.lo << " <> " << cppResult.lo);
      bool test_passed = (scaleDiff >= 20) || (scaleDiff == 0);
      REQUIRE(test_passed);
      
      // Call the Fortran version.
      double a_arr[2] = { rec.a_hi, rec.a_lo };
      double f_arr[2];
      ddabs_fortran(a_arr, f_arr);
      ddouble fortranResult(f_arr[0], f_arr[1]);
      int scaleDiffFortran = calculate_scale_difference(fortranResult, expected);
      INFO("Fortran ddabs scale difference: " << scaleDiffFortran << "; expected vs result:\n hi " << expected.hi << " <> " << fortranResult.hi << "; lo\n " << expected.lo << " <> " << fortranResult.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffFortran >= 20) || (scaleDiffFortran == 0);
      REQUIRE(test_passed);
   }
}

//
// ddsqrt Test Case
//
TEST_CASE("ddsqrt function test", "[ddfun]") {
   // Binary file format for ddsqrt.bin:
   //   1 mp input (2 doubles) followed by 1 mp output (2 doubles).
   struct ddsqrtRecord {
      double a_hi;
      double a_lo;
      double exp_hi;
      double exp_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddsqrt.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddsqrtRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddouble a(rec.a_hi, rec.a_lo);
      ddouble expected(rec.exp_hi, rec.exp_lo);
      
      // Compute the result using the C++ function.
      ddouble cppResult = a.sqrt();
      
      int scaleDiff = calculate_scale_difference(cppResult, expected);
      INFO("C++ ddsqrt scale difference: " << scaleDiff << "; expected vs result:\n hi " << expected.hi << " <> " << cppResult.hi << "; lo\n " << expected.lo << " <> " << cppResult.lo);
      bool test_passed = (scaleDiff >= 20) || (scaleDiff == 0);
      REQUIRE(test_passed);
      
      // Call the Fortran version.
      double a_arr[2] = { rec.a_hi, rec.a_lo };
      double f_arr[2];
      ddsqrt_fortran(a_arr, f_arr);
      ddouble fortranResult(f_arr[0], f_arr[1]);
      int scaleDiffFortran = calculate_scale_difference(fortranResult, expected);
      INFO("Fortran ddsqrt scale difference: " << scaleDiffFortran << "; expected vs result:\n hi " << expected.hi << " <> " << fortranResult.hi << "; lo\n " << expected.lo << " <> " << fortranResult.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffFortran >= 20) || (scaleDiffFortran == 0);
      REQUIRE(test_passed);
   }
}

//
// ddnpwr Test Case
//
TEST_CASE("ddnpwr function test", "[ddfun]") {
   // Binary file format for ddnpwr.bin:
   //   1 mp input (2 doubles) followed by 1 mp output (2 doubles).
   struct ddnpwrRecord {
      double a_hi;
      double a_lo;
      int n;
      int ph;
      double exp_hi;
      double exp_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddnpwr.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddnpwrRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddouble a(rec.a_hi, rec.a_lo);
      int n = rec.n;
      ddouble expected(rec.exp_hi, rec.exp_lo);
      
      // Compute the result using the C++ function.
      ddouble cppResult = a.npwr(n);
      
      int scaleDiff = calculate_scale_difference(cppResult, expected);
      INFO("C++ ddnpwr scale difference: " << scaleDiff << "; expected vs result:\n hi " << expected.hi << " <> " << cppResult.hi << "; lo\n " << expected.lo << " <> " << cppResult.lo);
      bool test_passed = (scaleDiff >= 20) || (scaleDiff == 0);
      REQUIRE(test_passed);
      
      // Call the Fortran version.
      double a_arr[2] = { rec.a_hi, rec.a_lo };
      double f_arr[2];
      ddnpwr_fortran(a_arr, n, f_arr);
      ddouble fortranResult(f_arr[0], f_arr[1]);
      int scaleDiffFortran = calculate_scale_difference(fortranResult, expected);
      INFO("Fortran ddnpwr scale difference: " << scaleDiffFortran << "; expected vs result:\n hi " << expected.hi << " <> " << fortranResult.hi << "; lo\n " << expected.lo << " <> " << fortranResult.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffFortran >= 20) || (scaleDiffFortran == 0);
      REQUIRE(test_passed);
   }
}

//
// ddpower Test Case
//
TEST_CASE("ddpower function test", "[ddfun]") {
   // Binary file format for ddpower.bin:
   //   2 mp input (4 doubles) followed by 1 mp output (2 doubles).
   struct ddpowerRecord {
      double a_hi;
      double a_lo;
      double n_hi;
      double n_lo;
      double exp_hi;
      double exp_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddpower.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddpowerRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddouble a(rec.a_hi, rec.a_lo);
      ddouble n(rec.n_hi, rec.n_lo);
      ddouble expected(rec.exp_hi, rec.exp_lo);
      
      // Compute the result using the C++ function.
      ddouble cppResult = a.power(n);
      
      int scaleDiff = calculate_scale_difference(cppResult, expected);
      INFO("C++ ddpower scale difference: " << scaleDiff << "; expected vs result:\n hi " << expected.hi << " <> " << cppResult.hi << "; lo\n " << expected.lo << " <> " << cppResult.lo);
      bool test_passed = (scaleDiff >= 20) || (scaleDiff == 0);
      REQUIRE(test_passed);
      
      // Call the Fortran version.
      double a_arr[2] = { rec.a_hi, rec.a_lo };
      double n_arr[2] = { rec.n_hi, rec.n_lo };
      double f_arr[2];
      ddpower_fortran(a_arr, n_arr, f_arr);
      ddouble fortranResult(f_arr[0], f_arr[1]);
      int scaleDiffFortran = calculate_scale_difference(fortranResult, expected);
      INFO("Fortran ddpower scale difference: " << scaleDiffFortran << "; expected vs result:\n hi " << expected.hi << " <> " << fortranResult.hi << "; lo\n " << expected.lo << " <> " << fortranResult.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffFortran >= 20) || (scaleDiffFortran == 0);
      REQUIRE(test_passed);
   }
}

// 
// ddacosh Test Case
//
TEST_CASE("ddacosh function test", "[ddfun]") {
   // Binary file format for ddacosh.bin:
   //   1 mp input (2 doubles) followed by 1 mp output (2 doubles).
   struct ddacoshRecord {
      double a_hi;
      double a_lo;
      double exp_hi;
      double exp_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddacosh.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddacoshRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddouble a(rec.a_hi, rec.a_lo);
      ddouble expected(rec.exp_hi, rec.exp_lo);
      
      // Compute the result using the C++ function.
      ddouble cppResult = a.acosh();
      
      int scaleDiff = calculate_scale_difference(cppResult, expected);
      INFO("C++ ddacosh scale difference: " << scaleDiff << "; expected vs result:\n hi " << expected.hi << " <> " << cppResult.hi << "; lo\n " << expected.lo << " <> " << cppResult.lo);
      bool test_passed = (scaleDiff >= 20) || (scaleDiff == 0);
      REQUIRE(test_passed);

      // Call the Fortran version.
      double a_arr[2] = { rec.a_hi, rec.a_lo };
      double f_arr[2];
      ddacosh_fortran(a_arr, f_arr);
      ddouble fortranResult(f_arr[0], f_arr[1]);
      int scaleDiffFortran = calculate_scale_difference(fortranResult, expected);
      INFO("Fortran ddacosh scale difference: " << scaleDiffFortran << "; expected vs result:\n hi " << expected.hi << " <> " << fortranResult.hi << "; lo\n " << expected.lo << " <> " << fortranResult.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffFortran >= 20) || (scaleDiffFortran == 0);
      REQUIRE(test_passed);
   }
}

//
// ddasinh Test Case
//
TEST_CASE("ddasinh function test", "[ddfun]") {
   // Binary file format for ddasinh.bin:
   //   1 mp input (2 doubles) followed by 1 mp output (2 doubles).
   struct ddasinhRecord {
      double a_hi;
      double a_lo;
      double exp_hi;
      double exp_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddasinh.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddasinhRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddouble a(rec.a_hi, rec.a_lo);
      ddouble expected(rec.exp_hi, rec.exp_lo);
      
      // Compute the result using the C++ function.
      ddouble cppResult = a.asinh();
      
      int scaleDiff = calculate_scale_difference(cppResult, expected);
      INFO("C++ ddasinh scale difference: " << scaleDiff << "; expected vs result:\n hi " << expected.hi << " <> " << cppResult.hi << "; lo\n " << expected.lo << " <> " << cppResult.lo);
      bool test_passed = (scaleDiff >= 20) || (scaleDiff == 0);
      REQUIRE(test_passed);
      
      // Call the Fortran version.
      double a_arr[2] = { rec.a_hi, rec.a_lo };
      double f_arr[2];
      ddasinh_fortran(a_arr, f_arr);
      ddouble fortranResult(f_arr[0], f_arr[1]);
      int scaleDiffFortran = calculate_scale_difference(fortranResult, expected);
      INFO("Fortran ddasinh scale difference: " << scaleDiffFortran << "; expected vs result:\n hi " << expected.hi << " <> " << fortranResult.hi << "; lo\n " << expected.lo << " <> " << fortranResult.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffFortran >= 20) || (scaleDiffFortran == 0);
      REQUIRE(test_passed);
   }
}

// 
// ddatanh Test Case
//
TEST_CASE("ddatanh function test", "[ddfun]") {
   // Binary file format for ddatanh.bin:
   //   1 mp input (2 doubles) followed by 1 mp output (2 doubles).
   struct ddatanhRecord {
      double a_hi;
      double a_lo;
      double exp_hi;
      double exp_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddatanh.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddatanhRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddouble a(rec.a_hi, rec.a_lo);
      ddouble expected(rec.exp_hi, rec.exp_lo);
      
      // Compute the result using the C++ function.
      ddouble cppResult = a.atanh();
      
      int scaleDiff = calculate_scale_difference(cppResult, expected);
      INFO("C++ ddatanh scale difference: " << scaleDiff << "; expected vs result:\n hi " << expected.hi << " <> " << cppResult.hi << "; lo\n " << expected.lo << " <> " << cppResult.lo);
      bool test_passed = (scaleDiff >= 20) || (scaleDiff == 0);
      REQUIRE(test_passed);

      // Call the Fortran version.
      double a_arr[2] = { rec.a_hi, rec.a_lo };
      double f_arr[2];
      ddatanh_fortran(a_arr, f_arr);
      ddouble fortranResult(f_arr[0], f_arr[1]);
      int scaleDiffFortran = calculate_scale_difference(fortranResult, expected);
      INFO("Fortran ddatanh scale difference: " << scaleDiffFortran << "; expected vs result:\n hi " << expected.hi << " <> " << fortranResult.hi << "; lo\n " << expected.lo << " <> " << fortranResult.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffFortran >= 20) || (scaleDiffFortran == 0);
      REQUIRE(test_passed);
   }
}

//
// ddcsshr Test Case
//
TEST_CASE("ddcsshr function test", "[ddfun]") {
   // Binary file format for ddcsshr.bin:
   //   1 mp input (2 doubles) followed by 2 mp outputs (4 doubles).
   struct ddcsshrRecord {
      double a_hi;
      double a_lo;
      double x_hi;
      double x_lo;
      double y_hi;
      double y_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddcsshr.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddcsshrRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddouble a(rec.a_hi, rec.a_lo);
      ddouble x(rec.x_hi, rec.x_lo);
      ddouble y(rec.y_hi, rec.y_lo);
      ddouble cpp_x, cpp_y;
      
      // Compute the result using the C++ function.
      ddcsshr(a, cpp_x, cpp_y);
      
      int scaleDiffX = calculate_scale_difference(cpp_x, x);
      int scaleDiffY = calculate_scale_difference(cpp_y, y);
      INFO("C++ ddcsshr scale difference x: " << scaleDiffX << "; expected vs result:\n hi " << x.hi << " <> " << cpp_x.hi << "; lo\n " << x.lo << " <> " << cpp_x.lo);
      INFO("C++ ddcsshr scale difference y: " << scaleDiffY << "; expected vs result:\n hi " << y.hi << " <> " << cpp_y.hi << "; lo\n " << y.lo << " <> " << cpp_y.lo);
      bool test_passed = (scaleDiffX >= 20) || (scaleDiffX == 0);
      REQUIRE(test_passed);
      test_passed = (scaleDiffY >= 20) || (scaleDiffY == 0);
      REQUIRE(test_passed);
      
      
      // Call the Fortran version.
      double a_arr[2] = { rec.a_hi, rec.a_lo };
      double x_arr[2];
      double y_arr[2];
      ddcsshr_fortran(a_arr, x_arr, y_arr);
      ddouble fortranResultX(x_arr[0], x_arr[1]);
      ddouble fortranResultY(y_arr[0], y_arr[1]);
      int scaleDiffFortranX = calculate_scale_difference(fortranResultX, x);
      int scaleDiffFortranY = calculate_scale_difference(fortranResultY, y);
      INFO("Fortran ddcsshr scale difference x: " << scaleDiffFortranX << "; expected vs result:\n hi " << x.hi << " <> " << fortranResultX.hi << "; lo\n " << x.lo << " <> " << fortranResultX.lo);
      INFO("Fortran ddcsshr scale difference y: " << scaleDiffFortranY << "; expected vs result:\n hi " << y.hi << " <> " << fortranResultY.hi << "; lo\n " << y.lo << " <> " << fortranResultY.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffFortranX >= 20) || (scaleDiffFortranX == 0);
      REQUIRE(test_passed);
      test_passed = (scaleDiffFortranY >= 20) || (scaleDiffFortranY == 0);
      REQUIRE(test_passed);
   }
}

//
// ddcssnr Test Case
//
TEST_CASE("ddcssnr function test", "[ddfun]") {
   // Binary file format for ddcsub.bin:
   //   2 mp inputs (4 doubles) followed by 2 mp outputs (4 doubles).
   struct ddcssnrRecord {
      double a_hi;
      double a_lo;
      double x_hi;
      double x_lo;
      double y_hi;
      double y_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddcssnr.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddcssnrRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddouble a(rec.a_hi, rec.a_lo);
      ddouble x(rec.x_hi, rec.x_lo);
      ddouble y(rec.y_hi, rec.y_lo);
      ddouble cpp_x, cpp_y;
      
      // Compute the result using the C++ function.
      ddcssnr(a, cpp_x, cpp_y);
      
      int scaleDiffX = calculate_scale_difference(cpp_x, x);
      int scaleDiffY = calculate_scale_difference(cpp_y, y);
      INFO("C++ ddcssnr scale difference x: " << scaleDiffX << "; expected vs result:\n hi " << x.hi << " <> " << cpp_x.hi << "; lo\n " << x.lo << " <> " << cpp_x.lo);
      INFO("C++ ddcssnr scale difference y: " << scaleDiffY << "; expected vs result:\n hi " << y.hi << " <> " << cpp_y.hi << "; lo\n " << y.lo << " <> " << cpp_y.lo);
      bool test_passed = (scaleDiffX >= 20) || (scaleDiffX == 0);
      REQUIRE(test_passed);
      test_passed = (scaleDiffY >= 20) || (scaleDiffY == 0);
      REQUIRE(test_passed);
      
      // Call the Fortran version.
      double a_arr[2] = { rec.a_hi, rec.a_lo };
      double x_arr[2];
      double y_arr[2];
      ddcssnr_fortran(a_arr, x_arr, y_arr);
      ddouble fortranResultX(x_arr[0], x_arr[1]);
      ddouble fortranResultY(y_arr[0], y_arr[1]);
      int scaleDiffFortranX = calculate_scale_difference(fortranResultX, x);
      int scaleDiffFortranY = calculate_scale_difference(fortranResultY, y);
      INFO("Fortran ddcssnr scale difference x: " << scaleDiffFortranX << "; expected vs result:\n hi " << x.hi << " <> " << fortranResultX.hi << "; lo\n " << x.lo << " <> " << fortranResultX.lo);
      INFO("Fortran ddcssnr scale difference y: " << scaleDiffFortranY << "; expected vs result:\n hi " << y.hi << " <> " << fortranResultY.hi << "; lo\n " << y.lo << " <> " << fortranResultY.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffFortranX >= 20) || (scaleDiffFortranX == 0);
      REQUIRE(test_passed);
      test_passed = (scaleDiffFortranY >= 20) || (scaleDiffFortranY == 0);
      REQUIRE(test_passed);
   }
}

//
// ddcadd Test Case
//
TEST_CASE("ddcadd function test", "[ddfun]") {
   // Binary file format for ddcadd.bin:
   //   2 mpc inputs (4 doubles each) followed by 1 mpc outputs (4 doubles).
   struct ddcaddRecord {
      double a_real_hi;
      double a_real_lo;
      double a_imag_hi;
      double a_imag_lo;
      double b_real_hi;
      double b_real_lo;
      double b_imag_hi;
      double b_imag_lo;
      double c_real_hi;
      double c_real_lo;
      double c_imag_hi;
      double c_imag_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddcadd.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddcaddRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddcomplex a(ddouble(rec.a_real_hi, rec.a_real_lo), ddouble(rec.a_imag_hi, rec.a_imag_lo));
      ddcomplex b(ddouble(rec.b_real_hi, rec.b_real_lo), ddouble(rec.b_imag_hi, rec.b_imag_lo));
      ddcomplex c(ddouble(rec.c_real_hi, rec.c_real_lo), ddouble(rec.c_imag_hi, rec.c_imag_lo));
      ddcomplex cpp_c;
      
      // Compute the result using the C++ function.
      cpp_c = a + b;
      
      int scaleDiffReal = calculate_scale_difference(cpp_c.real, c.real);
      int scaleDiffImag = calculate_scale_difference(cpp_c.imag, c.imag);
      INFO("C++ ddcadd scale difference real: " << scaleDiffReal << "; expected vs result:\n hi " << c.real.hi << " <> " << cpp_c.real.hi << "; lo\n " << c.real.lo << " <> " << cpp_c.real.lo);
      INFO("C++ ddcadd scale difference imag: " << scaleDiffImag << "; expected vs result:\n hi " << c.imag.hi << " <> " << cpp_c.imag.hi << "; lo\n " << c.imag.lo << " <> " << cpp_c.imag.lo);
      bool test_passed = (scaleDiffReal >= 20) || (scaleDiffReal == 0);
      REQUIRE(test_passed);
      test_passed = (scaleDiffImag >= 20) || (scaleDiffImag == 0);
      REQUIRE(test_passed);
      
      // Call the Fortran version.
      double a_arr[4] = { rec.a_real_hi, rec.a_real_lo, rec.a_imag_hi, rec.a_imag_lo };
      double b_arr[4] = { rec.b_real_hi, rec.b_real_lo, rec.b_imag_hi, rec.b_imag_lo };
      double c_arr[4];
      ddcadd_fortran(a_arr, b_arr, c_arr);
      ddcomplex fortranResult(ddouble(c_arr[0], c_arr[1]), ddouble(c_arr[2], c_arr[3]));
      scaleDiffReal = calculate_scale_difference(fortranResult.real, c.real);
      scaleDiffImag = calculate_scale_difference(fortranResult.imag, c.imag);
      INFO("Fortran ddcadd scale difference real: " << scaleDiffReal << "; expected vs result:\n hi " << c.real.hi << " <> " << fortranResult.real.hi << "; lo\n " << c.real.lo << " <> " << fortranResult.real.lo);
      INFO("Fortan ddcadd scale difference imag: " << scaleDiffImag << "; expected vs result:\n hi " << c.imag.hi << " <> " << fortranResult.imag.hi << "; lo\n " << c.imag.lo << " <> " << fortranResult.imag.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffReal >= 20) || (scaleDiffReal == 0);
      REQUIRE(test_passed);
      test_passed = (scaleDiffImag >= 20) || (scaleDiffImag == 0);
      REQUIRE(test_passed);

   }

}

//
// ddcsub Test Case
//
TEST_CASE("ddcsub function test", "[ddfun]") {
   // Binary file format for ddcsub.bin:
   //   2 mpc inputs (4 doubles each) followed by 1 mpc outputs (4 doubles).
   struct ddcsubRecord {
      double a_real_hi;
      double a_real_lo;
      double a_imag_hi;
      double a_imag_lo;
      double b_real_hi;
      double b_real_lo;
      double b_imag_hi;
      double b_imag_lo;
      double c_real_hi;
      double c_real_lo;
      double c_imag_hi;
      double c_imag_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddcsub.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddcsubRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddcomplex a(ddouble(rec.a_real_hi, rec.a_real_lo), ddouble(rec.a_imag_hi, rec.a_imag_lo));
      ddcomplex b(ddouble(rec.b_real_hi, rec.b_real_lo), ddouble(rec.b_imag_hi, rec.b_imag_lo));
      ddcomplex c(ddouble(rec.c_real_hi, rec.c_real_lo), ddouble(rec.c_imag_hi, rec.c_imag_lo));
      ddcomplex cpp_c;
      
      // Compute the result using the C++ function.
      cpp_c = a - b;
      
      int scaleDiffReal = calculate_scale_difference(cpp_c.real, c.real);
      int scaleDiffImag = calculate_scale_difference(cpp_c.imag, c.imag);
      INFO("C++ ddcsub scale difference real: " << scaleDiffReal << "; expected vs result:\n hi " << c.real.hi << " <> " << cpp_c.real.hi << "; lo\n " << c.real.lo << " <> " << cpp_c.real.lo);
      INFO("C++ ddcsub scale difference imag: " << scaleDiffImag << "; expected vs result:\n hi " << c.imag.hi << " <> " << cpp_c.imag.hi << "; lo\n " << c.imag.lo << " <> " << cpp_c.imag.lo);
      bool test_passed = (scaleDiffReal >= 20) || (scaleDiffReal == 0);
      REQUIRE(test_passed);
      test_passed = (scaleDiffImag >= 20) || (scaleDiffImag == 0);
      REQUIRE(test_passed);
      
      // Call the Fortran version.
      double a_arr[4] = { rec.a_real_hi, rec.a_real_lo, rec.a_imag_hi, rec.a_imag_lo };
      double b_arr[4] = { rec.b_real_hi, rec.b_real_lo, rec.b_imag_hi, rec.b_imag_lo };
      double c_arr[4];
      ddcsub_fortran(a_arr, b_arr, c_arr);
      ddcomplex fortranResult(ddouble(c_arr[0], c_arr[1]), ddouble(c_arr[2], c_arr[3]));
      scaleDiffReal = calculate_scale_difference(fortranResult.real, c.real);
      scaleDiffImag = calculate_scale_difference(fortranResult.imag, c.imag);
      INFO("Fortran ddcsub scale difference real: " << scaleDiffReal << "; expected vs result:\n hi " << c.real.hi << " <> " << fortranResult.real.hi << "; lo\n " << c.real.lo << " <> " << fortranResult.real.lo);
      INFO("Fortan ddcsub scale difference imag: " << scaleDiffImag << "; expected vs result:\n hi " << c.imag.hi << " <> " << fortranResult.imag.hi << "; lo\n " << c.imag.lo << " <> " << fortranResult.imag.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffReal >= 20) || (scaleDiffReal == 0);
      REQUIRE(test_passed);
      test_passed = (scaleDiffImag >= 20) || (scaleDiffImag == 0);
      REQUIRE(test_passed);
   }
}

//
// ddcmul Test Case
//
TEST_CASE("ddcmul function test", "[ddfun]") {
   // Binary file format for ddcmul.bin:
   //   2 mpc inputs (4 doubles each) followed by 1 mpc outputs (4 doubles).
   struct ddcmulRecord {
      double a_real_hi;
      double a_real_lo;
      double a_imag_hi;
      double a_imag_lo;
      double b_real_hi;
      double b_real_lo;
      double b_imag_hi;
      double b_imag_lo;
      double c_real_hi;
      double c_real_lo;
      double c_imag_hi;
      double c_imag_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddcmul.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddcmulRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddcomplex a(ddouble(rec.a_real_hi, rec.a_real_lo), ddouble(rec.a_imag_hi, rec.a_imag_lo));
      ddcomplex b(ddouble(rec.b_real_hi, rec.b_real_lo), ddouble(rec.b_imag_hi, rec.b_imag_lo));
      ddcomplex c(ddouble(rec.c_real_hi, rec.c_real_lo), ddouble(rec.c_imag_hi, rec.c_imag_lo));
      ddcomplex cpp_c;
      
      // Compute the result using the C++ function.
      cpp_c = a * b;
      
      int scaleDiffReal = calculate_scale_difference(cpp_c.real, c.real);
      int scaleDiffImag = calculate_scale_difference(cpp_c.imag, c.imag);
      INFO("C++ ddcmul scale difference real: " << scaleDiffReal << "; expected vs result:\n hi " << c.real.hi << " <> " << cpp_c.real.hi << "; lo\n " << c.real.lo << " <> " << cpp_c.real.lo);
      INFO("C++ ddcmul scale difference imag: " << scaleDiffImag << "; expected vs result:\n hi " << c.imag.hi << " <> " << cpp_c.imag.hi << "; lo\n " << c.imag.lo << " <> " << cpp_c.imag.lo);
      bool test_passed = (scaleDiffReal >= 20) || (scaleDiffReal == 0);
      REQUIRE(test_passed);
      test_passed = (scaleDiffImag >= 20) || (scaleDiffImag == 0);
      REQUIRE(test_passed);
      
      // Call the Fortran version.
      double a_arr[4] = { rec.a_real_hi, rec.a_real_lo, rec.a_imag_hi, rec.a_imag_lo };
      double b_arr[4] = { rec.b_real_hi, rec.b_real_lo, rec.b_imag_hi, rec.b_imag_lo };
      double c_arr[4];
      ddcmul_fortran(a_arr, b_arr, c_arr);
      ddcomplex fortranResult(ddouble(c_arr[0], c_arr[1]), ddouble(c_arr[2], c_arr[3]));
      scaleDiffReal = calculate_scale_difference(fortranResult.real, c.real);
      scaleDiffImag = calculate_scale_difference(fortranResult.imag, c.imag);
      INFO("Fortran ddcmul scale difference real: " << scaleDiffReal << "; expected vs result:\n hi " << c.real.hi << " <> " << fortranResult.real.hi << "; lo\n " << c.real.lo << " <> " << fortranResult.real.lo);
      INFO("Fortan ddcmul scale difference imag: " << scaleDiffImag << "; expected vs result:\n hi " << c.imag.hi << " <> " << fortranResult.imag.hi << "; lo\n " << c.imag.lo << " <> " << fortranResult.imag.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffReal >= 20) || (scaleDiffReal == 0);
      REQUIRE(test_passed);
      test_passed = (scaleDiffImag >= 20) || (scaleDiffImag == 0);
      REQUIRE(test_passed);
   }
}


//
// ddcdiv Test Case
//
TEST_CASE("ddcdiv function test", "[ddfun]") {
   // Binary file format for ddcdiv.bin:
   //   2 mpc inputs (4 doubles each) followed by 1 mpc outputs (4 doubles).
   struct ddcdivRecord {
      double a_real_hi;
      double a_real_lo;
      double a_imag_hi;
      double a_imag_lo;
      double b_real_hi;
      double b_real_lo;
      double b_imag_hi;
      double b_imag_lo;
      double c_real_hi;
      double c_real_lo;
      double c_imag_hi;
      double c_imag_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddcdiv.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddcdivRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddcomplex a(ddouble(rec.a_real_hi, rec.a_real_lo), ddouble(rec.a_imag_hi, rec.a_imag_lo));
      ddcomplex b(ddouble(rec.b_real_hi, rec.b_real_lo), ddouble(rec.b_imag_hi, rec.b_imag_lo));
      ddcomplex c(ddouble(rec.c_real_hi, rec.c_real_lo), ddouble(rec.c_imag_hi, rec.c_imag_lo));
      ddcomplex cpp_c;
      
      // Compute the result using the C++ function.
      cpp_c = a / b;
      
      int scaleDiffReal = calculate_scale_difference(cpp_c.real, c.real);
      int scaleDiffImag = calculate_scale_difference(cpp_c.imag, c.imag);
      INFO("C++ ddcdiv scale difference real: " << scaleDiffReal << "; expected vs result:\n hi " << c.real.hi << " <> " << cpp_c.real.hi << "; lo\n " << c.real.lo << " <> " << cpp_c.real.lo);
      INFO("C++ ddcdiv scale difference imag: " << scaleDiffImag << "; expected vs result:\n hi " << c.imag.hi << " <> " << cpp_c.imag.hi << "; lo\n " << c.imag.lo << " <> " << cpp_c.imag.lo);
      bool test_passed = (scaleDiffReal >= 20) || (scaleDiffReal == 0);
      REQUIRE(test_passed);
      test_passed = (scaleDiffImag >= 20) || (scaleDiffImag == 0);
      REQUIRE(test_passed);
      
      // Call the Fortran version.
      double a_arr[4] = { rec.a_real_hi, rec.a_real_lo, rec.a_imag_hi, rec.a_imag_lo };
      double b_arr[4] = { rec.b_real_hi, rec.b_real_lo, rec.b_imag_hi, rec.b_imag_lo };
      double c_arr[4];
      ddcdiv_fortran(a_arr, b_arr, c_arr);
      ddcomplex fortranResult(ddouble(c_arr[0], c_arr[1]), ddouble(c_arr[2], c_arr[3]));
      scaleDiffReal = calculate_scale_difference(fortranResult.real, c.real);
      scaleDiffImag = calculate_scale_difference(fortranResult.imag, c.imag);
      INFO("Fortran ddcdiv scale difference real: " << scaleDiffReal << "; expected vs result:\n hi " << c.real.hi << " <> " << fortranResult.real.hi << "; lo\n " << c.real.lo << " <> " << fortranResult.real.lo);
      INFO("Fortan ddcdiv scale difference imag: " << scaleDiffImag << "; expected vs result:\n hi " << c.imag.hi << " <> " << fortranResult.imag.hi << "; lo\n " << c.imag.lo << " <> " << fortranResult.imag.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffReal >= 20) || (scaleDiffReal == 0);
      REQUIRE(test_passed);
      test_passed = (scaleDiffImag >= 20) || (scaleDiffImag == 0);
      REQUIRE(test_passed);
   }
}

//
// ddcsqrt Test Case
//
TEST_CASE("ddcsqrt function test", "[ddfun]") {
   // Binary file format for ddcsqrt.bin:
   //   1 mpc input (4 doubles) followed by 1 mpc output (4 doubles).
   struct ddcsqrtRecord {
      double a_real_hi;
      double a_real_lo;
      double a_imag_hi;
      double a_imag_lo;
      double c_real_hi;
      double c_real_lo;
      double c_imag_hi;
      double c_imag_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddcsqrt.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddcsqrtRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddcomplex a(ddouble(rec.a_real_hi, rec.a_real_lo), ddouble(rec.a_imag_hi, rec.a_imag_lo));
      ddcomplex c(ddouble(rec.c_real_hi, rec.c_real_lo), ddouble(rec.c_imag_hi, rec.c_imag_lo));
      ddcomplex cpp_c;
      
      // Compute the result using the C++ function.
      cpp_c = ddcsqrt(a);
      
      int scaleDiffReal = calculate_scale_difference(cpp_c.real, c.real);
      int scaleDiffImag = calculate_scale_difference(cpp_c.imag, c.imag);
      INFO("C++ ddcsqrt scale difference real: " << scaleDiffReal << "; expected vs result:\n hi " << c.real.hi << " <> " << cpp_c.real.hi << "; lo\n " << c.real.lo << " <> " << cpp_c.real.lo);
      INFO("C++ ddcsqrt scale difference imag: " << scaleDiffImag << "; expected vs result:\n hi " << c.imag.hi << " <> " << cpp_c.imag.hi << "; lo\n " << c.imag.lo << " <> " << cpp_c.imag.lo);
      bool test_passed = (scaleDiffReal >= 20) || (scaleDiffReal == 0);
      REQUIRE(test_passed);
      test_passed = (scaleDiffImag >= 20) || (scaleDiffImag == 0);
      REQUIRE(test_passed);

      // Call the Fortran version.
      double a_arr[4] = { rec.a_real_hi, rec.a_real_lo, rec.a_imag_hi, rec.a_imag_lo };
      double c_arr[4];
      ddcsqrt_fortran(a_arr, c_arr);
      ddcomplex fortranResult(ddouble(c_arr[0], c_arr[1]), ddouble(c_arr[2], c_arr[3]));
      scaleDiffReal = calculate_scale_difference(fortranResult.real, c.real);
      scaleDiffImag = calculate_scale_difference(fortranResult.imag, c.imag);
      INFO("Fortran ddcsqrt scale difference real: " << scaleDiffReal << "; expected vs result:\n hi " << c.real.hi << " <> " << fortranResult.real.hi << "; lo\n " << c.real.lo << " <> " << fortranResult.real.lo);
      INFO("Fortran ddcsqrt scale difference imag: " << scaleDiffImag << "; expected vs result:\n hi " << c.imag.hi << " <> " << fortranResult.imag.hi << "; lo\n " << c.imag.lo << " <> " << fortranResult.imag.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffReal >= 20) || (scaleDiffReal == 0);
      REQUIRE(test_passed);
      test_passed = (scaleDiffImag >= 20) || (scaleDiffImag == 0);
      REQUIRE(test_passed);
   }
}

// 
// ddcpwr Test Case
//
TEST_CASE("ddcpwr function test", "[ddfun]") {
   // Binary file format for ddcpwr.bin:
   //   1 mpc input (4 doubles) and 1 int input, followed by 1 mpc output (4 doubles).
   struct ddcpwrRecord {
      double a_real_hi;
      double a_real_lo;
      double a_imag_hi;
      double a_imag_lo;
      int n;
      int ph;
      double c_real_hi;
      double c_real_lo;
      double c_imag_hi;
      double c_imag_lo;
   };

   std::ifstream infile(INPUT_FILES_DIR + "data/ddcpwr.bin", std::ios::binary);
   REQUIRE(infile.good());
   ddcpwrRecord rec;
   while (infile.read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
      ddcomplex a(ddouble(rec.a_real_hi, rec.a_real_lo), ddouble(rec.a_imag_hi, rec.a_imag_lo));
      ddcomplex c(ddouble(rec.c_real_hi, rec.c_real_lo), ddouble(rec.c_imag_hi, rec.c_imag_lo));
      ddcomplex cpp_c;
      
      // Compute the result using the C++ function.
      cpp_c = ddcpwr(a, rec.n);
      
      int scaleDiffReal = calculate_scale_difference(cpp_c.real, c.real);
      int scaleDiffImag = calculate_scale_difference(cpp_c.imag, c.imag);
      INFO("C++ ddcpwr scale difference real: " << scaleDiffReal << "; expected vs result:\n hi " << c.real.hi << " <> " << cpp_c.real.hi << "; lo\n " << c.real.lo << " <> " << cpp_c.real.lo);
      INFO("C++ ddcpwr scale difference imag: " << scaleDiffImag << "; expected vs result:\n hi " << c.imag.hi << " <> " << cpp_c.imag.hi << "; lo\n " << c.imag.lo << " <> " << cpp_c.imag.lo);
      bool test_passed = (scaleDiffReal >= 20) || (scaleDiffReal == 0);
      REQUIRE(test_passed);
      test_passed = (scaleDiffImag >= 20) || (scaleDiffImag == 0);
      REQUIRE(test_passed);

      // Call the Fortran version.
      double a_arr[4] = { rec.a_real_hi, rec.a_real_lo, rec.a_imag_hi, rec.a_imag_lo };
      double c_arr[4];
      ddcpwr_fortran(a_arr, rec.n, c_arr);
      ddcomplex fortranResult(ddouble(c_arr[0], c_arr[1]), ddouble(c_arr[2], c_arr[3]));
      scaleDiffReal = calculate_scale_difference(fortranResult.real, c.real);
      scaleDiffImag = calculate_scale_difference(fortranResult.imag, c.imag);
      INFO("Fortran ddcpwr scale difference real: " << scaleDiffReal << "; expected vs result:\n hi " << c.real.hi << " <> " << fortranResult.real.hi << "; lo\n " << c.real.lo << " <> " << fortranResult.real.lo);
      INFO("Fortran ddcpwr scale difference imag: " << scaleDiffImag << "; expected vs result:\n hi " << c.imag.hi << " <> " << fortranResult.imag.hi << "; lo\n " << c.imag.lo << " <> " << fortranResult.imag.lo);
      // We expect at least 20 digits of precision.
      test_passed = (scaleDiffReal >= 20) || (scaleDiffReal == 0);
      REQUIRE(test_passed);
      test_passed = (scaleDiffImag >= 20) || (scaleDiffImag == 0);
      REQUIRE(test_passed);
   }
}




