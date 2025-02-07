#pragma once
#include "Kokkos_Core.hpp"
#include "ddouble.h"
#include <cmath>
#include <iostream>
#include <iomanip>

// all function declarations:
KOKKOS_INLINE_FUNCTION ddouble ddabs(const ddouble& );
KOKKOS_INLINE_FUNCTION ddouble ddacosh(const ddouble& );
KOKKOS_INLINE_FUNCTION ddouble ddasinh(const ddouble& );
KOKKOS_INLINE_FUNCTION ddouble ddacosh(const ddouble& );
KOKKOS_INLINE_FUNCTION ddouble ddatanh(const ddouble& );
KOKKOS_INLINE_FUNCTION ddouble ddexp(const ddouble& );
KOKKOS_INLINE_FUNCTION ddouble ddlog(const ddouble& );
KOKKOS_INLINE_FUNCTION ddouble ddmuld(const ddouble& , const double& );
KOKKOS_INLINE_FUNCTION ddouble ddmuldd(const double&, const double& );
KOKKOS_INLINE_FUNCTION ddouble ddnint(const ddouble& );
KOKKOS_INLINE_FUNCTION ddouble ddnpwr(const ddouble& , const int& );
KOKKOS_INLINE_FUNCTION ddouble ddpower(const ddouble& , const ddouble& );
KOKKOS_INLINE_FUNCTION ddouble ddsqrt(const ddouble& );

/**
 * Print the double 'x' in both decimal and hex form.
 * The hex form shows its exact 64-bit pattern.
 */
KOKKOS_INLINE_FUNCTION void printDoubleBits(const char *label, double x, const int& ii)
{
   // We'll copy the double bits into a 64-bit integer.
   // A union is a common trick, or we can use memcpy.
   union {
      double d;
      uint64_t u;
   } conv;

   conv.d = x;

   // Use C99's PRIx64 for a portable 64-bit hex format.
   // %.16g prints up to 16 significant digits in decimal (just for reference).
   printf("[%d] %s: decimal=%.16g, bits=0x%016" PRIx64 "\n", ii,
               label, x, conv.u);
}

KOKKOS_INLINE_FUNCTION ddouble ddabs(const ddouble& a) {
   ddouble b;
   if (a.hi >= 0.0) {
      b.hi = a.hi;
      b.lo = a.lo;
   } else {
      b.hi = -a.hi;
      b.lo = -a.lo;
   }
   return b;
}

KOKKOS_INLINE_FUNCTION ddouble ddacosh(const ddouble& a) {
   if (a.hi < 1.0) {
      printf("DDACOSH: Argument is < 1.\n");
      return ddouble();
   }

   ddouble f1{1.0, 0.0};
   ddouble t1, t2, b;

   t1 = a * a;
   t2 = t1 - f1;
   t1 = ddsqrt(t2);
   t2 = a + t1;
   b = ddlog(t2);

   return b;
}

KOKKOS_INLINE_FUNCTION ddouble ddasinh(const ddouble& a) {
   ddouble f1, t1, t2, b;
   f1.hi = 1.0; f1.lo = 0.0;
   
   t1 = a *  a;
   t2 = t1 +  f1;
   t1 = ddsqrt(t2);
   t2 = a +  t1;
   b = ddlog(t2);
   
   return b;
}

KOKKOS_INLINE_FUNCTION ddouble ddatanh(const ddouble& a) {
   // Check if a <= -1 or a >= 1; error.
   if (Kokkos::abs(a.hi) >= 1.0) {
      printf("DDATANH: Argument is <= -1 or >= 1.\n");
      return ddouble();
   }

   ddouble f1 = {1.0, 0.0};
   ddouble t1, t2, t3;

   t1 = f1 +  a;
   t2 = f1 -  a;
   t3 = t1 /  t2;

   t1 = ddlog(t3);

   return ddmuld(t1, 0.5);
}

KOKKOS_INLINE_FUNCTION ddouble dddivd(const ddouble& dda, const double& db) {
   const double split = 134217729.0;
   double t1, t2, a1, a2, b1, b2, cona, conb, e, t11, t12, t21, t22;

   // Compute a DP approximation to the quotient.
   t1 = dda.hi / db;

   // This splits t1 and db into high-order and low-order words.
   cona = t1 * split;
   conb = db * split;
   a1 = cona - (cona - t1);
   b1 = conb - (conb - db);
   a2 = t1 - a1;
   b2 = db - b1;

   // Multiply t1 * db using Dekker's method.
   t12 = t1 * db;
   t22 = (((a1 * b1 - t12) + a1 * b2) + a2 * b1) + a2 * b2;

   // Compute dda - (t12, t22) using Knuth's trick.
   t11 = dda.hi - t12;
   e = t11 - dda.hi;
   t21 = ((-t12 - e) + (dda.hi - (t11 - e))) + dda.lo - t22;

   // Compute high-order word of (t11, t21) and divide by db.
   t2 = (t11 + t21) / db;

   // The result is t1 + t2, after normalization.
   ddouble ddc;
   ddc.hi = t1 + t2;
   ddc.lo = t2 - (ddc.hi - t1);
   return ddc;
}

KOKKOS_INLINE_FUNCTION ddouble ddexp(const ddouble& a) {
   const int nq = 6;
   int i, l1, nz;
   double eps = Kokkos::pow(10.0, -32.0);
   ddouble al2(0.69314718055994529, 2.3190468138462996e-17);
   ddouble f(1.0, 0.0);
   ddouble s0, s1, s2, s3, result;
   
   if (Kokkos::fabs(a.hi) >= 300.0) {
      if (a.hi > 0.0) {
         // Error: Argument is too large
         // Assume ddabrt handles this 
         printf("DDEXP: Argument is too large\n");
         return ddouble();
      } else {
         return ddouble();
      }
   }

   s0 = a / al2;
   s1 = ddnint(s0);

   double t1 = s1.hi;
   nz = t1 + std::copysign(1e-14, t1);

   s2 = al2 * s1;
   s0 = a - s2;

   if (s0.hi == 0.0) {
      s0 = ddouble(1.0, 0.0);
      l1 = 0;
   } else {
      s1 = dddivd(s0,ldexp(1.0, nq)); // replaced pow(2.0,nq)
      s2 = ddouble(1.0, 0.0);
      s3 = ddouble(1.0, 0.0);
      l1 = 0;

      do {
         l1 = l1 + 1;
         if (l1 == 100) {
               // Error: Iteration limit exceeded
               // Assume ddabrt handles this
               printf("DDEXP: Iteration limit exceeded\n");
               return ddouble();
         }

         double t2 = l1;
         s0 = s2 * s1;
         s2 = dddivd(s0,t2);
         s0 = s3 + s2;
         s3 = s0;

      } while (Kokkos::fabs(s2.hi) > eps * Kokkos::fabs(s3.hi));

      for (i = 0; i < nq; i++) {
         s1 = s0 * s0;
         s0 = s1;
      }
   }

   auto pow_nz = ldexp(1.0, nz); // replaced pow(2.0, nz)
   result = ddmuld(s0, pow_nz);
   
   return result;
}

KOKKOS_INLINE_FUNCTION ddouble ddlog(const ddouble &a) {
   if (a.hi <= 0.0) {
      printf("*** DDLOG: Argument is less than or equal to zero.\n");
      return ddouble();
   }

   ddouble b;
   double t1 = a.hi;
   double t2 = Kokkos::log(t1);
   b.hi = t2;
   b.lo = 0.0;

   for (int k = 1; k <= 3; ++k) {
      ddouble s0 = ddexp(b);
      ddouble s1 = a -  s0;
      ddouble s2 = s1 /  s0;
      ddouble s1_new = b +  s2;
      b.hi = s1_new.hi;
      b.lo = s1_new.lo;
   }

   return b;
}


// mixed type math operators
KOKKOS_INLINE_FUNCTION ddouble ddmuld(const ddouble& dda, const double& db){
   const double split = 134217729.0;
   double a1, a2, b1, b2, cona, conb, c11, c21, c2, e, t1, t2;

   // This splits this->hi and db into high-order and low-order words.
   cona = dda.hi * split;
   conb = db * split;
   a1 = cona - (cona - dda.hi);
   b1 = conb - (conb - db);
   a2 = dda.hi - a1;
   b2 = db - b1;

   // Multiply this->hi * db using Dekker's method.
   c11 = dda.hi * db;
   c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2;

   // Compute this->lo * db (only high-order word is needed).
   c2 = dda.lo * db;

   // Compute (c11, c21) + c2 using Knuth's trick.
   t1 = c11 + c2;
   e = t1 - c11;
   t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

   // The result is t1 + t2, after normalization.
   ddouble ddc;
   ddc.hi = t1 + t2;
   ddc.lo = t2 - (ddc.hi - t1);

   return ddc;
}


ddouble ddmuldd(const double& da, const double& db) {
    const double split = 134217729.0;
    double cona, conb, a1, a2, b1, b2, s1, s2;
    ddouble ddc;

    // On systems with a fused multiply add, such as IBM systems, it is faster to
    // uncomment the next two lines and comment out the following lines.
    // On other systems, do the opposite.

    // s1 = da * db;
    // s2 = da * db - s1;

    // This splits da and db into high-order and low-order words.

    cona = da * split;
    conb = db * split;
    a1 = cona - (cona - da);
    b1 = conb - (conb - db);
    a2 = da - a1;
    b2 = db - b1;

    // Multiply da * db using Dekker's method.

    s1 = da * db;
    s2 = (((a1 * b1 - s1) + a1 * b2) + a2 * b1) + a2 * b2;

    ddc.hi = s1;
    ddc.lo = s2;

    return ddc;
}



KOKKOS_INLINE_FUNCTION ddouble ddnint(const ddouble& a) {
   ddouble b = {0.0, 0.0};
   ddouble s0;
   // create some constants
   const double T105 = Kokkos::pow(2.0, 105);
   const double T52 = Kokkos::pow(2.0, 52);
   const ddouble CON(T105, T52);

   // Check if `a` is zero
   if (a.hi == 0.0) {
      return b;
   }

   // Check if `a` is too large
   if (a.hi >= T105) {
      printf("*** DDNINT: Argument is too large.\n");
      return b; // Unreachable if ddabrt terminates, added for safety
   }

   // Perform the rounding logic
   if (a.hi > 0.0) {
      s0 = a + CON;
      b = s0 - CON;
   } else {
      s0 = a - CON;
      b = s0 + CON;
   }

   return b;
}


KOKKOS_INLINE_FUNCTION ddouble ddnpwr(const ddouble& a, int n) {
   const double cl2 = 1.4426950408889633;
   ddouble s0, s1, s2 = {1.0, 0.0};
   double t1;
   int nn, mn, kn, kk;

   if (a.hi == 0.0) {
      if (n >= 0) {
         return {0.0, 0.0};
      } else {
         printf("*** DDNPWR: Argument is zero and N is negative or zero.\n");
         return {0.0, 0.0};
      }
   }

   nn = abs(n);
   if (nn == 0) {
      return {1.0, 0.0};
   } else if (nn == 1) {
      if(n > 0)
         return a;
      else
         // calculate reciprocal
         return ddouble(1.0,0.0) / a;
   } else if (nn == 2) {
      if(n > 0)
         return a *  a;
      else
         // calculate reciprocal
         return ddouble(1.0,0.0) / (a * a);
   }

   // Determine the least integer MN such that 2 ^ MN > NN.
   t1 = static_cast<double>(nn);
   mn = static_cast<int>(cl2 * log(t1) + 1.0 + 1.0e-14);

   s0 = a;
   kn = nn;

   // Compute B ^ N using the binary rule for exponentiation.
   for (int j = 1; j <= mn; ++j) {
      kk = kn / 2;
      if (kn != 2 * kk) {
         s1 = s2 *  s0;
         s2 = s1;
      }
      kn = kk;
      if (j < mn) {
         s1 = s0 *  s0;
         s0 = s1;
      }
   }

   // Compute reciprocal if N is negative.
   if (n < 0) {
      s1 = {1.0, 0.0};
      s0 = s1 /  s2;
      s2 = s0;
   }

   return s2;
}


// must address the variable size of the `ad` array before using on GPU
// KOKKOS_INLINE_FUNCTION ddouble ddpolyr(const int n, const Kokkos::View<ddouble*>& a, const ddouble& x0) {
//     const double eps = 1.0e-29;
//     ddouble ad[n+1];
//     ddouble x = x0;
//     double dt1;
   
//     for (int i = 0; i < n; ++i) {
//         dt1 = static_cast<double>(i + 1);
//         ad[i] = a[i+1] *  dt1;
//     }
//     ad[n].hi = 0.0;
//     ad[n].lo = 0.0;
   
//     for (int it = 0; it < 20; ++it) {
//         ddouble t1 = {0.0, 0.0};
//         ddouble t2 = {0.0, 0.0};
//         ddouble t3 = {1.0, 0.0};
//         ddouble t4, t5;
      
//         for (int i = 0; i <= n; ++i) {
//             t4 = a[i] *  t3;
//             t5 = t1 +  t4;
//             t1 = t5;
//             t4 = ad[i] *  t3;
//             t5 = t2 +  t4;
//             t2 = t5;
//             t4 = t3 *  x;
//             t3 = t4;
//         }
      
//         t3 = t1 /  t2;
//         t4 = x -  t3;
//         x = t4;
//         if (Kokkos::abs(t3.hi) <= eps) {
//             return x;
//         }
//     }
   
//     printf("DDROOT: failed to converge.\n");
//     return ddouble();
// }

KOKKOS_INLINE_FUNCTION ddouble ddpower(const ddouble& a, const ddouble& b) {
   if (a.hi <= 0.0) {
      printf("DDPOWER: A <= 0\n");
      return ddouble();
   }
   ddouble t1 = ddlog(a);
   ddouble t2 = t1 *  b;
   ddouble c = ddexp(t2);
   return c;
}



KOKKOS_INLINE_FUNCTION ddouble ddsqrt(const ddouble& a) {

   if (a.hi == 0.0) {
      return ddouble();
   }

   double t1 = 1.0 / sqrt(a.hi);
   double t2 = a.hi * t1;
   ddouble s0 = ddmuldd(t2, t2);
   ddouble s1 = a -  s0;
   double t3 = 0.5 * s1.hi * t1;
   s0.hi = t2;
   s0.lo = 0.0;
   s1.hi = t3;
   s1.lo = 0.0;
   ddouble b = s0 +  s1;

   return b;
}


// Function to calculate scale difference
inline int calculate_scale_difference(const ddouble &result, const ddouble &expected) {
   double error_hi_abs = std::abs(result.hi - expected.hi);

   if (error_hi_abs > 0.0) {
      double error_hi_exponent = std::log10(error_hi_abs);
      double expected_hi_exponent = std::log10(std::abs(expected.hi));
      int scale_difference = static_cast<int>(std::abs(error_hi_exponent - expected_hi_exponent));
      return scale_difference;
   }

   double error_lo_abs = std::abs(result.lo - expected.lo);
   if (error_lo_abs > 0.0) {
      double error_lo_exponent = std::log10(error_lo_abs);
      double expected_hi_exponent = std::log10(std::abs(expected.hi));
      int scale_difference = static_cast<int>(std::abs(error_lo_exponent - expected_hi_exponent));
      return scale_difference;
   }

   return 0;
}