// ddouble.h
#pragma once
#include <iostream>
#include "Kokkos_Core.hpp"

struct ddouble{
   double hi;
   double lo;

   // constructors
   KOKKOS_INLINE_FUNCTION ddouble(): hi(0.0), lo(0.0)  {};
   KOKKOS_INLINE_FUNCTION ddouble(double h) : hi(h), lo(0.0) {};
   KOKKOS_INLINE_FUNCTION ddouble(double h, double l) : hi(h), lo(l) {};

   // negation operator
   KOKKOS_INLINE_FUNCTION ddouble operator-() const{
      ddouble b;
      b.hi = -this->hi;
      b.lo = -this->lo;
      return b;
   }


   // assignment operator
   KOKKOS_INLINE_FUNCTION ddouble& operator=(const ddouble& other){
      hi = other.hi;
      lo = other.lo;
      return *this;
   }
   
   // math operators
   KOKKOS_INLINE_FUNCTION ddouble operator+(const ddouble& ddb) const {
      ddouble ddc;
      double e, t1, t2;

      // Compute dda + ddb using Knuth's trick.
      t1 = this->hi + ddb.hi;
      e = t1 - this->hi;
      t2 = ((ddb.hi - e) + (this->hi - (t1 - e))) + this->lo + ddb.lo;

      // The result is t1 + t2, after normalization.
      ddc.hi = t1 + t2;
      ddc.lo = t2 - (ddc.hi - t1);
      return ddc;
   }

   KOKKOS_INLINE_FUNCTION ddouble operator-(const ddouble& ddb) const {
      ddouble ddc;
      double e, t1, t2;

      // Compute dda + ddb using Knuth's trick.
      t1 = this->hi - ddb.hi;
      e = t1 - this->hi;
      t2 = ((-ddb.hi - e) + (this->hi - (t1 - e))) + this->lo - ddb.lo;

      // The result is t1 + t2, after normalization.
      ddc.hi = t1 + t2;
      ddc.lo = t2 - (ddc.hi - t1);
      return ddc;
   }  

   KOKKOS_INLINE_FUNCTION ddouble operator*(const ddouble& ddb) const {
         const double split = 134217729.0;
      double a1, a2, b1, b2, cona, conb, c11, c21, c2, e, t1, t2;
      ddouble ddc;

      // This splits this->hi and ddb.hi into high-order and low-order words.
      cona = this->hi * split;
      conb = ddb.hi * split;
      a1 = cona - (cona - this->hi);
      b1 = conb - (conb - ddb.hi);
      a2 = this->hi - a1;
      b2 = ddb.hi - b1;

      // Multiply this->hi * ddb.hi using Dekker's method.
      c11 = this->hi * ddb.hi;
      c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2;

      // Compute this->hi * ddb.lo + this->lo * ddb.hi (only high-order word is needed).
      c2 = this->hi * ddb.lo + this->lo * ddb.hi;

      // Compute (c11, c21) + c2 using Knuth's trick, also adding low-order product.
      t1 = c11 + c2;
      e = t1 - c11;
      t2 = ((c2 - e) + (c11 - (t1 - e))) + c21 + this->lo * ddb.lo;

      // The result is t1 + t2, after normalization.
      ddc.hi = t1 + t2;
      ddc.lo = t2 - (ddc.hi - t1);

      return ddc;
   }

   KOKKOS_INLINE_FUNCTION ddouble operator/(const ddouble& ddb) const {
      const double split = 134217729.0;
   
      // Compute a DDR approximation to the quotient.
      double s1 = this->hi / ddb.hi;
      
      // This splits s1 and ddb.hi into high-order and low-order words.
      double cona = s1 * split;
      double conb = ddb.hi * split;
      double a1 = cona - (cona - s1);
      double b1 = conb - (conb - ddb.hi);
      double a2 = s1 - a1;
      double b2 = ddb.hi - b1;

      // Multiply s1 * ddb.hi using Dekker's method.
      double c11 = s1 * ddb.hi;
      double c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2;

      // Compute s1 * ddb.lo (only high-order word is needed).
      double c2 = s1 * ddb.lo;

      // Compute (c11, c21) + c2 using Knuth's trick.
      double t1 = c11 + c2;
      double e = t1 - c11;
      double t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

      // The result is t1 + t2, after normalization.
      double t12 = t1 + t2;
      double t22 = t2 - (t12 - t1);

      // Compute dda - (t12, t22) using Knuth's trick.
      double t11 = this->hi - t12;
      e = t11 - this->hi;
      double t21 = ((-t12 - e) + (this->hi - (t11 - e))) + this->lo - t22;

      // Compute high-order word of (t11, t21) and divide by ddb.hi.
      double s2 = (t11 + t21) / ddb.hi;

      // The result is s1 + s2, after normalization.
      ddouble ddc;
      ddc.hi = s1 + s2;
      ddc.lo = s2 - (ddc.hi - s1);

      return ddc;
   }


};

// print operator
std::ostream& operator<<(std::ostream& os, const ddouble& dd) {
   os << "[ " <<dd.hi << ", " << dd.lo << " ]";
   return os;
}