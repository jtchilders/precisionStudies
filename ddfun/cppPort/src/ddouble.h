// ddouble.h
#pragma once
#include <iostream>

struct ddouble{
   double hi;
   double lo;

   // constructors
   ddouble() : hi(0.0), lo(0.0) {};
   ddouble(double h) : hi(h), lo(0.0) {};
   ddouble(double h, double l) : hi(h), lo(l) {};

   // assignment operator
   ddouble& operator=(const ddouble& other) {
      hi = other.hi;
      lo = other.lo;
      return *this;
   }

   // print operator
   friend std::ostream& operator<<(std::ostream& os, const ddouble& dd) {
      os << "[ " <<dd.hi << ", " << dd.lo << " ]";
      return os;
   }
};