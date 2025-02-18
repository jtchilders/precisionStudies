#pragma once
#include "ddouble.h"

namespace ddfun{
   

// ddcomplex data type which uses ddouble precision for real and imag
struct ddcomplex {
   ddouble real;
   ddouble imag;

   ddcomplex();
   ddcomplex(const ddouble& real);
   ddcomplex(const double& real);
   ddcomplex(const ddouble& real, const ddouble& imag);
   ddcomplex(const double& real, const double& imag);

   ddcomplex& operator=(const ddcomplex& other);
   ddcomplex& operator=(const ddouble& other);

   ddcomplex operator+(const ddcomplex& b) const;
   ddcomplex operator-(const ddcomplex& b) const;
   ddcomplex operator*(const ddcomplex& b) const;
   ddcomplex operator/(const ddcomplex& b) const;

};

}