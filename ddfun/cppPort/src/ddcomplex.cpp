#include "ddcomplex.h"

ddcomplex::ddcomplex() : real(0.0), imag(0.0) {};
ddcomplex::ddcomplex(const ddouble& real) : real(real), imag(0.0) {};
ddcomplex::ddcomplex(const double& real) : real(real), imag(0.0) {};
ddcomplex::ddcomplex(const ddouble& real, const ddouble& imag) : real(real), imag(imag) {};
ddcomplex::ddcomplex(const double& real, const double& imag) : real(real), imag(imag) {};

ddcomplex& ddcomplex::operator=(const ddcomplex& other) {
   real = other.real;
   imag = other.imag;
   return *this;
}

ddcomplex& ddcomplex::operator=(const ddouble& other) {
   real = other;
   imag = 0.0;
   return *this;
}

ddcomplex ddcomplex::operator+(const ddcomplex& b) const {
   return ddcomplex(real + b.real, imag + b.imag);
}

ddcomplex ddcomplex::operator-(const ddcomplex& b) const {
   return ddcomplex(real - b.real, imag - b.imag);
}

ddcomplex ddcomplex::operator*(const ddcomplex& b) const {
   
   ddouble s0, s1, s2, s3;

   s0 = real * b.real;
   s1 = imag * b.imag;
   s2 = real * b.imag;
   s3 = imag * b.real;

   return ddcomplex(s0 - s1, s2 + s3);
}

// following the ddcdiv function from the fortran code
ddcomplex ddcomplex::operator/(const ddcomplex& b) const {
   if (b.real.hi == 0.0 && b.imag.hi == 0.0) {
      std::cerr << "DDCOMPLEX: Division by zero." << std::endl;
      return ddcomplex();
   }

   ddouble f(1.0, 0.0), s0,s1,s2,s3,s4;
   s0 = real * b.real;
   s1 = imag * b.imag;
   s2 = s0 + s1;
   s3 = s0 - s1;
   s0  = real + imag;
   s1  = b.real - b.imag;
   s4  = s0 * s1;
   s1  = s4 - s3;
   s0  = b.real * b.real;
   s3  = b.imag * b.imag;
   s4  = (b.real * b.real) + (b.imag * b.imag);
   s0  = f / s4;

   return ddcomplex(s2 * s0, s1 * s0);
}