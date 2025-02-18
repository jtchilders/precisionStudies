// ddouble.h
#pragma once
#include <iostream>

namespace ddfun{
   
struct ddouble{
   double hi;
   double lo;

   // constructors
   ddouble();
   ddouble(double h);
   ddouble(double h, double l);

   // assignment operator
   ddouble& operator=(const ddouble& other);

   // negation operator
   ddouble operator-() const;

   // math operators
   ddouble operator+(const ddouble& ddb) const;
   ddouble operator-(const ddouble& ddb) const;
   ddouble operator*(const ddouble& ddb) const;
   ddouble operator/(const ddouble& ddb) const;

   // mixed type math operators
   ddouble operator*(const double& db) const;
   ddouble operator/(const double& db) const;

   // exp and log
   ddouble exp() const;
   ddouble log() const;

   // other math functions
   ddouble nint() const;
   ddouble abs() const;
   ddouble sqrt() const;
   ddouble npwr(int n) const;
   ddouble power(const ddouble& ddb) const;

   // hyperbolic functions
   ddouble acosh() const;
   ddouble asinh() const;
   ddouble atanh() const;

   // print operator
   friend std::ostream& operator<<(std::ostream& os, const ddouble& dd);
};


}
