#include "ddouble.h"
#include "ddmath.h"
#include <cmath>

namespace ddfun{
   

ddouble::ddouble() : hi(0.0), lo(0.0)  {};
ddouble::ddouble(double h) : hi(h), lo(0.0) {};
ddouble::ddouble(double h, double l) : hi(h), lo(l) {};

// assignment operator
ddouble& ddouble::operator=(const ddouble& other) {
   hi = other.hi;
   lo = other.lo;
   return *this;
}


// negation operator
ddouble ddouble::operator-() const{
   return ddneg(*this);
}


// math operators
ddouble ddouble::operator+(const ddouble& ddb) const {
   return ddadd(*this, ddb);
}

ddouble ddouble::operator-(const ddouble& ddb) const {
   return ddsub(*this, ddb);
}

ddouble ddouble::operator*(const ddouble& ddb) const {
   return ddmul(*this, ddb);
}

ddouble ddouble::operator/(const ddouble& ddb) const {
   return dddiv(*this, ddb);
}


// mixed type math operators
ddouble ddouble::operator*(const double& db) const{
   return ddmuld(*this, db);
}

ddouble ddouble::operator/(const double& db) const{
   return dddivd(*this, db);
}

// exp and log
ddouble ddouble::exp() const {
   return ddexp(*this);
}

ddouble ddouble::log() const {
   return ddlog(*this);
}

// other math operators
ddouble ddouble::nint() const {
   return ddnint(*this);
}

ddouble ddouble::abs() const {
   return ddabs(*this);
}

ddouble ddouble::sqrt() const {
   return ddsqrt(*this);
}

ddouble ddouble::npwr(int n) const{
   return ddnpwr(*this, n);
}

ddouble ddouble::power(const ddouble& ddb) const{
   return ddpower(*this, ddb);
}

// hyperbolic functions
ddouble ddouble::acosh() const {
   return ddacosh(*this);
}
ddouble ddouble::asinh() const{
   return ddasinh(*this);
}
ddouble ddouble::atanh() const{
   return ddatanh(*this);
}

// print operator
std::ostream& operator<<(std::ostream& os, const ddouble& dd) {
   os << "[ " <<dd.hi << ", " << dd.lo << " ]";
   return os;
}

}