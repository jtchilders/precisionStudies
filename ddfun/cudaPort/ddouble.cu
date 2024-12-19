#include "ddouble.h"

#ifdef __CUDACC__
#define DD_INLINE __host__ __device__ inline
#else
#define DD_INLINE inline
#endif

// Constants for splitting
static const double SPLIT = 134217729.0; // 2^27+1

DD_INLINE ddouble ddouble::quick_two_sum(double a, double b) {
    double s = a + b;
    double e = b - (s - a);
    return ddouble::from_two_doubles(s, e);
}

DD_INLINE ddouble ddouble::two_sum(double a, double b) {
    double s = a + b;
    double bb = s - a;
    double e = (a - (s - bb)) + (b - bb);
    return ddouble::from_two_doubles(s, e);
}

DD_INLINE ddouble ddouble::two_prod(double a, double b) {
    double p = a * b;
    double a_h = (double)((long long)(a * SPLIT));
    a_h = a_h - (a_h - a);
    double a_l = a - a_h;
    double b_h = (double)((long long)(b * SPLIT));
    b_h = b_h - (b_h - b);
    double b_l = b - b_h;
    double e = ((a_h * b_h - p) + a_h * b_l + a_l * b_h) + a_l * b_l;
    return ddouble::from_two_doubles(p, e);
}

// ddouble addition
DD_INLINE ddouble ddouble::operator+(const ddouble &b) const {
    // Using Knuth's TwoSum
    ddouble t = two_sum(hi, b.hi);
    double e = lo + b.lo;
    ddouble r = two_sum(t.hi, t.lo + e);
    return r;
}

// ddouble subtraction
DD_INLINE ddouble ddouble::operator-(const ddouble &b) const {
    ddouble t = two_sum(hi, -b.hi);
    double e = lo - b.lo;
    ddouble r = two_sum(t.hi, t.lo + e);
    return r;
}

// ddouble multiplication
DD_INLINE ddouble ddouble::operator*(const ddouble &b) const {
    // (hi,lo)*(b.hi,b.lo)
    ddouble p1 = two_prod(hi, b.hi);
    double p2 = hi * b.lo + lo * b.hi;
    ddouble s = two_sum(p1.lo, p2);
    double rh = p1.hi + s.hi;
    double rl = s.lo;
    ddouble r = ddouble::from_two_doubles(rh, rl);
    return r;
}

// ddouble division
DD_INLINE ddouble ddouble::operator/(const ddouble &b) const {
    // Rough idea: use double division for approximation, then refine
    double q1 = hi / b.hi;
    // Compute remainder: (a - q1*b)
    ddouble prod = ddouble(q1) * b;
    ddouble diff = (*this - prod);
    double q2 = diff.hi / b.hi;
    ddouble r = ddouble::from_two_doubles(q1 + q2, 0.0);
    return r;
}
