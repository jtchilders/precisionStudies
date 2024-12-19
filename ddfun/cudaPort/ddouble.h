#ifndef DDOUBLE_H
#define DDOUBLE_H

#include <cmath>
#include <cstdio>
#include <cassert>

#ifdef __CUDACC__
#define DD_INLINE __host__ __device__ inline
#else
#define DD_INLINE inline
#endif

class ddouble {
public:
    double hi;
    double lo;

    DD_INLINE ddouble() : hi(0.0), lo(0.0) {}
    DD_INLINE ddouble(double h, double l=0.0) : hi(h), lo(l) {}
    DD_INLINE ddouble(double x) : hi(x), lo(0.0) {}

    // Construct from two doubles representing a "double-double"
    static DD_INLINE ddouble from_two_doubles(double h, double l) {
        ddouble r;
        r.hi = h;
        r.lo = l;
        return r;
    }

    // Basic arithmetic operators
    DD_INLINE ddouble operator+(const ddouble &b) const;
    DD_INLINE ddouble operator-(const ddouble &b) const;
    DD_INLINE ddouble operator*(const ddouble &b) const;
    DD_INLINE ddouble operator/(const ddouble &b) const;

    // Comparison
    DD_INLINE bool operator==(const ddouble &b) const {
        return (hi == b.hi && lo == b.lo);
    }
    DD_INLINE bool operator!=(const ddouble &b) const {
        return !(*this == b);
    }
    DD_INLINE bool operator<(const ddouble &b) const {
        if (hi < b.hi) return true;
        if (hi > b.hi) return false;
        return (lo < b.lo);
    }
    DD_INLINE bool operator>(const ddouble &b) const {
        return b < *this;
    }
    DD_INLINE bool operator<=(const ddouble &b) const {
        return !(b < *this);
    }
    DD_INLINE bool operator>=(const ddouble &b) const {
        return !(*this < b);
    }

    // Utilities
    DD_INLINE static ddouble quick_two_sum(double a, double b);
    DD_INLINE static ddouble two_sum(double a, double b);
    DD_INLINE static ddouble two_prod(double a, double b);

    // Conversion
    DD_INLINE double to_double() const { return hi; }
};

#endif
