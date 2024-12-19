#ifndef DDMATH_H
#define DDMATH_H

#include "ddouble.h"
#include <cmath>

#ifdef __CUDACC__
#define DD_INLINE __host__ __device__ inline
#else
#define DD_INLINE inline
#endif

// Constants used in original code
// For demonstration, a few constants:
DD_INLINE ddouble dd_const_pi() {
    // Two-word representation from code:
    return ddouble::from_two_doubles(3.1415926535897931E+00, 1.2246467991473532E-16);
}

// Example: sqrt for ddouble
DD_INLINE ddouble dd_sqrt(const ddouble &a) {
    if (a.hi <= 0.0) {
        // error handling or return something
        // For now, let's assume a.hi>0
        return ddouble(0.0);
    }
    double x = 1.0 / std::sqrt(a.hi);
    double t = a.hi*x;
    // refinement similar to original code
    double s = 0.5 * ( (a.hi - t*t)*x );
    double hi = t + s;
    double lo = 0.0; // simplification; implement full correction from code if needed
    return ddouble::from_two_doubles(hi, lo);
}

// Example: exp for ddouble
// Implementing Newton iteration or direct series is lengthy; here is a placeholder
DD_INLINE ddouble dd_exp(const ddouble &a) {
    // Basic approach: dd_exp(a) ~ exp(a.hi) + ...
    // For simplicity, just do a double-based exp as approximation:
    // A more accurate approach would replicate the original code's steps.
    double eh = std::exp(a.hi);
    // Minimal refinement:
    // TODO: Add corrections from the original Fortran code.
    return ddouble(eh);
}

// Example: log for ddouble
DD_INLINE ddouble dd_log(const ddouble &a) {
    // Similarly, just double-based approximation:
    double lh = std::log(a.hi);
    return ddouble(lh);
}

// Add more functions (atanh, asinh, acosh, etc.) using the dd operations.

#endif
