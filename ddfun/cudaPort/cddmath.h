#ifndef CDDMATH_H
#define CDDMATH_H

#include "complexdd.h"
#include "ddmath.h"

#ifdef __CUDACC__
#define DD_INLINE __host__ __device__ inline
#else
#define DD_INLINE inline
#endif

// Example: complex log
DD_INLINE complexdd cdd_log(const complexdd &z) {
    // log(z) = log(|z|) + i arg(z)
    ddouble r = dd_sqrt(z.re*z.re + z.im*z.im); // magnitude
    ddouble theta; 
    // arg
    // For a robust arg, consider sign of re and im:
    // simplified:
    theta = ddouble(std::atan2(z.im.to_double(), z.re.to_double()));
    return complexdd(dd_log(r), theta);
}

// Similarly implement cdd_exp, cdd_sqrt, cdd_atanh, etc.

#endif
