#ifndef COMPLEXDD_H
#define COMPLEXDD_H

#include "ddouble.h"

#ifdef __CUDACC__
#define DD_INLINE __host__ __device__ inline
#else
#define DD_INLINE inline
#endif

class complexdd {
public:
    ddouble re;
    ddouble im;

    DD_INLINE complexdd() : re(0.0), im(0.0) {}
    DD_INLINE complexdd(const ddouble &r, const ddouble &i) : re(r), im(i) {}
    DD_INLINE complexdd(double r, double i=0.0) : re(r), im(i) {}

    DD_INLINE complexdd operator+(const complexdd &b) const {
        return complexdd(re + b.re, im + b.im);
    }

    DD_INLINE complexdd operator-(const complexdd &b) const {
        return complexdd(re - b.re, im - b.im);
    }

    DD_INLINE complexdd operator*(const complexdd &b) const {
        // (x+yi)*(u+vi) = (xu - yv) + (xv+yu)i
        ddouble xu = re * b.re;
        ddouble yv = im * b.im;
        ddouble xv = re * b.im;
        ddouble yu = im * b.re;
        return complexdd(xu - yv, xv + yu);
    }

    DD_INLINE complexdd operator/(const complexdd &b) const {
        // (a+bi)/(c+di) = [ (a+bi)(c-di) ] / (c^2+d^2)
        ddouble denom = b.re*b.re + b.im*b.im;
        complexdd numerator(re*b.re + im*b.im, im*b.re - re*b.im);
        return complexdd(numerator.re/denom, numerator.im/denom);
    }
};

#endif
