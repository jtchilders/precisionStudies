#ifndef DDOUBLE_H
#define DDOUBLE_H

#include <cmath>
#include <cstdio>
#include <cassert>

#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>

#ifdef __CUDACC__
#define DD_INLINE __host__ __device__ inline
#else
#define DD_INLINE inline
#endif


class ddouble {
public:
    double hi;
    double lo;

    // Constants for splitting
    static constexpr double SPLIT = 134217729.0; // 2^27+1

    DD_INLINE ddouble() : hi(0.0), lo(0.0) {}
    DD_INLINE ddouble(double h, double l) : hi(h), lo(l) {}
    DD_INLINE ddouble(double x) : hi(x), lo(0.0) {}

    // Construct from two doubles representing a "double-double"
    DD_INLINE static ddouble from_two_doubles(double h, double l) {
        ddouble r;
        r.hi = h;
        r.lo = l;
        return r;
    }

    DD_INLINE static ddouble quick_two_sum(double a, double b) {
        double s = a + b;
        double e = b - (s - a);
        return ddouble::from_two_doubles(s, e);
    }

    DD_INLINE static ddouble two_sum(double a, double b) {
        printf("two_sum: a = %.17g, b = %.17g\n", a, b);
        double s = a + b;
        printf("two_sum: s = %.17g\n", s);
        double bb = s - a;
        printf("two_sum: bb = %.17g\n", bb);
        double e = (a - (s - bb)) + (b - bb);
        printf("two_sum: e = %.17g\n", e);
        return ddouble::from_two_doubles(s, e);
    }

    DD_INLINE static ddouble two_prod(double a, double b) {
        double p = a * b;
        double a_h = (double)((long long)(a * SPLIT));
        a_h = a_h - (a_h - a);
        double a_l = a - a_h;
        double b_h = (double)((long long)(b * SPLIT));
        b_h = b_h - (b_h - b);
        double b_l = b - b_h;
        double e = ((a_h * b_h - p) + a_h * b_l + a_l * b_h) + a_l * b_l;
        return from_two_doubles(p, e);
    }

    // ddouble addition
    DD_INLINE ddouble operator+(const ddouble &b) const {
        printf("ddouble addition: hi = %.17g, lo = %.17g, b.hi = %.17g, b.lo = %.17g\n", hi, lo, b.hi, b.lo);
        // Using Knuth's TwoSum
        ddouble t = ddouble::two_sum(hi, b.hi);
        printf("ddouble addition: t.hi = %.17g, t.lo = %.17g\n", t.hi, t.lo);
        double e = lo + b.lo;
        printf("ddouble addition: e = %.17g\n", e);
        ddouble r = ddouble::two_sum(t.hi, t.lo + e);
        printf("ddouble addition: r.hi = %.17g, r.lo = %.17g\n", r.hi, r.lo);
        return r;
    }

    DD_INLINE static ddouble ddadd(const ddouble &dda,const ddouble &ddb) {
        auto t1 = dda(1) + ddb(1);
        auto e = t1 - dda(1);
        auto t2 = ((ddb(1) - e) + (dda(1) - (t1 - e))) + dda(2) + ddb(2);

        //   The result is t1 + t2, after normalization.
        auto combo = t1 + t2;
        return ddouble::from_two_doubles(combo, t2 - (combo - t1) );
    }

    // ddouble subtraction
    DD_INLINE ddouble operator-(const ddouble &b) const {
        printf("ddouble subtraction: hi = %.17g, lo = %.17g, b.hi = %.17g, b.lo = %.17g\n", hi, lo, b.hi, b.lo);
        // 1) main difference of the high parts
        ddouble s1 = two_sum(hi, -b.hi); // call that result (hi1, lo1)
        printf("ddouble subtraction: s1.hi = %.17g, s1.lo = %.17g\n", s1.hi, s1.lo);

        // 2) combine leftover: (a.lo - b.lo)
        double e = lo - b.lo;
        printf("ddouble subtraction: e = %.17g\n", e);

        // 3) s2 = two_sum(s1.lo, e) merges leftover from step1 with e
        ddouble s2 = two_sum(s1.lo, e);
        printf("ddouble subtraction: s2.hi = %.17g, s2.lo = %.17g\n", s2.hi, s2.lo);

        // 4) final = two_sum(s1.hi, s2.hi)
        ddouble final_1 = two_sum(s1.hi, s2.hi);
        printf("ddouble subtraction: final_1.hi = %.17g, final_1.lo = %.17g\n", final_1.hi, final_1.lo);

        // leftover final_2 = s2.lo
        double leftover = final_1.lo + s2.lo;
        printf("ddouble subtraction: leftover = %.17g\n", leftover);

        // 5) final = two_sum(final_1.hi, leftover)
        ddouble final_result = two_sum(final_1.hi, leftover);
        printf("ddouble subtraction: final_result.hi = %.17g, final_result.lo = %.17g\n", final_result.hi, final_result.lo);

        return final_result;
    }

    // ddouble multiplication
    DD_INLINE ddouble operator*(const ddouble &b) const {
        printf("ddouble multiplication: hi = %.17g, lo = %.17g, b.hi = %.17g, b.lo = %.17g\n", hi, lo, b.hi, b.lo);
        // Step 1: main product hi parts
        ddouble p1 = two_prod(hi, b.hi); 
        // => p1.hi ~ hi*b.hi (main part), p1.lo ~ leftover from that product
        printf("ddouble multiplication: p1.hi = %.17g, p1.lo = %.17g\n", p1.hi, p1.lo);

        // Step 2: cross terms
        double cross = (hi * b.lo + lo * b.hi);
        printf("ddouble multiplication: cross = %.17g\n", cross);
        ddouble sum2 = two_sum(p1.lo, cross);
        // => sum2.hi ~ p1.lo + cross, sum2.lo ~ leftover
        printf("ddouble multiplication: sum2.hi = %.17g, sum2.lo = %.17g\n", sum2.hi, sum2.lo);

        // Step 3: partial result
        ddouble partial_hi = two_sum(p1.hi, sum2.hi);
        printf("ddouble multiplication: partial_hi.hi = %.17g, partial_hi.lo = %.17g\n", partial_hi.hi, partial_hi.lo);

        // Step 4: Merge leftover from partial_hi and sum2
        double final_mid = partial_hi.lo + sum2.lo;
        printf("ddouble multiplication: final_mid = %.17g\n", final_mid);

        // Step 5: multiply lo*lo if we want full coverage
        ddouble small = two_prod(lo, b.lo);
        printf("ddouble multiplication: small.hi = %.17g, small.lo = %.17g\n", small.hi, small.lo);
        // Now we add small.hi into final_mid carefully
        ddouble final_mid_dd = two_sum(final_mid, small.hi);
        printf("ddouble multiplication: final_mid_dd.hi = %.17g, final_mid_dd.lo = %.17g\n", final_mid_dd.hi, final_mid_dd.lo);

        // Step 6: Merge final_mid_dd.hi with partial_hi.hi
        ddouble final_hi_dd = two_sum(partial_hi.hi, final_mid_dd.hi);
        // That merges all the 'lo' parts in a raw double. Next, we do a two_sum with final_hi_dd.hi:
        printf("ddouble multiplication: final_hi_dd.hi = %.17g, final_hi_dd.lo = %.17g\n", final_hi_dd.hi, final_hi_dd.lo);

        // Step 7: Summation of leftover pieces
        // We have leftover: final_hi_dd.lo, final_mid_dd.lo, small.lo
        // Let's sum them all in one place:
        double leftover1 = final_hi_dd.lo + final_mid_dd.lo;
        double leftover2 = leftover1 + small.lo;
        printf("ddouble multiplication: leftover1 = %.17g, leftover2 = %.17g\n", leftover1, leftover2);

        // Step 8: Now combine leftover2 with final_hi_dd.hi carefully
        ddouble result = two_sum(final_hi_dd.hi, leftover2);
        printf("ddouble multiplication: result.hi = %.17g, result.lo = %.17g\n", result.hi, result.lo);

        return result;
    }

    // Overload operator* to multiply ddouble by double
    ddouble operator*(const double other) {
        const double split = 134217729.0; // 2^27 + 1, used for Dekker's splitting method

        // Split hi and other into high-order and low-order parts
        double cona = hi * split;
        double conb = other * split;
        double a1 = cona - (cona - hi);
        double b1 = conb - (conb - other);
        double a2 = hi - a1;
        double b2 = other - b1;

        // Multiply hi * other using Dekker's method
        double c11 = hi * other;
        double c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2;

        // Compute lo * other (only high-order part is needed)
        double c2 = lo * other;

        // Combine (c11, c21) + c2 using Knuth's trick
        double t1 = c11 + c2;
        double e = t1 - c11;
        double t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

        // Normalize the result
        double result_hi = t1 + t2;
        double result_lo = t2 - (result_hi - t1);

        return ddouble(result_hi, result_lo);
    }


    // ddouble division
    DD_INLINE ddouble operator/(const ddouble &b) const {
        // Rough idea: use double division for approximation, then refine
        double q1 = hi / b.hi;
        // Compute remainder: (a - q1*b)
        ddouble prod = ddouble(q1) * b;
        ddouble diff = (*this - prod);
        double q2 = diff.hi / b.hi;
        ddouble r = ddouble::from_two_doubles(q1 + q2, 0.0);
        return r;
    }

    // ddouble greater than
    DD_INLINE bool operator>(const ddouble &b) const {
        return hi > b.hi || (hi == b.hi && lo > b.lo);
    }

    // ddouble less than
    DD_INLINE bool operator<(const ddouble &b) const {
        return hi < b.hi || (hi == b.hi && lo < b.lo);
    }

    // ddouble greater than or equal to
    DD_INLINE bool operator>=(const ddouble &b) const {
        return hi >= b.hi || (hi == b.hi && lo >= b.lo);
    }

    // ddouble less than or equal to
    DD_INLINE bool operator<=(const ddouble &b) const {
        return hi <= b.hi || (hi == b.hi && lo <= b.lo);
    }

    DD_INLINE ddouble npwr(int n) const {
        const double cl2 = 1.4426950408889633; // log2(e)
        ddouble b;
        if (hi == 0.0 && lo == 0.0) {
            if (n >= 0) {
                b.hi = 0.0;
                b.lo = 0.0;
            } else {
                printf("ddnpwr: Argument is zero and N is negative.");
            }
            return b;
        }

        int nn = std::abs(n);

        if (nn == 0) {
            b.hi = 1.0;
            b.lo = 0.0;
            return b;
        } else if (nn == 1) {
            b.hi = hi;
            b.lo = lo;
            return b;
        } else if (nn == 2) {
            b = *this * *this;
            return b;
        }

        // Determine the least integer MN such that 2^MN > NN.
        double t1 = static_cast<double>(nn);
        int mn = static_cast<int>(cl2 * std::log(t1) + 1.0 + 1.0e-14);

        ddouble s0 = *this;
        ddouble s2(1.0, 0.0);
        int kn = nn;

        // Compute a^n using the binary method for exponentiation.
        for (int j = 1; j <= mn; ++j) {
            int kk = kn / 2;
            if (kn != 2 * kk) {
                s2 = s2 * s0;
            }
            kn = kk;
            if (j < mn) {
                s0 =  s0 * s0;
            }
        }

        // Compute reciprocal if N is negative.
        if (n < 0) {
            ddouble one(1.0, 0.0);
            s2 = one / s2;
        }

        b = s2;
    }

    DD_INLINE double operator()(int n) const { return 1 ? hi : lo; }

    // Conversion
    DD_INLINE double to_double() const { return hi; }

    std::string str() const;
    std::string strs() const;
};


#endif
