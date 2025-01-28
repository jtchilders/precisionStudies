#include "ddmath.h"
#include <cmath>
#include <iostream>
#include <iomanip>

ddouble ddabs(const ddouble& a) {
    ddouble b;
    if (a.hi >= 0.0) {
        b.hi = a.hi;
        b.lo = a.lo;
    } else {
        b.hi = -a.hi;
        b.lo = -a.lo;
    }
    return b;
}

ddouble ddacosh(const ddouble& a) {
  if (a.hi < 1.0) {
    std::cerr << "DDACOSH: Argument is < 1." << std::endl;
    std::abort();
  }

  ddouble f1{1.0, 0.0};
  ddouble t1, t2, b;

  t1 = ddmul(a, a);
  t2 = ddsub(t1, f1);
  t1 = ddsqrt(t2);
  t2 = ddadd(a, t1);
  b = ddlog(t2);

  return b;
}

ddouble ddadd(const ddouble& dda, const ddouble& ddb) {
    ddouble ddc;
    double e, t1, t2;

    // Compute dda + ddb using Knuth's trick.
    t1 = dda.hi + ddb.hi;
    e = t1 - dda.hi;
    t2 = ((ddb.hi - e) + (dda.hi - (t1 - e))) + dda.lo + ddb.lo;

    // The result is t1 + t2, after normalization.
    ddc.hi = t1 + t2;
    ddc.lo = t2 - (ddc.hi - t1);
    return ddc;
}

ddouble ddasinh(const ddouble& a) {
    ddouble f1, t1, t2, b;
    f1.hi = 1.0; f1.lo = 0.0;
    
    t1 = ddmul(a, a);
    t2 = ddadd(t1, f1);
    t1 = ddsqrt(t2);
    t2 = ddadd(a, t1);
    b = ddlog(t2);
    
    return b;
}

ddouble ddatanh(const ddouble& a) {
    // Check if a <= -1 or a >= 1; error.
    if (std::abs(a.hi) >= 1.0) {
        std::cerr << "DDATANH: Argument is <= -1 or >= 1." << std::endl;
        std::abort();
    }

    ddouble f1 = {1.0, 0.0};
    ddouble t1, t2, t3;

    t1 = ddadd(f1, a);
    t2 = ddsub(f1, a);
    t3 = dddiv(t1, t2);
    t1 = ddlog(t3);
    return ddmuld(t1, 0.5);
}

ddouble dddiv(const ddouble& dda, const ddouble& ddb) {
    const double split = 134217729.0;
    
    // Compute a DDR approximation to the quotient.
    double s1 = dda.hi / ddb.hi;
    
    // This splits s1 and ddb.hi into high-order and low-order words.
    double cona = s1 * split;
    double conb = ddb.hi * split;
    double a1 = cona - (cona - s1);
    double b1 = conb - (conb - ddb.hi);
    double a2 = s1 - a1;
    double b2 = ddb.hi - b1;

    // Multiply s1 * ddb.hi using Dekker's method.
    double c11 = s1 * ddb.hi;
    double c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2;

    // Compute s1 * ddb.lo (only high-order word is needed).
    double c2 = s1 * ddb.lo;

    // Compute (c11, c21) + c2 using Knuth's trick.
    double t1 = c11 + c2;
    double e = t1 - c11;
    double t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

    // The result is t1 + t2, after normalization.
    double t12 = t1 + t2;
    double t22 = t2 - (t12 - t1);

    // Compute dda - (t12, t22) using Knuth's trick.
    double t11 = dda.hi - t12;
    e = t11 - dda.hi;
    double t21 = ((-t12 - e) + (dda.hi - (t11 - e))) + dda.lo - t22;

    // Compute high-order word of (t11, t21) and divide by ddb.hi.
    double s2 = (t11 + t21) / ddb.hi;

    // The result is s1 + s2, after normalization.
    ddouble ddc;
    ddc.hi = s1 + s2;
    ddc.lo = s2 - (ddc.hi - s1);

    return ddc;
}

ddouble dddivd(const ddouble& dda, const double db) {
    const double split = 134217729.0;
    double t1, t2, a1, a2, b1, b2, cona, conb, e, t11, t12, t21, t22;

    // Compute a DP approximation to the quotient.
    t1 = dda.hi / db;

    // This splits t1 and db into high-order and low-order words.
    cona = t1 * split;
    conb = db * split;
    a1 = cona - (cona - t1);
    b1 = conb - (conb - db);
    a2 = t1 - a1;
    b2 = db - b1;

    // Multiply t1 * db using Dekker's method.
    t12 = t1 * db;
    t22 = (((a1 * b1 - t12) + a1 * b2) + a2 * b1) + a2 * b2;

    // Compute dda - (t12, t22) using Knuth's trick.
    t11 = dda.hi - t12;
    e = t11 - dda.hi;
    t21 = ((-t12 - e) + (dda.hi - (t11 - e))) + dda.lo - t22;

    // Compute high-order word of (t11, t21) and divide by db.
    t2 = (t11 + t21) / db;

    // The result is t1 + t2, after normalization.
    ddouble ddc;
    ddc.hi = t1 + t2;
    ddc.lo = t2 - (ddc.hi - t1);
    return ddc;
}

ddouble ddexp(const ddouble& a) {
    const int nq = 6;
    int i, l1, nz;
    double eps = std::pow(10.0, -32.0);
    ddouble al2 = {0.69314718055994529, 2.3190468138462996e-17};
    ddouble f = {1.0, 0.0};
    ddouble s0, s1, s2, s3, result;
    
    if (fabs(a.hi) >= 300.0) {
        if (a.hi > 0.0) {
            // Error: Argument is too large
            // Assume ddabrt handles this 
            std::cerr << "DDEXP: Argument is too large" << std::endl;
            return ddouble();
        } else {
            return ddouble();
        }
    }

    s0 = dddiv(a, al2);
    s1 = ddnint(s0);
    double t1 = s1.hi;
    nz = t1 + copysign(1e-14, t1);
    s2 = ddmul(al2, s1);
    s0 = ddsub(a, s2);

    if (s0.hi == 0.0) {
        s0 = {1.0, 0.0};
        l1 = 0;
    } else {
        s1 = dddivd(s0, pow(2.0, nq));
        s2 = {1.0, 0.0};
        s3 = {1.0, 0.0};
        l1 = 0;

        do {
            l1 = l1 + 1;
            if (l1 == 100) {
                // Error: Iteration limit exceeded
                // Assume ddabrt handles this
                std::cerr << "DDEXP: Iteration limit exceeded" << std::endl;
                return ddouble(0.0);
            }

            double t2 = l1;
            s0 = ddmul(s2, s1);
            s2 = dddivd(s0, t2);
            s0 = ddadd(s3, s2);
            s3 = s0;

        } while (fabs(s2.hi) > eps * fabs(s3.hi));

        for (i = 0; i < nq; i++) {
            s1 = ddmul(s0, s0);
            s0 = s1;
        }
    }

    result = ddmuld(s0, pow(2.0, nz));
    
    return result;
}

ddouble ddlog(const ddouble &a) {
    if (a.hi <= 0.0) {
        std::cerr << "*** DDLOG: Argument is less than or equal to zero." << std::endl;
        return ddouble();
    }

    ddouble b;
    double t1 = a.hi;
    double t2 = std::log(t1);
    b.hi = t2;
    b.lo = 0.0;

    for (int k = 1; k <= 3; ++k) {
        ddouble s0 = ddexp(b);
        ddouble s1 = ddsub(a, s0);
        ddouble s2 = dddiv(s1, s0);
        ddouble s1_new = ddadd(b, s2);
        b.hi = s1_new.hi;
        b.lo = s1_new.lo;
    }

    return b;
}


ddouble ddmul(const ddouble& dda, const ddouble& ddb) {
    const double split = 134217729.0;
    double a1, a2, b1, b2, cona, conb, c11, c21, c2, e, t1, t2;
    ddouble ddc;

    // This splits dda.hi and ddb.hi into high-order and low-order words.
    cona = dda.hi * split;
    conb = ddb.hi * split;
    a1 = cona - (cona - dda.hi);
    b1 = conb - (conb - ddb.hi);
    a2 = dda.hi - a1;
    b2 = ddb.hi - b1;

    // Multiply dda.hi * ddb.hi using Dekker's method.
    c11 = dda.hi * ddb.hi;
    c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2;

    // Compute dda.hi * ddb.lo + dda.lo * ddb.hi (only high-order word is needed).
    c2 = dda.hi * ddb.lo + dda.lo * ddb.hi;

    // Compute (c11, c21) + c2 using Knuth's trick, also adding low-order product.
    t1 = c11 + c2;
    e = t1 - c11;
    t2 = ((c2 - e) + (c11 - (t1 - e))) + c21 + dda.lo * ddb.lo;

    // The result is t1 + t2, after normalization.
    ddc.hi = t1 + t2;
    ddc.lo = t2 - (ddc.hi - t1);

    return ddc;
}


ddouble ddmuld(const ddouble& dda, double db) {
    const double split = 134217729.0;
    double a1, a2, b1, b2, cona, conb, c11, c21, c2, e, t1, t2;

    // This splits dda.hi and db into high-order and low-order words.
    cona = dda.hi * split;
    conb = db * split;
    a1 = cona - (cona - dda.hi);
    b1 = conb - (conb - db);
    a2 = dda.hi - a1;
    b2 = db - b1;

    // Multiply dda.hi * db using Dekker's method.
    c11 = dda.hi * db;
    c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2;

    // Compute dda.lo * db (only high-order word is needed).
    c2 = dda.lo * db;

    // Compute (c11, c21) + c2 using Knuth's trick.
    t1 = c11 + c2;
    e = t1 - c11;
    t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;

    // The result is t1 + t2, after normalization.
    ddouble ddc;
    ddc.hi = t1 + t2;
    ddc.lo = t2 - (ddc.hi - t1);

    return ddc;
}

ddouble ddmuldd(const double da, const double db) {
    const double split = 134217729.0;
    double cona, conb, a1, a2, b1, b2, s1, s2;
    ddouble ddc;

    // On systems with a fused multiply add, such as IBM systems, it is faster to
    // uncomment the next two lines and comment out the following lines.
    // On other systems, do the opposite.

    // s1 = da * db;
    // s2 = da * db - s1;

    // This splits da and db into high-order and low-order words.

    cona = da * split;
    conb = db * split;
    a1 = cona - (cona - da);
    b1 = conb - (conb - db);
    a2 = da - a1;
    b2 = db - b1;

    // Multiply da * db using Dekker's method.

    s1 = da * db;
    s2 = (((a1 * b1 - s1) + a1 * b2) + a2 * b1) + a2 * b2;

    ddc.hi = s1;
    ddc.lo = s2;

    return ddc;
}

ddouble ddneg(const ddouble& a)
{
    ddouble b;
    b.hi = -a.hi;
    b.lo = -a.lo;
    return b;
}

const double T105 = pow(2.0, 105);
const double T52 = pow(2.0, 52);
const ddouble CON = {T105, T52};

ddouble ddnint(const ddouble& a) {
    ddouble b = {0.0, 0.0};
    ddouble s0;

    // Check if `a` is zero
    if (a.hi == 0.0) {
        return b;
    }

    // Check if `a` is too large
    if (a.hi >= T105) {
        std::cerr << "*** DDNINT: Argument is too large.\n";
        return b; // Unreachable if ddabrt terminates, added for safety
    }

    // Perform the rounding logic
    if (a.hi > 0.0) {
        s0 = ddadd(a, CON);
        b = ddsub(s0, CON);
    } else {
        s0 = ddsub(a, CON);
        b = ddadd(s0, CON);
    }

    return b;
}


ddouble ddnpwr(const ddouble& a, int n) {
    const double cl2 = 1.4426950408889633;
    ddouble s0, s1, s2 = {1.0, 0.0};
    double t1;
    int nn, mn, kn, kk;

    if (a.hi == 0.0) {
        if (n >= 0) {
            return {0.0, 0.0};
        } else {
            throw std::runtime_error("*** DDNPWR: Argument is zero and N is negative or zero.");
        }
    }

    nn = abs(n);
    if (nn == 0) {
        return {1.0, 0.0};
    } else if (nn == 1) {
        if(n > 0)
            return a;
        else
            // calculate reciprocal
            return dddiv(ddouble(1.0,0.0),a);
    } else if (nn == 2) {
        if(n > 0)
            return ddmul(a, a);
        else
            // calculate reciprocal
            return dddiv(ddouble(1.0,0.0),ddmul(a, a));
    }

    // Determine the least integer MN such that 2 ^ MN > NN.
    t1 = static_cast<double>(nn);
    mn = static_cast<int>(cl2 * log(t1) + 1.0 + 1.0e-14);

    s0 = a;
    kn = nn;

    // Compute B ^ N using the binary rule for exponentiation.
    for (int j = 1; j <= mn; ++j) {
        kk = kn / 2;
        if (kn != 2 * kk) {
            s1 = ddmul(s2, s0);
            s2 = s1;
        }
        kn = kk;
        if (j < mn) {
            s1 = ddmul(s0, s0);
            s0 = s1;
        }
    }

    // Compute reciprocal if N is negative.
    if (n < 0) {
        s1 = {1.0, 0.0};
        s0 = dddiv(s1, s2);
        s2 = s0;
    }

    return s2;
}


ddouble ddpolyr(const int n, const std::vector<ddouble>& a, const ddouble& x0) {
    const double eps = 1.0e-29;
    std::vector<ddouble> ad(n+1);
    ddouble x = x0;
    double dt1;
    
    for (int i = 0; i < n; ++i) {
        dt1 = static_cast<double>(i + 1);
        ad[i] = ddmuld(a[i+1], dt1);
    }
    ad[n].hi = 0.0;
    ad[n].lo = 0.0;
    
    for (int it = 0; it < 20; ++it) {
        ddouble t1 = {0.0, 0.0};
        ddouble t2 = {0.0, 0.0};
        ddouble t3 = {1.0, 0.0};
        ddouble t4, t5;
        
        for (int i = 0; i <= n; ++i) {
            t4 = ddmul(a[i], t3);
            t5 = ddadd(t1, t4);
            t1 = t5;
            t4 = ddmul(ad[i], t3);
            t5 = ddadd(t2, t4);
            t2 = t5;
            t4 = ddmul(t3, x);
            t3 = t4;
        }
        
        t3 = dddiv(t1, t2);
        t4 = ddsub(x, t3);
        x = t4;
        if (std::abs(t3.hi) <= eps) {
            return x;
        }
    }
    
    std::cerr << "DDROOT: failed to converge." << std::endl;
    return ddouble();
}

ddouble ddpower(const ddouble& a, const ddouble& b) {
    if (a.hi <= 0.0) {
        std::cerr << "DDPOWER: A <= 0" << std::endl;
        return ddouble();
    }
    ddouble t1 = ddlog(a);
    ddouble t2 = ddmul(t1, b);
    ddouble c = ddexp(t2);
    return c;
}



ddouble ddsqrt(const ddouble& a) {

    if (a.hi == 0.0) {
        return ddouble();
    }

    double t1 = 1.0 / sqrt(a.hi);
    double t2 = a.hi * t1;
    ddouble s0 = ddmuldd(t2, t2);
    ddouble s1 = ddsub(a, s0);
    double t3 = 0.5 * s1.hi * t1;
    s0.hi = t2;
    s0.lo = 0.0;
    s1.hi = t3;
    s1.lo = 0.0;
    ddouble b = ddadd(s0, s1);

    return b;
}

ddouble ddsub(const ddouble& dda, const ddouble& ddb) {
    ddouble ddc;
    double e, t1, t2;

    // Compute dda - ddb using Knuth's trick.
    t1 = dda.hi - ddb.hi;
    e = t1 - dda.hi;
    t2 = ((-ddb.hi - e) + (dda.hi - (t1 - e))) + dda.lo - ddb.lo;

    // The result is t1 + t2, after normalization.
    ddc.hi = t1 + t2;
    ddc.lo = t2 - (ddc.hi - t1);

    return ddc;
}

