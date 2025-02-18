#pragma once
#include "ddouble.h"
#include "ddcomplex.h"
#include <vector>

namespace ddfun{
   
ddouble ddabs(const ddouble& );
ddouble ddacosh(const ddouble& );
ddouble ddadd(const ddouble&, const ddouble& );
ddouble ddasinh(const ddouble& );
ddouble ddatanh(const ddouble& );
void ddcsshr(const ddouble &, ddouble &, ddouble &);
void ddcssnr(const ddouble &, ddouble &, ddouble &);
ddouble dddiv(const ddouble&, const ddouble& );
ddouble dddivd(const ddouble&, const double );
ddouble ddexp(const ddouble& );
ddouble ddlog(const ddouble& );
ddouble ddmul(const ddouble&, const ddouble& );
ddouble ddmuld(const ddouble&, const double );
ddouble ddmuldd(const double, const double);
ddouble ddneg(const ddouble& );
ddouble ddnint(const ddouble& );
ddouble ddnpwr(const ddouble&, int );
ddouble ddpolyr(const int, const std::vector<ddouble>&, const ddouble& );
ddouble ddpower(const ddouble& , const ddouble& );
ddouble ddsqrt(const ddouble& );
ddouble ddsub(const ddouble&, const ddouble& );

ddcomplex ddcpwr(const ddcomplex& , const int& );
ddcomplex ddcsqrt(const ddcomplex& );


}