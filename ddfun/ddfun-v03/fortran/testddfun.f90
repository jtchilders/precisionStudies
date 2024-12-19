program testddfun

!   This briefly tests most individual DDFUN operations and functions
!   (including all mixed mode arithmetic and comparison operations), by
!   comparing each result with results in the reference file testddfun.ref.txt.
!   This is clearly not an exhaustive test, but it often detects software bugs
!   and compiler issues.

!   David H Bailey  27 Feb 2023

use ddmodule
implicit none
integer, parameter:: ndp = 32, neps = -31, nfile = 11, n1 = 40, n2 = 32, &
  np = 2, nq = 3, nbe = 2.d0 * ndp, nrr = 10
integer i, i1, i2, i3, i4
logical l1, l2, l3, l4, l5, l6, l7, l8, l9, l10
character(120) chr120
character(32) chr32
character(1) chr1(120)
real (ddknd) d1, d2, d3, d4, e1, e2
complex (ddknd) dc1, dc2, dc3, dc4, ec1, ec2
real (ddknd) eps, err, errmx
type (dd_real) aa(np), bb(nq), be(nbe), rr(nrr), t1, t2, t3, t4, zero, one
type (dd_complex) z1, z2, z3, z4

! t1 = '0.5772156649015328606065120900824024310422'
! write (6, '(1p,2d28.18)') t1
! stop

open (nfile, file = 'testddfun.ref.txt')
rewind nfile

write (6, '(a)') 'DDFUN quick check of operations and functions'

!   Define a few sample data values.

eps = ddreal (10.d0) ** neps
zero = 0.d0
one = 1.d0
errmx = zero
aa(1) = 0.75d0
aa(2) = 1.25d0
bb(1) = 1.d0
bb(2) = 1.5d0
bb(3) = 2.d0
call ddberne (nbe, be)
call polylog_ini (-nrr, rr)

t1 = ddpi ()
t2 = - ddlog2 ()
d1 = t1
d2 = t2
e1 = 3141.d0 / 8192.d0
e2 = 6931.d0 / 8192.d0
z1 = ddcmplx (0.5d0 * ddpi (), exp (ddreal (0.5d0)))
z2 = ddcmplx (- sqrt (ddpi ()), cos (ddreal (1.d0)))
dc1 = z1
dc2 = z2
ec1 = cmplx (e1, e2, ddknd)
ec2 = cmplx (-e2, e1, ddknd)
i1 = 5
i2 = -3
chr120 = ' '

write (6, '(/a/)') 'Test data:'
write (6, '(a)') 't1 = pi:'
call ddwrite (6, n1, n2, t1)
call checkdd (nfile, 5, t1, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 't2 = -log(2):'
call ddwrite (6, n1, n2, t2)
call checkdd (nfile, 1, t2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'z1 = (0.5*pi, exp(0.5)):'
call ddwrite (6, n1, n2, z1)
call checkddc (nfile, 1, z1, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'z2 = (-Gamma(0.5), Cos(1)):'
call ddwrite (6, n1, n2, z2)
call checkddc (nfile, 1, z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'e1 = 3141/8192:'
write (6, '(1p,d45.35)') e1
write (6, '(a)') 'e2 = 6931 / 8192:'
write (6, '(1p,d45.35)') e2
write (6, '(a)')  'ec1 = (3141/8192, 6931/8192)'
write (6, '(1p,d45.35)') ec1
write (6, '(a)') 'ec2 = (6931/8192, 3141/8192):'
write (6, '(1p,d45.35)') ec2

write (6, '(/a/)') 'Real data operations:'

write (6, '(a)') 'addition: t1+t2 ='
call ddwrite (6, n1, n2, t1 + t2)
call checkdd (nfile, 14, t1 + t2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'addition: t1+e2 ='
call ddwrite (6, n1, n2, t1 + e2)
call checkdd (nfile, 1, t1 + e2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'addition: e1+t2 ='
call ddwrite (6, n1, n2, e1 + t2)
call checkdd (nfile, 1, e1 + t2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: t1-t2 ='
call ddwrite (6, n1, n2, t1 - t2)
call checkdd (nfile, 1, t1 - t2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: t1-e2 ='
call ddwrite (6, n1, n2, t1 - e2)
call checkdd (nfile, 1, t1 - e2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: e1-t2 ='
call ddwrite (6, n1, n2, e1 - t2)
call checkdd (nfile, 1, e1 - t2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: t1*t2 ='
call ddwrite (6, n1, n2, t1 * t2)
call checkdd (nfile, 1, t1 * t2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: t1*e2 ='
call ddwrite (6, n1, n2, t1 * e2)
call checkdd (nfile, 1, t1 * e2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: e1*t2 ='
call ddwrite (6, n1, n2, e1 * t2)
call checkdd (nfile, 1, e1 * t2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: t1/t2 ='
call ddwrite (6, n1, n2, t1 / t2)
call checkdd (nfile, 1, t1 / t2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: t1/e2 ='
call ddwrite (6, n1, n2, t1 / e2)
call checkdd (nfile, 1, t1 / e2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: e1/t2 ='
call ddwrite (6, n1, n2, e1 / t2)
call checkdd (nfile, 1, e1 / t2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'exponentiation: t1**i1 ='
call ddwrite (6, n1, n2, t1 ** i1)
call checkdd (nfile, 1, t1 ** i1, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'exponentiation: t1**t2 ='
call ddwrite (6, n1, n2, t1 ** t2)
call checkdd (nfile, 1, t1 ** t2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'equal test: t1 == t2, e1 == t2, t1 == e2'
l1 = t1 == t2; l2 = d1 == d2; l3 = e1 == t2
write (6, '(6l4)') l1, l2, l3
if (l1 .eqv. .false. .and. l2 .eqv. .false. .and. l3 .eqv. .false.) then
  err = zero
else
  err = one
endif
errmx = max (err, errmx)

write (6, '(a)') 'not-equal test: t1 /= t2, e1 /= t2, t1 =/ e2'
l1 = t1 /= t2; l2 = d1 /= d2; l3 = e1 /= t2
write (6, '(6l4)') l1, l2, l3
if (l1 .eqv. .true. .and. l2 .eqv. .true. .and. l3 .eqv. .true.) then
  err = zero
else
  err = one
endif
errmx = max (err, errmx)

write (6, '(a)') 'less-than-or-equal test: t1 <= t2, e1 <= t2, t1 <= e2'
l1 = t1 <= t2; l2 = d1 <= d2; l3 = e1 <= t2
write (6, '(6l4)') l1, l2, l3
if (l1 .eqv. .false. .and. l2 .eqv. .false. .and. l3 .eqv. .false.) then
  err = zero
else
  err = one
endif
errmx = max (err, errmx)

write (6, '(a)') 'greater-than-or-equal test: t1 >= t2, e1 >= t2, t1 >= e2'
l1 = t1 >= t2; l2 = d1 >= d2; l3 = e1 >= t2
write (6, '(6l4)') l1, l2, l3
if (l1 .eqv. .true. .and. l2 .eqv. .true. .and. l3 .eqv. .true.) then
  err = zero
else
  err = one
endif
errmx = max (err, errmx)

write (6, '(a)') 'less-than test: t1 < t2, e1 < t2, t1 < e2'
l1 = t1 < t2; l2 = d1 < d2; l3 = e1 < t2
write (6, '(6l4)') l1, l2, l3
if (l1 .eqv. .false. .and. l2 .eqv. .false. .and. l3 .eqv. .false.) then
  err = zero
else
  err = one
endif
errmx = max (err, errmx)

write (6, '(a)') 'greater-than test: t1 > t2, e1 > t2, t1 > e2'
l1 = t1 > t2; l2 = d1 > d2; l3 = e1 > t2
write (6, '(6l4)') l1, l2, l3
if (l1 .eqv. .true. .and. l2 .eqv. .true. .and. l3 .eqv. .true.) then
  err = zero
else
  err = one
endif
errmx = max (err, errmx)

write (6, '(a)') 'abs(t2) ='
call ddwrite (6, n1, n2, abs (t2))
call checkdd (nfile, 13, abs (t2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'acos(t2) ='
call ddwrite (6, n1, n2, acos (t2))
call checkdd (nfile, 1, acos (t2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'acosh(t1) ='
call ddwrite (6, n1, n2, acosh (t1))
call checkdd (nfile, 1, acosh (t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'agm(t1,abs(t2)) ='
call ddwrite (6, n1, n2, agm (t1, abs (t2)))
call checkdd (nfile, 1, agm (t1, abs (t2)), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'aint(t1) ='
call ddwrite (6, n1, n2, aint (t1))
call checkdd (nfile, 1, aint (t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'anint(t1) ='
call ddwrite (6, n1, n2, anint (t1))
call checkdd (nfile, 1, anint (t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'asin(t2) ='
call ddwrite (6, n1, n2, asin (t2))
call checkdd (nfile, 1, asin (t2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'asinh(t1) ='
call ddwrite (6, n1, n2, asinh (t1))
call checkdd (nfile, 1, asinh (t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'atan(t1) ='
call ddwrite (6, n1, n2, atan (t1))
call checkdd (nfile, 1, atan (t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'atan2(t1,t2) ='
call ddwrite (6, n1, n2, atan2 (t1,t2))
call checkdd (nfile, 1, atan2 (t1,t2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'atanh(t2) ='
call ddwrite (6, n1, n2, atanh (t2))
call checkdd (nfile, 1, atanh (t2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_i(t1,-t2) ='
call ddwrite (6, n1, n2, bessel_i (t1, -t2))
call checkdd (nfile, 1, bessel_i (t1, -t2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_in(3,t1) ='
call ddwrite (6, n1, n2, bessel_in (3, t1))
call checkdd (nfile, 1, bessel_in (3, t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_in(3,50*t1) ='
call ddwrite (6, n1, n2, bessel_in (3, 50.d0*t1))
call checkdd (nfile, 1, bessel_in (3, 50.d0*t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_j(t1,-t2) ='
call ddwrite (6, n1, n2, bessel_j (t1, -t2))
call checkdd (nfile, 1, bessel_j (t1, -t2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_jn(3,t1) ='
call ddwrite (6, n1, n2, bessel_jn (3, t1))
call checkdd (nfile, 1, bessel_jn (3, t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_jn(3,50*t1) ='
call ddwrite (6, n1, n2, bessel_jn (3, 50.d0*t1))
call checkdd (nfile, 1, bessel_jn (3, 50.d0*t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_k(t1,-t2) ='
call ddwrite (6, n1, n2, bessel_k (t1, -t2))
call checkdd (nfile, 1, bessel_k (t1, -t2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_kn(3,t1) ='
call ddwrite (6, n1, n2, bessel_kn (3, t1))
call checkdd (nfile, 1, bessel_kn (3, t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_kn(3,50*t1) ='
call ddwrite (6, n1, n2, bessel_kn (3, 50.d0*t1))
call checkdd (nfile, 1, bessel_kn (3, 50.d0*t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_y(t1,-t2) ='
call ddwrite (6, n1, n2, bessel_y (t1, -t2))
call checkdd (nfile, 1, bessel_y (t1, -t2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_yn(3,t1) ='
call ddwrite (6, n1, n2, bessel_yn (3, t1))
call checkdd (nfile, 1, bessel_yn (3, t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'bessel_yn(3,50*t1) ='
call ddwrite (6, n1, n2, bessel_yn (3, 50.d0*t1))
call checkdd (nfile, 1, bessel_yn (3, 50.d0*t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'cos(t2) ='
call ddwrite (6, n1, n2, cos (t2))
call checkdd (nfile, 1, cos (t2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'cosh(t1) ='
call ddwrite (6, n1, n2, cosh (t1))
call checkdd (nfile, 1, cosh (t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'digamma_be(nbe,be,t1) ='
call ddwrite (6, n1, n2, digamma_be (nbe,be,t1))
call checkdd (nfile, 1, digamma_be (nbe,be,t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'ddcssh(t1) ='
call ddcssh (t1, t3, t4)
call ddwrite (6, n1, n2, t3, t4)
call checkdd (nfile, 1, t3, eps, err)
errmx = max (err, errmx)
call checkdd (nfile, 0, t4, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'ddcssn(t2) ='
call ddcssn (t2, t3, t4)
call ddwrite (6, n1, n2, t3, t4)
call checkdd (nfile, 1, t3, eps, err)
errmx = max (err, errmx)
call checkdd (nfile, 0, t4, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'ddeform(t1,n1,n2,chr1) ='
call ddeform (t1, n1, n2, chr1)
write (6, '(80a1)') (chr1(i), i = 1, n1)
do i = 1, n1; chr120(i:i) = chr1(i); enddo; t4 = chr120(1:n1)
call checkdd (nfile, 1, t4, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'ddfform(t1,n1,n2,chr1) ='
call ddfform (t1, n1, n2, chr1)
write (6, '(80a1)') (chr1(i), i = 1, n1)
do i = 1, n1; chr120(i:i) = chr1(i); enddo; t4 = chr120(1:n1)
call checkdd (nfile, 1, t4, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'ddnrt(t1,i1) ='
call ddwrite (6, n1, n2, ddnrt (t1,i1))
call checkdd (nfile, 1, ddnrt(t1,i1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'ddreal (chr120) ='
call ddwrite (6, n1, n2, ddreal (chr120))
call checkdd (nfile, 1, t1, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'erf(t1) ='
call ddwrite (6, n1, n2, erf (t1))
call checkdd (nfile, 1, erf (t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'erfc(t1) ='
call ddwrite (6, n1, n2, erfc (t1))
call checkdd (nfile, 1, erfc (t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'exp(t1) ='
call ddwrite (6, n1, n2, exp (t1))
call checkdd (nfile, 1, exp (t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'expint(t1) ='
call ddwrite (6, n1, n2, expint (t1))
call checkdd (nfile, 1, expint (t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'gamma(t1) ='
call ddwrite (6, n1, n2, gamma (t1))
call checkdd (nfile, 1, gamma (t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'hurwitz_zetan(3,1/t1) ='
call ddwrite (6, n1, n2, hurwitz_zetan (3,1.d0/t1))
call checkdd (nfile, 1, hurwitz_zetan (3,1.d0/t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'hurwitz_zetan_be(nbe,be,5,t1) ='
call ddwrite (6, n1, n2, hurwitz_zetan_be (nbe,be,5,t1))
call checkdd (nfile, 1, hurwitz_zetan_be (nbe,be,5,t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'hypergeom_pfq(np,nq,aa,bb,t1) ='
call ddwrite (6, n1, n2, hypergeom_pfq (np,nq,aa,bb,t1))
call checkdd (nfile, 1, hypergeom_pfq (np,nq,aa,bb,t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'incgamma(t1,t2+2) ='
call ddwrite (6, n1, n2, incgamma (t1, t2+2.d0))
call checkdd (nfile, 1, incgamma (t1, t2+2.d0), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'incgamma(-3,t1) ='
t3 = -3.d0
call ddwrite (6, n1, n2, incgamma (t3, t1))
call checkdd (nfile, 1, incgamma (t3, t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'log(t1) ='
call ddwrite (6, n1, n2, log (t1))
call checkdd (nfile, 1, log (t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'log10(t1) ='
call ddwrite (6, n1, n2, log10 (t1))
call checkdd (nfile, 1, log10 (t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'max(t1,t2) ='
call ddwrite (6, n1, n2, max (t1,t2))
call checkdd (nfile, 1, max (t1,t2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'min(t1,t2) ='
call ddwrite (6, n1, n2, min (t1,t2))
call checkdd (nfile, 1, min (t1,t2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'mod(t1,t2) ='
call ddwrite (6, n1, n2, mod(t1,t2))
call checkdd (nfile, 1, mod (t1, t2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'polygamma(3,1/t1) ='
call ddwrite (6, n1, n2, polygamma (3,1.d0/t1))
call checkdd (nfile, 1, polygamma (3,1.d0/t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'polygamma_be(nbe,be,5,t1) ='
call ddwrite (6, n1, n2, polygamma_be (nbe,be,5,t1))
call checkdd (nfile, 1, polygamma_be (nbe,be,5,t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'polylog_neg (-10, rr, -t1) ='
call ddwrite (6, n1, n2, polylog_neg (-10, rr, -t1))
call checkdd (nfile, 1, polylog_neg (-10, rr, -t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'polylog_pos (10, 1/t1) ='
call ddwrite (6, n1, n2, polylog_pos (10, 1.d0/t1))
call checkdd (nfile, 1, polylog_pos (10, 1.d0/t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'sign(t1,t2) ='
call ddwrite (6, n1, n2, sign (t1,t2))
call checkdd (nfile, 1, sign (t1,t2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'sin(t2) ='
call ddwrite (6, n1, n2, sin (t2))
call checkdd (nfile, 1, sin (t2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'sinh(t1) ='
call ddwrite (6, n1, n2, sinh (t1))
call checkdd (nfile, 1, sinh (t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'sqrt(t1) ='
call ddwrite (6, n1, n2, sqrt (t1))
call checkdd (nfile, 1, sqrt (t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'struve_hn(10,t1) ='
call ddwrite (6, n1, n2, struve_hn (10, t1))
call checkdd (nfile, 1, struve_hn (10,t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'tan(t2) ='
call ddwrite (6, n1, n2, tan (t2))
call checkdd (nfile, 1, tan (t2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'tanh(t1) ='
call ddwrite (6, n1, n2, tanh (t1))
call checkdd (nfile, 1, tanh (t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'zeta(t1) ='
call ddwrite (6, n1, n2, zeta (t1))
call checkdd (nfile, 1, zeta (t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'zeta_be(nbe,be,t1) ='
call ddwrite (6, n1, n2, zeta_be (nbe,be,t1))
call checkdd (nfile, 1, zeta_be (nbe,be,t1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'zeta_int(10) ='
call ddwrite (6, n1, n2, zeta_int (10))
call checkdd (nfile, 1, zeta_int (10), eps, err)
errmx = max (err, errmx)

write (6, '(/a/)') 'Complex data operations:'

write (6, '(a)') 'addition: z1+z2 ='
call ddwrite (6, n1, n2, z1 + z2)
call checkddc (nfile, 4, z1 + z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'addition: z1+e2 ='
call ddwrite (6, n1, n2, z1 + e2)
call checkddc (nfile, 1, z1 + e2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'addition: e1+z2 ='
call ddwrite (6, n1, n2, e1 + z2)
call checkddc (nfile, 1, e1 + z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'addition: z1+ec2 ='
call ddwrite (6, n1, n2, z1 + ec2)
call checkddc (nfile, 1, z1 + ec2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'addition: ec1+z2 ='
call ddwrite (6, n1, n2, ec1 + z2)
call checkddc (nfile, 1, ec1 + z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'addition: z1+t2 ='
call ddwrite (6, n1, n2, z1 + t2)
call checkddc (nfile, 1, z1 + t2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'addition: t1+z2 ='
call ddwrite (6, n1, n2, t1 + z2)
call checkddc (nfile, 1, t1 + z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: z1-z2 ='
call ddwrite (6, n1, n2, z1 - z2)
call checkddc (nfile, 1, z1 - z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: z1-e2 ='
call ddwrite (6, n1, n2, z1 - e2)
call checkddc (nfile, 1, z1 - e2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: e1-z2 ='
call ddwrite (6, n1, n2, e1 - z2)
call checkddc (nfile, 1, e1 - z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: z1-ec2 ='
call ddwrite (6, n1, n2, z1 - ec2)
call checkddc (nfile, 1, z1 - ec2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: ec1-z2 ='
call ddwrite (6, n1, n2, ec1 - z2)
call checkddc (nfile, 1, ec1 - z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: z1-t2 ='
call ddwrite (6, n1, n2, z1 - t2)
call checkddc (nfile, 1, z1 - t2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'subtraction: t1-z2 ='
call ddwrite (6, n1, n2, t1 - z2)
call checkddc (nfile, 1, t1 - z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: z1*z2 ='
call ddwrite (6, n1, n2, z1 * z2)
call checkddc (nfile, 1, z1 * z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: z1*e2 ='
call ddwrite (6, n1, n2, z1 * e2)
call checkddc (nfile, 1, z1 * e2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: e1*z2 ='
call ddwrite (6, n1, n2, e1 * z2)
call checkddc (nfile, 1, e1 * z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: z1*ec2 ='
call ddwrite (6, n1, n2, z1 * ec2)
call checkddc (nfile, 1, z1 * ec2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: ec1*z2 ='
call ddwrite (6, n1, n2, ec1 * z2)
call checkddc (nfile, 1, ec1 * z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: z1*t2 ='
call ddwrite (6, n1, n2, z1 * t2)
call checkddc (nfile, 1, z1 * t2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: t1*z2 ='
call ddwrite (6, n1, n2, t1 * z2)
call checkddc (nfile, 1, t1 * z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: z1*e2 ='
call ddwrite (6, n1, n2, z1 * e2)
call checkddc (nfile, 1, z1 * e2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'multiplication: e1*z2 ='
call ddwrite (6, n1, n2, e1 * z2)
call checkddc (nfile, 1, e1 * z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: z1/z2 ='
call ddwrite (6, n1, n2, z1 / z2)
call checkddc (nfile, 1, z1 / z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: z1/e2 ='
call ddwrite (6, n1, n2, z1 / e2)
call checkddc (nfile, 1, z1 / e2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: e1/z2 ='
call ddwrite (6, n1, n2, e1 / z2)
call checkddc (nfile, 1, e1 / z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: z1/ec2 ='
call ddwrite (6, n1, n2, z1 / ec2)
call checkddc (nfile, 1, z1 / ec2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: ec1/z2 ='
call ddwrite (6, n1, n2, ec1 / z2)
call checkddc (nfile, 1, ec1 / z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: z1/t2 ='
call ddwrite (6, n1, n2, z1 / t2)
call checkddc (nfile, 1, z1 / t2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'division: t1/z2 ='
call ddwrite (6, n1, n2, t1 / z2)
call checkddc (nfile, 1, t1 / z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'exponentiation: z1**i1 ='
call ddwrite (6, n1, n2, z1 ** i1)
call checkddc (nfile, 1, z1 ** i1, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'exponentiation: z1**z2 ='
call ddwrite (6, n1, n2, z1 ** z2)
call checkddc (nfile, 1, z1 ** z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'exponentiation: t1**z2 ='
call ddwrite (6, n1, n2, t1 ** z2)
call checkddc (nfile, 1, t1 ** z2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'exponentiation: z1**t2 ='
call ddwrite (6, n1, n2, z1 ** t2)
call checkddc (nfile, 1, z1 ** t2, eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'equal test: z1 == z2, e1 == z2, z1 == e2, ec1 == z2, z1 == ec2'
l1 = z1 == z2; l2 = dc1 == dc2; l3 = e1 == z2; l4 = e1 == dc2; l5 = z1 == e2
write (6, '(10l4)') l1, l2, l3, l4, l5
if (l1 .eqv. .false. .and. l2 .eqv. .false. .and. l3 .eqv. .false. &
  .and. l4 .eqv. .false. .and. l5 .eqv. .false.) then
  err = zero
else
  err = one
endif
errmx = max (err, errmx)

write (6, '(a)') 'not-equal test: z1 /= z2, e1 /= z2, z1 /= e2, ec1 /= z2, z1 /= ec2'
l1 = z1 /= z2; l2 = dc1 /= dc2; l3 = e1 /= z2; l4 = e1 /= dc2; l5 = z1 /= e2
write (6, '(10l4)') l1, l2, l3, l4, l5
if (l1 .eqv. .true. .and. l2 .eqv. .true. .and. l3 .eqv. .true. &
  .and. l4 .eqv. .true. .and. l5 .eqv. .true.) then
  err = zero
else
  err = one
endif
errmx = max (err, errmx)

write (6, '(a)') 'abs(z2) ='
call ddwrite (6, n1, n2, abs (z2))
call checkdd (nfile, 5, abs (z2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'aimag(z1) ='
call ddwrite (6, n1, n2, aimag (z1))
call checkdd (nfile, 1, aimag (z1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'conjg(z1) ='
call ddwrite (6, n1, n2, conjg (z1))
call checkddc (nfile, 1, conjg (z1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'cos(z2) ='
call ddwrite (6, n1, n2, cos (z2))
call checkddc (nfile, 1, cos (z2), eps, err)
errmx = max (err, errmx)

! write (6, '(a)') 'qcmplx(z1) ='
! write (6, '(1p,d45.35)') qcmplx (z1)
! call checkddc (nfile, 1, qcmplx (z1), eps, err)
! errmx = max (err, errmx)

write (6, '(a)') 'exp(z1) ='
call ddwrite (6, n1, n2, exp (z1))
call checkddc (nfile, 1, exp (z1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'log(z1) ='
call ddwrite (6, n1, n2, log (z1))
call checkddc (nfile, 1, log (z1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'ddcmplx(t1,t2) ='
call ddwrite (6, n1, n2, ddcmplx (t1, t2))
call checkddc (nfile, 1, ddcmplx (t1, t2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'ddreal(z1) ='
call ddwrite (6, n1, n2, ddreal (z1))
call checkdd (nfile, 1, ddreal (z1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'sin(z2) ='
call ddwrite (6, n1, n2, sin (z2))
call checkddc (nfile, 1, sin (z2), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'sqrt(z1) ='
call ddwrite (6, n1, n2, sqrt (z1))
call checkddc (nfile, 1, sqrt (z1), eps, err)
errmx = max (err, errmx)

write (6, '(a)') 'sqrt(z2) ='
call ddwrite (6, n1, n2, sqrt (z2))
call checkddc (nfile, 1, sqrt (z2), eps, err)
errmx = max (err, errmx)

write (6, 9) errmx
9 format (/'Max relative error =',1p,d15.6)

if (abs (errmx) < eps) then
  write (6, '(a)') 'ALL TESTS PASSED'
else
  write (6, '(a)') 'ONE OR MORE TESTS FAILED'
endif

stop
end program testddfun

subroutine checkdd (nfile, i1, t1, eps, err)
use ddmodule
implicit none
integer, intent(in):: nfile, i1
type (dd_real), intent(in):: t1
real (ddknd), intent(in):: eps
real (ddknd), intent(out):: err
type (dd_real) t2
integer i
character(120) c1

do i = 1, i1
  read (nfile, '(a)') c1
enddo

call ddread (nfile, t2)
err = abs ((t1 - t2) / t2)

if (abs (err) > eps) write (6, 1) err
1 format ('ERROR:',1pd15.6)

return
end subroutine checkdd

subroutine checkddc (nfile, i1, z1, eps, err)
use ddmodule
implicit none
integer, intent(in):: nfile, i1
type (dd_complex), intent(in):: z1
real (ddknd), intent(in):: eps
real (ddknd), intent(out):: err
type (dd_complex) z2
integer i
character(64) c1

do i = 1, i1
  read (nfile, '(a)') c1
enddo

call ddread (nfile, z2)
err = abs ((z1 - z2) / z2)

if (abs (err) > eps) write (6, 1) err
1 format ('ERROR:',1pd15.6)

return
end subroutine checkddc

