program tquaddd

!   This program tests the tanh-sinh, exp-sinh and sinh-sinh quadrature algorithms
!   on a suite of test numerical integration problems.
!   This version uses the DDFUN double-double precision package.

!   David H Bailey   7 Jan 2023

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2023 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!   This program demonstrates three quadrature routines:

!   quadts:  Implements the tanh-sinh quadrature scheme of Takahashi and Mori,
!     for functions on a finite interval such as (0,1).
!   quades:  Implements the exp-sinh scheme, for functions on a
!     semi-infinite interval such as (0, infinity).
!   quadss:  Implements the sinh-sinh scheme, for functions on the entire
!     real line.

!   These schemes have two very desirable properties:  (a) the cost of
!   precomputing abscissa-eight pairs only increases linearly (instead of
!   quadratically, as with Gaussian quadrature), and (b) the schemes often
!   work well when the function has a blow-up singularity at one or both
!   endpoints.  Each of these routines proceeds level-by-level, with the
!   computational cost (and, often, the accuracy) approximately doubling with
!   each iteration, until a target accuracy (approx. 70-digit precision),
!   based on an accuracy estimate computed by the program, is achieved.  

!   These schemes are all variants on the tanh-sinh scheme, which is described here:

!   David H. Bailey, Xiaoye S. Li and Karthik Jeyabalan, "A comparison of
!   three high-precision quadrature schemes," Experimental Mathematics, 
!   vol. 14 (2005), no. 3, pg. 317-329, preprint available at
!   http://www.davidhbailey.com/dhbpapers/quadrature-em.pdf.

!   The function to be integrated must be defined in an external function
!   subprogram (see samples below), and the name of the function must be
!   included in a "type (dd_real)" and an "external" statement.  Prior to
!   calling the quadrature routine, the corresponding initialization routine
!   must be called to initialize the abscissa and weight arrays:  initqts,
!   initqes and initqss, corresponding to the quadrature routines quadts,
!   quades and quadss, respectively.

!   All of these routines are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   Here are some specific instructions for the individual quadrature schemes:

!   quadts:  

!   For this routine, the endpoints are specified in the variables x1 and
!   x2. In the initialization routine for quadts, abscissa-weight pairs are
!   calculated until the weights are smaller than 10^(neps2).

!   quades:

!   For this routine, it is assumed that the left endpoint is x1, and the
!   right endpoint is +infinity.

!   quadss:

!   No endpoints are specified here -- the integral is performed over the
!   entire real line.

!   These inputs are set in the parameter statement below:

!   ndp1   Primary ("low") precision in digits; this is the target accuracy
!          of quadrature results; default = 32.
!   ndp2   Secondary ("high") precision in digits; default = ndp1.
!   neps1  Log10 of the primary tolerance. By default, neps1 = - ndp1.
!   neps2  Log10 of the secondary tolerance. By default, neps2 = - ndp2.
!   nq1    Max number of phases in quadrature routine; adding 1 increases
!          (possibly doubles) the number of accurate digits in the result,
!          but also roughly doubles the run time. nq1 must be at least 3.
!   nq2    Space parameter for wk and xk arrays in the calling program.  By
!          default it is set to 12 * 2^nq1. Increase nq2 if directed by a 
!          message produced in initqts, inites or initss.

!   The endpoints x1 and x2 are set in executable statements below.

use ddmodule
implicit none
integer, parameter:: m1 = 40, m2 = 32, ndp1 = 32, ndp2 = 2 * ndp1, neps1 = -ndp1, &
  neps2 = -ndp2, nq1 = 11, nq2 = 12 * 2 ** nq1
real (kind (0.d0)), external:: second
type (dd_real), external:: quades, quadss, quadts, fun01, fun02, fun03, fun04, &
  fun05, fun06, fun07, fun08, fun09, fun10, fun11, fun12, fun13, fun14, &
  fun15, fun16, fun17, fun18
integer i
real (kind (0.d0)) tm0, tm1, tm2
real (ddknd) err, errmx
type (dd_real) one, t1, t2, t3, t4, wkes(-1:nq2), xkes(-1:nq2), &
  wkss(-1:nq2), xkss(-1:nq2), wkts(-1:nq2), xkts(-1:nq2), zero
type (dd_real) mppic, mpl02, x1, x2

!   Compute pi and log(2) to high precision (nwds2 words).

one = ddreal (1.d0)
zero = ddreal (0.d0)
mppic = ddpi ()
mpl02 = ddlog2 ()
errmx = zero

write (6, 1) nq1, ndp1, ndp2, neps1, neps2
1 format ('Quadts test:  Quadlevel =',i6/'Digits1 =',i6,'  Digits2 =',i6, &
  '  Epsilon1 =',i6,'  Epsilon2 =',i6)

!   Initialize quadrature tables wk and xk (weights and abscissas).

tm0 = second ()
call initqes (nq1, nq2, neps2, wkes, xkes)
call initqss (nq1, nq2, neps2, wkss, xkss)
call initqts (nq1, nq2, neps2, wkts, xkts)
tm1 = second ()
tm2 = tm1 - tm0
write (6, 2) tm1 - tm0
2 format ('Quadrature initialization completed: cpu time =',f12.6)

!   Begin quadrature tests.

write (6, 11)
11 format (/'Continuous functions on finite intervals:'//&
  'Problem 1: Int_0^1 t*log(1+t) dt = 1/4'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadts (fun01, x1, x2, nq1, nq2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
3 format ('Quadrature completed: CPU time =',f12.6/'Result =')
call ddwrite (6, m1, m2, t1)
t2 = ddreal (0.25d0)
err = t2 - t1
errmx = max (abs (err), errmx)
write (6, 4) err
4 format ('Actual error =',1p,d15.6)

write (6, 12)
12 format (/'Problem 2: Int_0^1 t^2*arctan(t) dt = (pi - 2 + 2*log(2))/12'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadts (fun02, x1, x2, nq1, nq2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call ddwrite (6, m1, m2, t1)
t2 = (mppic - 2.d0 + 2.d0 * mpl02) / 12.d0
err = t2 - t1
errmx = max (abs (err), errmx)
write (6, 4) err

write (6, 13)
13 format (/'Problem 3: Int_0^(pi/2) e^t*cos(t) dt = 1/2*(e^(pi/2) - 1)'/)
x1 = zero
x2 = 0.5d0 * mppic
tm0 = second ()
t1 = quadts (fun03, x1, x2, nq1, nq2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call ddwrite (6, m1, m2, t1)
t2 = 0.5d0 * (exp (0.5d0 * mppic) - 1.d0)
err = t2 - t1
errmx = max (abs (err), errmx)
write (6, 4) err

write (6, 14)
14 format (/ &
  'Problem 4: Int_0^1 arctan(sqrt(2+t^2))/((1+t^2)sqrt(2+t^2)) dt = 5*Pi^2/96'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadts (fun04, x1, x2, nq1, nq2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call ddwrite (6, m1, m2, t1)
t2 = 5.d0 * mppic**2 / 96.d0
err = t2 - t1
errmx = max (abs (err), errmx)
write (6, 4) err

write (6, 15)
15 format (/&
  'Continuous functions on finite intervals, but non-diff at an endpoint:'// &
  'Problem 5: Int_0^1 sqrt(t)*log(t) dt = -4/9'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadts (fun05, x1, x2, nq1, nq2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call ddwrite (6, m1, m2, t1)
t2 = ddreal (-4.d0) / 9.d0
err = t2 - t1
errmx = max (abs (err), errmx)
write (6, 4) err

write (6, 16)
16 format (/'Problem 6: Int_0^1 sqrt(1-t^2) dt = pi/4'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadts (fun06, x1, x2, nq1, nq2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call ddwrite (6, m1, m2, t1)
t2 = 0.25d0 * mppic
err = t2 - t1
errmx = max (abs (err), errmx)
write (6, 4) err

write (6, 17)
17 format (/&
  'Functions on finite intervals with integrable singularity at an endpoint:'//&
  'Problem 7: Int_0^1 sqrt(t)/sqrt(1-t^2) dt = 2*sqrt(pi)*gamma(3/4)/gamma(1/4)'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadts (fun07, x1, x2, nq1, nq2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call ddwrite (6, m1, m2, t1)

! t2 = 2.d0 * sqrt (ddreal (mppic)) * gamma (ddreal (0.75d0)) &
!  / gamma (ddreal (0.25d0))

t3 = '1.2254167024651776451290983033628905268512392481080706112301189382898228884267984'
t4 = '3.6256099082219083119306851558676720029951676828800654674333779995699192435387291'
t2 = 2.d0 * sqrt (ddreal (mppic)) * t3 / t4
err = t2 - t1
errmx = max (abs (err), errmx)
write (6, 4) err

write (6, 18)
18 format (/'Problem 8: Int_0^1 log(t)^2 dt = 2'/)
x1 = zero
x2 = one
tm0 = second ()
t1 = quadts (fun08, x1, x2, nq1, nq2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call ddwrite (6, m1, m2, t1)
t2 = ddreal (2.d0)
err = t2 - t1
errmx = max (abs (err), errmx)
write (6, 4) err

write (6, 19)
19 format (/'Problem 9: Int_0^(pi/2) log(cos(t)) dt = -pi*log(2)/2'/)
x1 = zero
x2 = 0.5d0 * mppic
tm0 = second ()
t1 = quadts (fun09, x1, x2, nq1, nq2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call ddwrite (6, m1, m2, t1)
t2 = -0.5d0 * mppic * mpl02
err = t2 - t1
errmx = max (abs (err), errmx)
write (6, 4) err

write (6, 20)
20 format (/'Problem 10: Int_0^(pi/2) sqrt(tan(t)) dt = pi*sqrt(2)/2'/)
x1 = zero
x2 = 0.5d0 * mppic
tm0 = second ()
t1 = quadts (fun10, x1, x2, nq1, nq2, neps1, wkts, xkts)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call ddwrite (6, m1, m2, t1)
t2 = 0.5d0 * mppic * sqrt (ddreal (2.d0))
err = t2 - t1
errmx = max (abs (err), errmx)
write (6, 4) err

write (6, 21)
21 format (/&
  'Functions on a semi-infinite interval:'//&
  'Problem 11: Int_0^inf 1/(1+t^2) dt = pi/2'/)
x1 = zero
tm0 = second ()
t1 = quades (fun11, x1, nq1, nq2, neps1, wkes, xkes)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call ddwrite (6, m1, m2, t1)
t2 = 0.5d0 * mppic
err = t2 - t1
errmx = max (abs (err), errmx)
write (6, 4) err

write (6, 22)
22 format (/'Problem 12: Int_0^inf e^(-t)/sqrt(t) dt = sqrt(pi)'/)
x1 = zero
tm0 = second ()
t1 = quades (fun12, x1, nq1, nq2, neps1, wkes, xkes)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call ddwrite (6, m1, m2, t1)
t2 = sqrt (mppic)
err = t2 - t1
errmx = max (abs (err), errmx)
write (6, 4) err

write (6, 23)
23 format (/'Problem 13: Int_0^inf e^(-t^2/2) dt = sqrt(pi/2)'/)
x1 = zero
tm0 = second ()
t1 = quades (fun13, x1, nq1, nq2, neps1, wkes, xkes)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call ddwrite (6, m1, m2, t1)
t2 = sqrt (0.5d0 * mppic)
err = t2 - t1
errmx = max (abs (err), errmx)
write (6, 4) err

write (6, 24)
24 format (/'Problem 14: Int_0^inf e^(-t)*cos(t) dt = 1/2'/)
x1 = zero
tm0 = second ()
t1 = quades (fun14, x1, nq1, nq2, neps1, wkes, xkes)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call ddwrite (6, m1, m2, t1)
t2 = ddreal (0.5d0)
err = t2 - t1
errmx = max (abs (err), errmx)
write (6, 4) err

write (6, 25)
25 format (/ &
   'Functions on the entire real line:'// &
   'Problem 15: Int_-inf^inf 1/(1+t^2) dt = Pi'/)
tm0 = second ()
t1 = quadss (fun15, nq1, nq2, neps1, wkss, xkss)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call ddwrite (6, m1, m2, t1)
t2 = mppic
err = t2 - t1
errmx = max (abs (err), errmx)
write (6, 4) err

write (6, 26)
26 format (/'Problem 16: Int_-inf^inf 1/(1+t^4) dt = Pi/Sqrt(2)'/)
tm0 = second ()
t1 = quadss (fun16, nq1, nq2, neps1, wkss, xkss)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call ddwrite (6, m1, m2, t1)
t2 = mppic / sqrt (ddreal (2.d0))
err = t2 - t1
errmx = max (abs (err), errmx)
write (6, 4) err

write (6, 27)
27 format (/'Problem 17: Int_-inf^inf e^(-t^2/2) dt = sqrt (2*Pi)'/)
tm0 = second ()
t1 = quadss (fun17, nq1, nq2, neps1, wkss, xkss)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call ddwrite (6, m1, m2, t1)
t2 = sqrt (2.d0 * mppic)
err = t2 - t1
errmx = max (abs (err), errmx)
write (6, 4) err

write (6, 28)
28 format (/'Problem 18: Int_-inf^inf e^(-t^2/2) cos(t) dt = sqrt (2*Pi/e)'/)
tm0 = second ()
t1 = quadss (fun18, nq1, nq2, neps1, wkss, xkss)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
write (6, 3) tm1 - tm0
call ddwrite (6, m1, m2, t1)
t2 = sqrt (2.d0 * mppic / exp (ddreal (1.d0)))
err = t2 - t1
errmx = max (abs (err), errmx)
write (6, 4) err

write (6, 91) tm2, errmx
91 format (/'Total CPU time =',f12.6/'Max abs error =',1p,d15.6)
if (errmx < 10.q0**(3+neps1/2)) then
  write (6, '(a)') 'ALL TESTS PASSED'
else
  write (6, '(a)') 'ONE OR MORE TESTS FAILED'
endif

stop
end program tquaddd

type (dd_real) function fun01 (t)

!   fun01(t) = t * log(1+t)

use ddmodule
implicit none
type (dd_real), intent(in):: t
type (dd_real) t1

t1 = ddreal (t)
fun01 = t1 * log (1.d0 + t1)
return
end function fun01

type (dd_real) function fun02 (t)

!   fun02(t) = t^2 * arctan(t)

use ddmodule
implicit none
type (dd_real), intent(in):: t
type (dd_real) t1

t1 = ddreal (t)
fun02 = t1 ** 2 * atan (t1)
return
end function fun02

type (dd_real) function fun03 (t)

!   fun03(t) = e^t * cos(t)

use ddmodule
implicit none
type (dd_real), intent(in):: t
type (dd_real) t1

t1 = ddreal (t)
fun03 = exp(t1) * cos(t1)
return
end function fun03

type (dd_real) function fun04 (t)

!   fun04(t) = arctan(sqrt(2+t^2))/((1+t^2)sqrt(2+t^2))

use ddmodule
implicit none
type (dd_real), intent(in):: t
type (dd_real) t1, t2

t1 = ddreal (t)
t2 = sqrt (2.d0 + t1**2)
fun04 = atan(t2) / ((1.d0 + t1**2) * t2)
return
end function fun04

type (dd_real) function fun05 (t)

!    fun05(t) = sqrt(t)*log(t)

use ddmodule
implicit none
type (dd_real), intent(in):: t
type (dd_real) t1

t1 = ddreal (t)
fun05 = sqrt (t1) * log (t1)
return
end function fun05

type (dd_real) function fun06 (t)

!    fun06(t) = sqrt(1-t^2)

use ddmodule
implicit none
type (dd_real), intent(in):: t
type (dd_real) t1, t2

t1 = ddreal (t)
t2 = ddreal (1.d0 - t**2)
fun06 = sqrt (t2)
return
end function fun06

type (dd_real) function fun07 (t)

!   fun07(t) = sqrt (t) / sqrt(1-t^2)

use ddmodule
implicit none
type (dd_real), intent(in):: t
type (dd_real) t1, t2

!   The subtraction to compute t2 must be performed using high precision
!   (nwds2), but after the subtraction its low precision value is fine.

t1 = ddreal (t)
t2 = ddreal (1.d0 - t)
fun07 = sqrt (t1) / sqrt (t2 * (1.d0 + t1))
return
end function fun07

type (dd_real) function fun08 (t)

!   fun08(t) = log(t)^2

use ddmodule
implicit none
type (dd_real), intent(in):: t
type (dd_real) t1

t1 = ddreal (t)
fun08 = log (t1) ** 2
return
end function fun08

type (dd_real) function fun09 (t)

!   fun09(t) = log (cos (t))

use ddmodule
implicit none
type (dd_real), intent(in):: t
type (dd_real) pi, t1, t2, t3, t4

t1 = ddreal (t)
pi = ddpi ()
t3 = ddreal (0.25d0 * pi)
t2 = ddreal (0.5d0 * pi - t)

if (t1 < t3) then
  t4 = cos (t1)
else
  t4 = sin (t2)
endif
fun09 = log (t4)

return
end function fun09

type (dd_real) function fun10 (t)

!   fun10(t) = sqrt(tan(t))

use ddmodule
implicit none
type (dd_real), intent(in):: t
type (dd_real) pi, t1, t2, t3

t1 = ddreal (t)
pi = ddpi ()
t3 = ddreal (0.25d0 * pi)
t2 = ddreal (0.5d0 * pi - t)

if (t1 < t3) then
  fun10 = sqrt (tan (t1))
else
  fun10 = 1.d0 / sqrt (tan (t2))
endif
return
end function fun10

type (dd_real) function fun11 (t)

!   1/(1 + t^2) on (0, Inf).

use ddmodule
implicit none
type (dd_real), intent(in):: t
type (dd_real) t1

t1 = ddreal (t)
fun11 = 1.d0 / (1.d0 + t1 ** 2)
return
end function fun11

type (dd_real) function fun12 (t)

!   e^(-t)/sqrt(t) on (0, inf).

use ddmodule
implicit none
type (dd_real), intent(in):: t
type (dd_real) t1

!   The subtraction to compute t2 must be performed using high precision
!   (nwds2), but after the subtraction its low precision value is fine.

t1 = ddreal (t)
fun12 = exp (-t1) / sqrt (t1)
return
end function fun12

type (dd_real) function fun13 (t)

!   e^(-t^2/2) on (0, inf).

use ddmodule
implicit none
type (dd_real), intent(in):: t
type (dd_real) t1

t1 = ddreal (t)
fun13 = exp (-0.5d0 * t1 ** 2)
return
end function fun13

type (dd_real) function fun14 (t)

!   e^(-t) cos(t) on (0, inf).

use ddmodule
implicit none
type (dd_real), intent(in):: t
type (dd_real) t1

t1 = ddreal (t)
fun14 = exp (-t1) * cos (t1)
return
end function fun14

type (dd_real) function fun15 (t)
use ddmodule
implicit none
type (dd_real), intent(in):: t

fun15 = 1.d0 / (1.d0 + t**2)
return
end function fun15

type (dd_real) function fun16 (t)
use ddmodule
implicit none
type (dd_real), intent(in):: t

fun16 = 1.d0 / (1.d0 + t**4)
return
end function fun16

type (dd_real) function fun17 (t)
use ddmodule
implicit none
type (dd_real), intent(in):: t

fun17 = exp (-0.5d0 * t**2)
return
end function fun17

type (dd_real) function fun18 (t)
use ddmodule
implicit none
type (dd_real), intent(in):: t

fun18 = exp (-0.5d0 * t**2) * cos (t)
return
end function fun18

subroutine initqes (nq1, nq2, neps2, wk, xk)

!   This subroutine initializes the quadrature arays xk and wk using the
!   function x(t) = exp (pi/2*sinh(t)). The argument nq2 is the space
!   allocated for wk and xk in the calling program. By default it is set to 
!   12 * 2^nq1.  Increase nq2 if directed by a message produced below.
!   neps2 controls termination of the loop below, which ends when 
!   wk(k) * 10^(neps2) > 1. If quadts outputs the message "Increase 
!   Quadlevel" or "Terms too large", adjust nq1 and neps2 as necessary.

!   Both initqes and quades are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   David H Bailey   7 Jan 2023

use ddmodule
implicit none
integer, intent(in):: nq1, nq2, neps2
type (dd_real), intent(out)::  wk(-1:nq2), xk(-1:nq2)
integer, parameter:: iprint = 1024, ndebug = 2
integer i, j, k, k1
type (dd_real) eps2, h, p2, t1, t2, t3, t4, u1, u2
 
write (6, 1)
1 format ('initqes: Exp-sinh quadrature initialization')

eps2 = ddreal (10.d0) ** neps2
p2 = 0.5d0 * ddpi ()
h = ddreal (0.5d0 ** nq1)
wk(-1) = ddreal (real (nq1, ddknd))

do k = 0, nq2
  if (ndebug >= 2 .and. mod (k, iprint) == 0) write (6, *) k, nq2
    t1 = ddreal (real (k, ddknd) * h)

!   xk(k) = exp (u1)
!   wk(k) = exp (u1) * u2
!   where u1 = pi/2 * sinh (t1) and u2 = pi/2 * cosh (t1)

  t2 = exp (t1)
  u1 = 0.5d0 * p2 * (t2 - 1.d0 / t2)
  u2 = 0.5d0 * p2 * (t2 + 1.d0 / t2)
  xk(k) = exp (u1)
  wk(k) = xk(k) * u2

  if (wk(k) * eps2 > 1.d0) goto 100
enddo

write (6, 2) nq2
2 format ('initqes: Table space parameter is too small; value =',i8)
stop

100 continue

xk(-1) = ddreal (real (k, ddknd))
if (ndebug >= 2) then
  write (6, 3) k
3 format ('initqes: Table spaced used =',i8)
endif

return
end subroutine initqes

type (dd_real) function quades (fun, x1, nq1, nq2, neps1, wk, xk)

!   This routine computes the integral of the function fun on the interval
!   (x1, inf) with a target tolerance of 10^neps1.  The quadrature level is
!   progressively increased (approximately doubling the work with each level)
!   until level nq1 has been performed or the target tolerance has been met.
!   nq2 is the size of the wk and xk arrays, which must first be initialized
!   by calling initqes.  If quades outpues the message "Increase Quadlevel"
!   or "Terms too large", adjust nq1 and neps2 as necessary in the call to
!   initqes.

!   Both initqes and quades are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   David H Bailey  7 Jan 2023

use ddmodule
implicit none
type (dd_real), external:: fun
type (dd_real), intent(in):: x1
integer, intent(in):: nq1, nq2, neps1
type (dd_real), intent(in)::  wk(-1:nq2), xk(-1:nq2)
integer, parameter:: izx = 5, m1 = 40, m2 = 32, ndebug = 2
real (ddknd), external:: dplog10q
integer i, ierror, ip(0:100), iz1, iz2, j, k, k1, k2, n, nds, nqq1
logical log1
real (ddknd) d1, d2, d3, d4
type (dd_real) c10, eps1, eps2, epsilon1, err, h, &
  tsum, s1, s2, s3, t1, t2, t3, tw1, tw2, twi1, twi2, twmx 
type (dd_real) xx1, xx2

epsilon1 = ddreal (10.d0) ** neps1
tsum = ddreal (0.d0)
s1 = ddreal (0.d0)
s2 = ddreal (0.d0)
h = ddreal (1.d0)
c10 = ddreal (10.d0)

if (wk(-1) < real (nq1, ddknd)) then
  write (6, 1) nq1
1 format ('quades: quadrature arrays have not been initialized; nq1 =',i6)
  goto 140
endif
nqq1 = dble (wk(-1))
n = dble (xk(-1))

do k = 0, nqq1
  ip(k) = 2 ** k
enddo

do k = 1, nq1
  h = 0.5d0 * h
  s3 = s2
  s2 = s1
  k1 = ip(nqq1-k)
  k2 = ip(nqq1-k+1)
  iz1 = 0
  iz2 = 0
  twmx = ddreal (0.d0)

!   Evaluate function at level k in x, avoiding unnecessary computation.

  do i = 0, n, k1
    if (mod (i, k2) /= 0 .or. k == 1) then

!   These next few lines, which scale the abscissas, must be performed in
!   high precision (nwds2 words) to ensure full accuracy in the quadrature
!   results, even though the abscissas xk(i) were computed in low precision.

      xx1 = x1 + ddreal (xk(i))
      xx2 = x1 + 1.d0 / ddreal (xk(i))
      log1 = xx1 > x1
  
!   The remaining computations are performed in low precision (nwds1 words).

      if (iz1 < izx) then
        t1 = fun (xk(i))
        tw1 = t1 * wk(i)
        twi1 = abs (tw1)
        if (twi1 < epsilon1) then
          iz1 = iz1 + 1
        else
          iz1 = 0
        endif
      else
        t1 = ddreal (0.d0)
        tw1 = ddreal (0.d0)
      endif

      if (i > 0 .and. log1 .and. iz2 < izx) then
        t2 = fun (xx2)
        tw2 = t2 * wk(i) / xk(i)**2
        twi2 = abs (tw2)
        if (twi2 < epsilon1) then
          iz2 = iz2 + 1
        else
          iz2 = 0
        endif
      else
        t2 = ddreal (0.d0)
        tw2 = ddreal (0.d0)
      endif

      tsum = tsum + tw1 + tw2
      twmx = max (twmx, abs (tw1), abs (tw2))
    endif
  enddo

!   Compute s1 = current integral approximation and err = error estimate.
!   Tsum is the sum of all tw1 and tw2 from the loop above.
!   Twmx is the largest absolute value of tw1 and tw2 from the loop above.
!   Twi1 and twi2 are the final nonzero values of abs(tw1) and abs(tw2).

  s1 =  h * tsum
  eps1 = twmx * epsilon1
  eps2 = max (twi1, twi2)
  d1 = dplog10q (abs (s1 - s2))
  d2 = dplog10q (abs (s1 - s3))
  d3 = dplog10q (eps1) - 1.d0
  d4 = dplog10q (eps2) - 1.d0

  if (k <= 2) then
    err = ddreal (1.d0)
  elseif (d1 .eq. -999999.d0) then
    err = ddreal (0.d0)
  else
    err = c10 ** nint (min (0.d0, max (d1 ** 2 / d2, 2.d0 * d1, d3, d4)))
  endif

!   Output current integral approximation and error estimate, to 60 digits.

  if (ndebug >= 2) then
    write (6, 2) k, nq1, nint (dplog10q (abs (err)))
2   format ('quades: Iteration',i3,' of',i3,'; est error = 10^',i7, &
      '; approx value =')
    call ddwrite (6, m1, m2, s1)
  endif

  if (k >= 3 .and. iz1 == 0 .and. iz2 == 0) then
    write (6, 3)
3   format ('quades: Terms too large -- adjust neps2 in call to initqss.')
    goto 140
  endif

  if (k >= 3 .and. err < eps1) then
    write (6, 4) nint (dplog10q (abs (err)))
4   format ('quades: Estimated error = 10^',i7)
    goto 140
  endif

  if (k >= 3 .and. err < eps2) then
    write (6, 5) nint (dplog10q (abs (err)))
5   format ('quades: Estimated error = 10^',i7/&
    'Adjust nq1 and neps2 in initqss for greater accuracy.')
    goto 140
  endif
enddo

140 continue

quades = s1
return
end function quades

subroutine initqss (nq1, nq2, neps2, wk, xk)

!   This subroutine initializes the quadrature arays xk and wk using the
!   function x(t) = sinh (pi/2*sinh(t)).  The argument nq2 is the space
!   allocated for wk and xk in the calling program.  By default it is set to 
!   12 * 2^nq1.  Increase nq2 if directed by a message produced below.
!   neps2 controls termination of the loop below, which ends when 
!   wk(k) * 10^(neps2) > 1. If quadss outputs the message "Increase 
!   Quadlevel" or "Terms too large", adjust nq1 and neps2 as necessary.

!   Both initqss and quadss are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   David H Bailey   7 Jan 2023

use ddmodule
implicit none
integer, intent(in):: nq1, nq2, neps2
type (dd_real), intent(out):: wk(-1:nq2), xk(-1:nq2)
integer, parameter:: iprint = 1024, ndebug = 2
integer i, ierror, j, k, k1
type (dd_real) eps2, h, p2, t1, t2, t3, t4, u1, u2

write (6, 1)
1 format ('initqss: Sinh-sinh quadrature initialization')

eps2 = ddreal (10.d0) ** neps2
p2 = 0.5d0 * ddpi ()
h = ddreal (0.5d0 ** nq1)
wk(-1) = ddreal (real (nq1, ddknd))

do k = 0, nq2
  if (ndebug >= 2 .and. mod (k, iprint) == 0) write (6, *) k, nq2
    t1 = ddreal (real (k, ddknd) * h)

!   xk(k) = sinh (u1)
!   wk(k) = cosh (u1) * u2
!   where u1 = pi/2 * sinh (t1) and u2 = pi/2 * cosh (t1)

  t2 = exp (t1)
  u1 = 0.5d0 * p2 * (t2 - 1.d0 / t2)
  u2 = 0.5d0 * p2 * (t2 + 1.d0 / t2)
  t3 = exp (u1)
  xk(k) = 0.5d0 * (t3 - 1.d0 / t3)
  wk(k) = 0.5d0 * (t3 + 1.d0 / t3) * u2

  if (wk(k) * eps2 > 1.d0) goto 100
enddo

write (6, 2) nq2
2 format ('initqss: Table space parameter is too small; value =',i8)
stop

100 continue

xk(-1) = ddreal (real (k, ddknd))
if (ndebug >= 2) then
  write (6, 3) k
3 format ('initqss: Table spaced used =',i8)
endif

return
end subroutine initqss

type (dd_real) function quadss (fun, nq1, nq2, neps1, wk, xk)

!   This routine computes the integral of the function fun on the interval
!   (-inf, inf) with a target tolerance of 10^neps1.  The quadrature level is
!   progressively increased (approximately doubling the work with each level)
!   until level nq1 has been performed or the target tolerance has been met.
!   nq2 is the size of the wk and xk arrays, which must first be initialized
!   by calling initqss.  If quadss outputs the message "Increase Quadlevel"
!   or "Terms too large", adjust nq1 and neps2 as necessary in the call to
!   initqss.

!   Both initqss and quadss are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   David H Bailey  7 Jan 2023

use ddmodule
implicit none
type (dd_real), external:: fun
integer, intent(in):: nq1, nq2, neps1
type (dd_real), intent(in):: wk(-1:nq2), xk(-1:nq2)
integer, parameter:: izx = 5, m1 = 40, m2 = 32, ndebug = 2
real (ddknd), external:: dplog10q
integer i, ierror, ip(0:100), iz1, iz2, j, k, k1, k2, n, nds, nqq1
real (ddknd) d1, d2, d3, d4
type (dd_real) c10, eps1, eps2, epsilon1, err, h, &
  tsum, s1, s2, s3, t1, t2, t3, tw1, tw2, twi1, twi2, twmx

epsilon1 = ddreal (10.d0) ** neps1
tsum = ddreal (0.d0)
s1 = ddreal (0.d0)
s2 = ddreal (0.d0)
h = ddreal (1.d0)
c10 = ddreal (10.d0)

if (wk(-1) < real (nq1, ddknd)) then
  write (6, 1) nq1
1 format ('quadss: quadrature arrays have not been initialized; nq1 =',i6)
  goto 140
endif
nqq1 = dble (wk(-1))
n = dble (xk(-1))

do k = 0, nqq1
  ip(k) = 2 ** k
enddo

do k = 1, nq1
  h = 0.5d0 * h
  s3 = s2
  s2 = s1
  k1 = ip(nqq1-k)
  k2 = ip(nqq1-k+1)
  iz1 = 0
  iz2 = 0
  twmx = ddreal (0.d0)

!   Evaluate function at level k in x, avoiding unnecessary computation.

  do i = 0, n, k1
    if (mod (i, k2) /= 0 .or. k == 1) then
      if (iz1 < izx) then
        t1 = fun (xk(i))
        tw1 = t1 * wk(i)
        twi1 = abs (tw1)
        if (twi1 < epsilon1) then
          iz1 = iz1 + 1
        else
          iz1 = 0
        endif
      else
        t1 = ddreal (0.d0)
        tw1 = ddreal (0.d0)
      endif

      if (i > 0 .and. iz2 < izx) then
        t2 = fun (-xk(i))
        tw2 = t2 * wk(i)
        twi2 = abs (tw2)
        if (twi2 < epsilon1) then
          iz2 = iz2 + 1
        else
          iz2 = 0
        endif
      else
        t2 = ddreal (0.d0)
        tw2 = ddreal (0.d0)
      endif

      tsum = tsum + tw1 + tw2
      twmx = max (twmx, abs (tw1), abs (tw2))
    endif
  enddo

!   Compute s1 = current integral approximation and err = error estimate.
!   Tsum is the sum of all tw1 and tw2 from the loop above.
!   Twmx is the largest absolute value of tw1 and tw2 from the loop above.
!   Twi1 and twi2 are the final nonzero values of abs(tw1) and abs(tw2).

  s1 =  h * tsum
  eps1 = twmx * epsilon1
  eps2 = max (twi1, twi2)
  d1 = dplog10q (abs (s1 - s2))
  d2 = dplog10q (abs (s1 - s3))
  d3 = dplog10q (eps1) - 1.d0
  d4 = dplog10q (eps2) - 1.d0

  if (k <= 2) then
    err = ddreal (1.d0)
  elseif (d1 .eq. -999999.d0) then
    err = ddreal (0.d0)
  else
    err = c10 ** nint (min (0.d0, max (d1 ** 2 / d2, 2.d0 * d1, d3, d4)))
  endif

!   Output current integral approximation and error estimate, to 60 digits.

  if (ndebug >= 2) then
    write (6, 2) k, nq1, nint (dplog10q (abs (err)))
2   format ('quadss: Iteration',i3,' of',i3,'; est error = 10^',i7, &
      '; approx value =')
    call ddwrite (6, m1, m2, s1)
  endif

  if (k >= 3 .and. iz1 == 0 .and. iz2 == 0) then
    write (6, 3)
3   format ('quadss: Terms too large -- adjust neps2 in call to initqss.')
    goto 140
  endif

  if (k >= 3 .and. err < eps1) then
    write (6, 4) nint (dplog10q (abs (err)))
4   format ('quadss: Estimated error = 10^',i7)
    goto 140
  endif

  if (k >= 3 .and. err < eps2) then
    write (6, 5) nint (dplog10q (abs (err)))
5   format ('quadss: Estimated error = 10^',i7/&
    'Adjust nq1 and neps2 in initqss for greater accuracy.')
    goto 140
  endif
enddo

140 continue

quadss = s1
return
end function quadss

subroutine initqts (nq1, nq2, neps2, wk, xk)

!   This subroutine initializes the quadrature arays xk and wk using the
!   function x(t) = tanh (pi/2*sinh(t)).  The argument nq2 is the space
!   allocated for wk and xk in the calling program.  By default it is set to 
!   12 * 2^nq1.  Increase nq2 if directed by a message produced below.
!   Upon completion, wk(-1) = nq1, and xk(-1) = n, the maximum space parameter
!   for these arrays.  In other words, the arrays occupy (wk(i), i = -1 to n)
!   and (xk(i), i = -1 to n), where n = xk(-1).   The array x_k contains 
!   1 minus the abscissas; the wk array contains the weights at these abscissas.
!   If quadts outputs the message "Increase Quadlevel" or "Terms too large",
!   adjust nq1 and neps2 as necessary in the call to initqss.

!   Both initqts and quadts are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   David H Bailey   7 Jan 2023

use ddmodule
implicit none
integer, intent(in):: nq1, nq2, neps2
type (dd_real), intent(out):: wk(-1:nq2), xk(-1:nq2)
integer, parameter:: iprint = 1024, ndebug = 2
integer i, ierror, j, k, k1
type (dd_real) eps2, h, p2, t1, t2, t3, t4, u1, u2

write (6, 1)
1 format ('initqts: Tanh-sinh quadrature initialization')

eps2 = ddreal (10.d0) ** neps2
p2 = 0.5d0 * ddpi ()
h = ddreal (0.5d0 ** nq1)
wk(-1) = ddreal (real (nq1, ddknd))

do k = 0, nq2
  if (ndebug >= 2 .and. mod (k, iprint) == 0) write (6, *) k, nq2
  t1 = ddreal (real (k, ddknd) * h)

!   xk(k) = 1 - tanh (u1) = 1 /(e^u1 * cosh (u1))
!   wk(k) = u2 / cosh (u1)^2
!   where u1 = pi/2 * cosh (t1), u2 = pi/2 * sinh (t1)

  t2 = exp (t1)
  u1 = 0.5d0 * p2 * (t2 + 1.d0 / t2)
  u2 = 0.5d0 * p2 * (t2 - 1.d0 / t2)
  t3 = exp (u2)
  t4 = 0.5d0 * (t3 + 1.d0 / t3)
  xk(k) = 1.d0 / (t3 * t4)
  wk(k) = u1 / t4 ** 2

  if (wk(k) < eps2) goto 100
enddo

write (6, 2) nq2
2 format ('initqts: Table space parameter is too small; value =',i8)
stop

100 continue

xk(-1) = ddreal (real (k, ddknd))
if (ndebug >= 2) then
  write (6, 3) k
3 format ('initqts: Table spaced used =',i8)
endif

return
end subroutine initqts

type (dd_real) function quadts (fun, x1, x2, nq1, nq2, neps1, wk, xk)

!   This routine computes the integral of the function fun on the interval
!   (x1, x2) with a target tolerance of 10^neps1.  The quadrature level is
!   progressively increased (approximately doubling the work with each level)
!   until level nq1 has been performed or the target tolerance has been met.
!   nq2 is the size of the wk and xk arrays, which must first be initialized
!   by calling initqts. The function fun is not evaluated at the endpoints
!   x1 and x2.  If quadts outputs the message "Increase Quadlevel" or "Terms
!   too large", adjust nq1 and neps2 as necessary in the call to initqss.

!   Both initqts and quadts are 100% THREAD SAFE -- all requisite parameters
!   and arrays are passed through subroutine arguments. 

!   David H Bailey  7 Jan 2023

use ddmodule
implicit none
type (dd_real), external:: fun
type (dd_real), intent(in):: x1, x2
integer, intent(in):: nq1, nq2, neps1
type (dd_real), intent(in):: wk(-1:nq2), xk(-1:nq2)
integer, parameter:: izx = 5, m1 = 40, m2 = 32, ndebug = 2
real (ddknd), external:: dplog10q
integer i, ierror, ip(0:100), iz1, iz2, j, k, k1, k2, n, nds, nqq1
logical log1, log2
real (ddknd) d1, d2, d3, d4
type (dd_real) c10, eps1, eps2, epsilon1, err, h, &
  tsum, s1, s2, s3, t1, t2, t3, tw1, tw2, twi1, twi2, twmx
type (dd_real) ax, bx, xki, xt1, xx1, xx2

!  These two lines are performed in high precision (nwds2 words).

ax = 0.5d0 * (x2 - x1)
bx = 0.5d0 * (x2 + x1)

!  The remaining initialization is performed in low precision (nwds1 words).

epsilon1 = ddreal (10.d0) ** neps1
tsum = ddreal (0.d0)
s1 = ddreal (0.d0)
s2 = ddreal (0.d0)
h = ddreal (1.d0)
c10 = ddreal (10.d0)

if (wk(-1) < real (nq1, ddknd)) then
  write (6, 1) nq1
1 format ('quadts: Quadrature arrays have not been initialized; nq1 =',i6)
  goto 140
endif
nqq1 = dble (wk(-1))
n = dble (xk(-1))

do k = 0, nqq1
  ip(k) = 2 ** k
enddo

do k = 1, nq1
  h = 0.5d0 * h
  s3 = s2
  s2 = s1
  k1 = ip(nqq1-k)
  k2 = ip(nqq1-k+1)
  iz1 = 0
  iz2 = 0
  twmx = ddreal (0.d0)

!   Evaluate function at level k in x, avoiding unnecessary computation.

  do i = 0, n, k1  
    if (mod (i, k2) /= 0 .or. k == 1) then

!   These next few lines, which scale the abscissas, must be performed in
!   high precision (nwds2 words) to ensure full accuracy in the quadrature
!   results, even though the abscissas xk(i) were computed in low precision.

      xki = xk(i)
      xt1 = 1.d0 - ddreal (xki)
      xx1 = - ax * xt1 + bx
      xx2 = ax * xt1 + bx
      log1 = xx1 > x1
      log2 = xx2 < x2      

!   The remaining computations are performed in low precision (nwds1 words).

      if (log1 .and. iz1 < izx) then
        t1 = fun (xx1)
        tw1 = t1 * wk(i)
        twi1 = abs (tw1)
        if (twi1 < epsilon1) then
          iz1 = iz1 + 1
        else
          iz1 = 0
        endif
      else
        t1 = ddreal (0.d0)
        tw1 = ddreal (0.d0)
      endif

      if (i > 0 .and. log2 .and. iz2 < izx) then
        t2 = fun (xx2)
        tw2 = t2 * wk(i)
        twi2 = abs (tw2)
        if (twi2 < epsilon1) then
          iz2 = iz2 + 1
        else
          iz2 = 0
        endif
      else
        t2 = ddreal (0.d0)
        tw2 = ddreal (0.d0)
      endif

      tsum = tsum + tw1 + tw2
      twmx = max (twmx, abs (tw1), abs (tw2))
    endif
  enddo

!   Compute s1 = current integral approximation and err = error estimate.
!   Tsum is the sum of all tw1 and tw2 from the loop above.
!   Twmx is the largest absolute value of tw1 and tw2 from the loop above.
!   Twi1 and twi2 are the final nonzero values of abs(tw1) and abs(tw2).

  s1 =  ddreal (ax) * h * tsum
  eps1 = twmx * epsilon1
  eps2 = max (twi1, twi2)
  d1 = dplog10q (abs (s1 - s2))
  d2 = dplog10q (abs (s1 - s3))
  d3 = dplog10q (eps1) - 1.d0
  d4 = dplog10q (eps2) - 1.d0

  if (k <= 2) then
    err = ddreal (1.d0)
  elseif (d1 .eq. -999999.d0) then
    err = ddreal (0.d0)
  else
    err = c10 ** nint (min (0.d0, max (d1 ** 2 / d2, 2.d0 * d1, d3, d4)))
  endif

!   Output current integral approximation and error estimate, to 60 digits.

  if (ndebug >= 2) then
    write (6, 2) k, nq1, nint (dplog10q (abs (err)))
2   format ('quadts: Iteration',i3,' of',i3,'; est error = 10^',i7, &
      '; approx value =')
    call ddwrite (6, m1, m2, s1)
  endif

  if (k >= 3 .and. iz1 == 0 .and. iz2 == 0) then
    write (6, 3)
3   format ('quadts: Terms too large -- adjust neps2 in call to initqss.')
    goto 140
  endif

  if (k >= 3 .and. err < eps1) then
    write (6, 4) nint (dplog10q (abs (err)))
4   format ('quadts: Estimated error = 10^',i7)
    goto 140
  endif

  if (k >= 3 .and. err < eps2) then
    write (6, 5) nint (dplog10q (abs (err)))
5   format ('quadts: Estimated error = 10^',i7/&
    'Adjust nq1 and neps2 in initqts for greater accuracy.')
    goto 140
  endif
enddo

140 continue

quadts = s1
return
end function quadts

real (ddknd) function dplog10q (a)

!   For input DD value a, this routine returns a DP approximation to log10 (a).

use ddmodule
implicit none
type (dd_real), intent(in):: a
integer ia
real (ddknd) da

da = dble (a)
if (da .eq. 0.d0) then
  dplog10q = -999999.d0
else
  dplog10q = log10 (abs (da)) + ia * log10 (2.d0)
endif

100 continue
return
end function dplog10q

