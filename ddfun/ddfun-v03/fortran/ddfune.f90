!*****************************************************************************

!  DDFUN: A double-double package with special functions
!  Special functions module (module DDFUNE)

!  Revision date:  16 Mar 2023

!  AUTHOR:
!    David H. Bailey
!    Lawrence Berkeley National Lab (retired)
!    Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2023 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF PACKAGE:
!    This package permits one to perform floating-point computations (real and
!    complex) to arbitrarily high numeric precision, by making only relatively
!    minor changes to existing Fortran-90 programs. All basic arithmetic
!    operations and transcendental functions are supported, together with numerous
!    special functions.

!    In addition to fast execution times, one key feature of this package is a
!    100% THREAD-SAFE design, which means that user-level applications can be
!    easily converted for parallel execution, say using a threaded parallel
!    environment such as OpenMP.

!  DOCUMENTATION:
!    A detailed description of this package, and instructions for compiling
!    and testing this program on various specific systems are included in the
!    README file accompanying this package, and, in more detail, in the
!    following technical paper:

!    David H. Bailey, "MPFUN2020: A new thread-safe arbitrary precision package,"
!    available at http://www.davidhbailey.com/dhbpapers/mpfun2020.pdf.

!  DESCRIPTION OF THIS MODULE (DDFUNE):
!    This module contains subroutines to perform special functions. Additional
!    functions will be added as they are completed.

!  NOTE ON PROGRAMMING CONVENTION FOR THIS MODULE:
!    This module is designed to facilitate easy translation (using the program
!    convmpfr.f90), for use in the MPFUN-MPFR package, and also for translation
!    (using the program convmpdq.f90), for use in the DQFUN double-quad package.

module ddfune
use ddfuna

contains

!  These routines perform simple operations on the MP data structure. Those
!  listed between !> and !>> are for MPFUN20-Fort; those between !>> and !>>>
!  are for MPFUN20-MPFR. Those between !>>> and !>>>> are for DQFUN. The
!  translator selects the proper set.
!>
! 
! subroutine mpinitwds (ra, mpnw)
! implicit none
! integer (mpiknd), intent(out):: ra(0:)
! integer, intent(in):: mpnw
! 
! ra(0) = mpnw + 6
! ra(1) = mpnw
! ra(2) = 0; ra(3) = 0; ra(4) = 0
! return
! end subroutine mpinitwds
! 
! integer function mpwprecr (ra)
! implicit none
! integer (mpiknd), intent(in):: ra(0:)
! 
! mpwprecr = ra(1)
! return
! end function mpwprecr
! 
! integer function mpspacer (ra)
! implicit none
! integer (mpiknd), intent(in):: ra(0:)
! 
! mpspacer = ra(0)
! return
! end function mpspacer
! 
!>>

! subroutine mpabrt (ier)
! implicit none
! integer, intent(in):: ier
! write (mpldb, 1) ier
! 1 format ('*** MPABRT: Execution terminated, error code =',i4)
! stop
! end subroutine mpabrt

! subroutine mpinitwds (ra, mpnw)
! implicit none
! integer (mpiknd), intent(out):: ra(0:)
! integer, intent(in):: mpnw
!
! ra(0) = mpnw + 6
! ra(1) = mpnw * mpnbt
! ra(2) = 1
! ra(3) = mpnan
! ra(4) = loc (ra(4)) + 8
! ra(mpnw+5) = 0
! return
! end subroutine mpinitwds

! subroutine mpfixlocr (ia)
! implicit none
! integer (mpiknd), intent(out):: ia(0:)
!
! ia(4) = loc (ia(4)) + 8
! return
! end subroutine

! integer function mpwprecr (ra)
! implicit none
! integer (mpiknd), intent(in):: ra(0:)
!
! mpwprecr = ra(1) / mpnbt
! return
! end function mpwprecr

! integer function mpspacer (ra)
! implicit none
! integer (mpiknd), intent(in):: ra(0:)
!
! mpspacer = ra(0)
! return
! end function mpspacer

!>>>

! integer function dqwprecr (ra)
! implicit none
! real (dqknd), intent(in):: ra(2)
!
! dqwprecr = 2
! return
! end function dqwprecr

! integer function dqspacer (ra)
! implicit none
! real (dqknd), intent(in):: ra(2)
!
! dqspacer = 2
! return
! end function dqspacer

!>>>>
integer function ddwprecr (ra)
implicit none
real (ddknd), intent(in):: ra(2)

ddwprecr = 2
return
end function ddwprecr

integer function ddspacer (ra)
implicit none
real (ddknd), intent(in):: ra(2)

ddspacer = 2
return
end function ddspacer
!>>>>>

subroutine ddberner (nb2, berne)

!   This returns the array berne, containing Bernoulli numbers indexed 2*k for
!   k = 1 to n, to mpnw words precision. This is done by first computing
!   zeta(2*k), based on the following known formulas:

!   coth (pi*x) = cosh (pi*x) / sinh (pi*x)

!            1      1 + (pi*x)^2/2! + (pi*x)^4/4! + ...
!        =  ---- * -------------------------------------
!           pi*x    1 + (pi*x)^2/3! + (pi*x)^4/5! + ...

!        = 1/(pi*x) * (1 + (pi*x)^2/3 - (pi*x)^4/45 + 2*(pi*x)^6/945 - ...)

!        = 2/(pi*x) * Sum_{k >= 1} (-1)^(k+1) * zeta(2*k) * x^{2*k}

!   The strategy is to calculate the coefficients of the series by polynomial
!   operations. Polynomial division is performed by computing the reciprocal
!   of the denominator polynomial, by a polynomial Newton iteration, as follows.
!   Let N(x) be the polynomial approximation to the numerator series; let D(x) be
!   a polynomial approximation to the numerator numerator series; and let Q_k(x)
!   be polynomial approximations to R(x) = 1/D(x). Then iterate:

!   Q_{k+1} = Q_k(x) + [1 - D(x)*Q_k(x)]*Q_k(x)

!   In these iterations, both the degree of the polynomial Q_k(x) and the
!   precision level in words are initially set to 4. When convergence is
!   achieved at this precision level, the degree is doubled, and iterations are
!   continued, etc., until the final desired degree is achieved. Then the
!   precision level is doubled and iterations are performed in a similar way,
!   until the final desired precision level is achieved. The reciprocal polynomial
!   R(x) produced by this process is then multiplied by the numerator polynomial
!   N(x) to yield an approximation to the quotient series. The even zeta values
!   are then the coefficients of this series, scaled according to the formula above.

!   Once the even integer zeta values have been computed in this way, the even
!   Bernoulli numbers are computed via the formula (for n > 0):

!   B(2*n) = (-1)^(n-1) * 2 * (2*n)! * zeta(2*n) / (2*pi)^(2*n)

!   Note: The notation in the description above is not the same as in the code below.

implicit none
integer, intent(in):: nb2
real (ddknd), intent(out):: berne(1:2,nb2)
integer, parameter:: ibz = 6, idb = 0
real (ddknd), parameter:: alog102 = 0.30102999566398119d0, pi = 3.1415926535897932385d0
integer i, i1, ic1, j, kn, n, n1, nn1
real (ddknd) d1, dd1, dd2, dd3
real (ddknd) c1(1:2,0:nb2), cp2(1:2), p1(1:2,0:nb2), &
  p2(1:2,0:nb2), q(1:2,0:nb2), q1(1:2), &
  r(1:2,0:nb2), s(1:2,0:nb2), t1(1:2), t2(1:2), &
  t3(1:2), t4(1:2), eps(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

n = nb2
ddnw1 = min (ddnw + 1, ddnwx)
ddnw2 = ddnw1

if (idb > 0) write (ddldb, 3) n
3 format ('Even Bernoulli number calculation; n =',2i6)


do i = 0, nb2
enddo

call ddmul (ddpicon, ddpicon, cp2)
call dddmc (1.d0, 0, c1(1:2,0))
call dddmc (1.d0, 0, p1(1:2,0))
call dddmc (1.d0, 0, p2(1:2,0))
call dddmc (1.d0, 0, q(1:2,0))

!   Construct numerator and denominator polynomials.

do i = 1, n
  call dddmc (0.d0, 0, c1(1:2,i))
  dd1 = 2.d0 * (i + 1) - 3.d0
  dd2 = dd1 + 1.d0
  dd3 = dd2 + 1.d0
  call ddmul (cp2, p1(1:2,i-1), t1)
  call dddivd (t1, dd1 * dd2, p1(1:2,i))
  call ddmul (cp2, p2(1:2,i-1), t1)
  call dddivd (t1, dd2 * dd3, p2(1:2,i))
  call dddmc (0.d0, 0, q(1:2,i))
enddo

kn = 4
ddnw2 = min (4, ddnwx)

call dddmc (2.d0, 0, t1)
call ddnpwr (t1, ibz - ddnw2 * ddnbt, eps)
if (idb > 0) then
  call ddmdc (eps, dd1, nn1)
  write (ddldb, 4) ddnw2, nint (alog102*nn1)
4 format ('mpnw2, log10eps =',2i6)
endif
call dddmc (0.d0, 0, q1)

!   Perform Newton iterations with dynamic precision levels, using an
!   iteration formula similar to that used to evaluate reciprocals.

do j = 1, 10000
  if (idb > 0) write (ddldb, 5) j, kn
5 format ('j, kn =',3i6)

  call ddpolymul (kn, p2, q, r)
  call ddpolysub (kn, c1, r, s)
  call ddpolymul (kn, s, q, r)
  call ddpolyadd (kn, q, r, q)
  call ddsub (q(1:2,kn), q1, t1)

  if (idb > 0) then
    call ddmdc (t1, dd1, nn1)
    if (dd1 == 0.d0) then
      write (ddldb, 6)
6     format ('Newton error = 0')
    else
      write (ddldb, 7) nint (alog102*nn1)
7     format ('Newton error = 10^',i6)
    endif
  endif

  call ddabs (t1, t2)
  call ddcpr (t2, eps, ic1)
  if (ic1 < 0) then
    if (kn == n .and. ddnw2 == ddnw1) goto 100
    if (kn < n) then
      kn = min (2 * kn, n)
      call dddmc (0.d0, 0, q1)
    elseif (ddnw2 < ddnw1) then
      ddnw2 = min (2 * ddnw2, ddnw1)
      call dddmc (2.d0, 0, t1)
      call ddnpwr (t1, ibz - ddnw2 * ddnbt, eps)
      call dddmc (0.d0, 0, q1)
      if (idb > 0) then
        call ddmdc (eps, dd1, nn1)
        write (ddldb, 4) ddnw2, nint (alog102*nn1)
      endif
    endif
  else
    call ddeq (q(1:2,kn), q1)
  endif
enddo

write (ddldb, 8)
8 format ('*** DDBERNER: Loop end error')
call ddabrt

100 continue

if (idb > 0) write (ddldb, 9)
9 format ('Even zeta computation complete')

!   Multiply numerator polynomial by reciprocal of denominator polynomial.

call ddpolymul (n, p1, q, r)

!   Apply formula to produce Bernoulli numbers.

call dddmc (-2.d0, 0, t1)
call dddmc (1.d0, 0, t2)

do i = 1, n
  d1 = - dble (2*i-1) * dble (2*i)
  call ddmuld (t1, d1, t3)
  call ddeq (t3, t1)
  call ddmuld (cp2, 4.d0, t3)
  call ddmul (t3, t2, t4)
  call ddeq (t4, t2)
  call ddmuld (t1, 0.5d0, t3)
  call dddiv (t3, t2, t4)
  call ddabs (r(1:2,i), t3)
  call ddmul (t4, t3, berne(1:2,i))
enddo

if (idb > 0) write (ddldb, 10)
10 format ('Bernoulli number computation complete')
return
end subroutine ddberner

subroutine ddpolyadd (n, a, b, c)

!   This adds two polynomials, as is required by mpberne.
!   The output array C may be the same as A or B.

implicit none
integer, intent(in):: n
real (ddknd), intent(in):: a(1:2,0:n), b(1:2,0:n)
real (ddknd), intent(out):: c(1:2,0:n)
integer k
real (ddknd) t1(1:2), t2(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx


do k = 0, n
  call ddeq (a(1:2,k), t1)
  call ddeq (b(1:2,k), t2)
  call ddadd (t1, t2, c(1:2,k))
enddo

return
end subroutine ddpolyadd

subroutine ddpolysub (n, a, b, c)

!   This adds two polynomials, as is required by mpberne.
!   The output array C may be the same as A or B.

implicit none
integer, intent(in):: n
real (ddknd), intent(in):: a(1:2,0:n), b(1:2,0:n)
real (ddknd), intent(out):: c(1:2,0:n)
integer k
real (ddknd) t1(1:2), t2(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx


do k = 0, n
  call ddeq (a(1:2,k), t1)
  call ddeq (b(1:2,k), t2)
  call ddsub (t1, t2, c(1:2,k))
enddo

return
end subroutine ddpolysub

subroutine ddpolymul (n, a, b, c)

!   This adds two polynomials (ignoring high-order terms), as is required
!   by mpberne. The output array C may not be the same as A or B.

implicit none
integer, intent(in):: n
real (ddknd), intent(in):: a(1:2,0:n), b(1:2,0:n)
real (ddknd), intent(out):: c(1:2,0:n)
integer j, k
real (ddknd) t0(1:2), t1(1:2), t2(1:2), t3(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx


do k = 0, n
  call dddmc (0.d0, 0, t0)

  do j = 0, k
    call ddeq (a(1:2,j), t1)
    call ddeq (b(1:2,k-j), t2)
    call ddmul (t1, t2, t3)
    call ddadd (t0, t3, t2)
    call ddeq (t2, t0)
  enddo

  call ddeq (t0, c(1:2,k))
enddo

return
end subroutine ddpolymul

subroutine ddbesselinr (nu, rr, ss)

!   This evaluates the modified Bessel function BesselI (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.25.2 for modest RR,
!   and DLMF 10.40.1 for large RR, relative to precision.

implicit none
integer, intent(in):: nu
real (ddknd), intent(in):: rr(1:2)
real (ddknd), intent(out):: ss(1:2)
integer, parameter:: itrmax = 1000000
real (ddknd), parameter:: dfrac = 1.5d0, pi = 3.1415926535897932385d0
integer ic1, k, nua, n1
real (ddknd) d1
real (ddknd) f1(1:2), f2(1:2), sum(1:2), td(1:2), &
  tn(1:2), t1(1:2), t2(1:2), t3(1:2), t4(1:2), &
  rra(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

!   Check for RR = 0.

if (ddsgn (rr) == 0) then
  write (ddldb, 2)
2 format ('*** DDBESSELINR: Second argument is zero')
  call ddabrt
endif

ddnw1 = min (ddnw + 1, ddnwx)
nua = abs (nu)
call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw1*ddnbt, eps)
call ddabs (rr, rra)
call ddmdc (rra, d1, n1)
d1 = d1 * 2.d0 ** n1

if (d1 < dfrac * ddnw1 * dddpw) then
  call dddmc (1.d0, 0, tn)
  call dddmc (1.d0, 0, f1)
  call dddmc (1.d0, 0, f2)
  call ddmul (rra, rra, t2)
  call ddmuld (t2, 0.25d0, t1)

  do k = 1, nua
    call ddmuld (f2, dble (k), t2)
    call ddeq (t2, f2)
  enddo

  call ddmul (f1, f2, td)
  call dddiv (tn, td, t2)
  call ddeq (t2, sum)

  do k = 1, itrmax
    call ddmuld (f1, dble (k), t2)
    call ddeq (t2, f1)
    call ddmuld (f2, dble (k + nua), t2)
    call ddeq (t2, f2)
    call ddmul (t1, tn, t2)
    call ddeq (t2, tn)
    call ddmul (f1, f2, td)
    call dddiv (tn, td, t2)
    call ddadd (sum, t2, t3)
    call ddeq (t3, sum)

    call ddabs (t2, tc1)
    call ddmul (eps, sum, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 100
  enddo

  write (ddldb, 4)
  4 format ('*** DDBESSELINR: Loop end error 1')
  call ddabrt

100 continue

  call ddmuld (rra, 0.5d0, t1)
  call ddnpwr (t1, nua, t2)
  call ddmul (sum, t2, t3)
else
!  sum1 = mpreal (1.d0, mpnw)
!  t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
!  t2 = mpreal (1.d0, mpnw)
!  t3 = mpreal (1.d0, mpnw)

  call dddmc (1.d0, 0, sum)
  d1 = 4.d0 * dble (nua)**2
  call dddmc (d1, 0, t1)
  call dddmc (1.d0, 0, tn)
  call dddmc (1.d0, 0, td)

  do k = 1, itrmax
!  t2 = -t2 * (t1 - (2.d0*k - 1.d0)**2)
!  t3 = t3 * dble (k) * 8.d0 * xa
!  t4 = t2 / t3
!  sum1 = sum1 + t4

    d1 = 2.d0 * k - 1.d0
    call dddmc (d1, 0, t2)
    call ddmul (t2, t2, t3)
    call ddsub (t1, t3, t2)
    call ddmul (tn, t2, t3)
    call ddneg (t3, tn)
    call ddmuld (rra, 8.d0 * k, t2)
    call ddmul (td, t2, t3)
    call ddeq (t3, td)
    call dddiv (tn, td, t4)
    call ddadd (sum, t4, t3)
    call ddeq (t3, sum)

!   if (abs (t4) / abs (sum1) < eps) goto 110

    call ddabs (t4, tc1)
    call ddmul (eps, sum, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 110
  enddo

write (ddldb, 5)
5 format ('*** DDBESSELINR: Loop end error 2')
call ddabrt

110 continue

! t1 = exp (xa) / sqrt (2.d0 * mppi (mpnw) * xa)
! besseli = t1 * sum1

  call ddexp (rra, t1)
  call ddmuld (ddpicon, 2.d0, t2)
  call ddmul (t2, rra, t3)
  call ddsqrt (t3, t4)
  call dddiv (t1, t4, t2)
  call ddmul (t2, sum, t3)
endif

! if (x < 0.d0 .and. mod (nu, 2) /= 0) besseli = - besseli

if (ddsgn (rr) < 0 .and. mod (nu, 2) /= 0) then
  call ddneg (t3, t4)
  call ddeq (t4, t3)
endif

call ddeq (t3, ss)

return
end subroutine ddbesselinr

subroutine ddbesselir (qq, rr, ss)

!   This evaluates the modified Bessel function BesselI (QQ,RR) for QQ and RR
!   both MPR. The algorithm is DLMF formula 10.25.2 for modest RR, and
!   DLMF 10.40.1 for large RR, relative to precision.

implicit none
real (ddknd), intent(in):: qq(1:2), rr(1:2)
real (ddknd), intent(out):: ss(1:2)
integer ic1, i1, i2, k, n1
integer, parameter:: itrmax = 1000000
real (ddknd), parameter:: dfrac = 1.5d0, pi = 3.1415926535897932385d0
real (ddknd) d1
real (ddknd) f1(1:2), f2(1:2), sum(1:2), td(1:2), &
  tn(1:2), t1(1:2), t2(1:2), t3(1:2), t4(1:2), &
  rra(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

ddnw1 = min (ddnw + 1, ddnwx)
call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw1*ddnbt, eps)

!   If QQ is integer, call mpbesselinr; if qq < 0 and rr <= 0, then error.

call ddinfr (qq, t1, t2)
i1 = ddsgn (qq)
i2 = ddsgn (rr)
if (ddsgn (t2) == 0) then
  call ddmdc (qq, d1, n1)
  d1 = d1 * 2.d0**n1
  n1 = nint (d1)
  call ddbesselinr (n1, rr, t3)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (ddldb, 3)
3 format ('*** DDBESSELIR: First argument < 0 and second argument <= 0')
  call ddabrt
endif

call ddabs (rr, rra)
call ddmdc (rra, d1, n1)
d1 = d1 * 2.d0 ** n1

if (d1 < dfrac * ddnw1 * dddpw) then
  call dddmc (1.d0, 0, tn)
  call dddmc (1.d0, 0, f1)
  call ddadd (qq, f1, t1)
  call ddgammar (t1, f2)
  call ddmul (rra, rra, t2)
  call ddmuld (t2, 0.25d0, t1)

  call ddmul (f1, f2, td)
  call dddiv (tn, td, t2)
  call ddeq (t2, sum)

  do k = 1, itrmax
    call ddmuld (f1, dble (k), t2)
    call ddeq (t2, f1)
    call dddmc (dble (k), 0, t3)
    call ddadd (qq, t3, t4)
    call ddmul (f2, t4, t3)
    call ddeq (t3, f2)
    call ddmul (t1, tn, t2)
    call ddeq (t2, tn)
    call ddmul (f1, f2, td)
    call dddiv (tn, td, t2)
    call ddadd (sum, t2, t3)
    call ddeq (t3, sum)

    call ddabs (t2, tc1)
    call ddmul (eps, sum, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 100
  enddo

  write (ddldb, 4)
  4 format ('*** DDBESSELIR: Loop end error 1')
  call ddabrt

100 continue

  call ddmuld (rra, 0.5d0, t1)
  call ddpower (t1, qq, t2)
  call ddmul (sum, t2, t3)
else
!  sum1 = mpreal (1.d0, mpnw)
!  t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
!  t2 = mpreal (1.d0, mpnw)
!  t3 = mpreal (1.d0, mpnw)

  call dddmc (1.d0, 0, sum)
  call ddmul (qq, qq, t2)
  call ddmuld (t2, 4.d0, t1)
  call dddmc (1.d0, 0, tn)
  call dddmc (1.d0, 0, td)

  do k = 1, itrmax
!  t2 = -t2 * (t1 - (2.d0*k - 1.d0)**2)
!  t3 = t3 * dble (k) * 8.d0 * xa
!  t4 = t2 / t3
!  sum1 = sum1 + t4

    d1 = 2.d0 * k - 1.d0
    call dddmc (d1, 0, t2)
    call ddmul (t2, t2, t3)
    call ddsub (t1, t3, t2)
    call ddmul (tn, t2, t3)
    call ddneg (t3, tn)
    call ddmuld (rra, 8.d0 * k, t2)
    call ddmul (td, t2, t3)
    call ddeq (t3, td)
    call dddiv (tn, td, t4)
    call ddadd (sum, t4, t3)
    call ddeq (t3, sum)

!   if (abs (t4) / abs (sum1) < eps) goto 110

    call ddabs (t4, tc1)
    call ddmul (eps, sum, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 110
  enddo

write (ddldb, 5)
5 format ('*** DDBESSELIR: Loop end error 2')
call ddabrt

110 continue

! t1 = exp (xa) / sqrt (2.d0 * mppi (mpnw) * xa)
! besseli = t1 * sum1

  call ddexp (rra, t1)
  call ddmuld (ddpicon, 2.d0, t2)
  call ddmul (t2, rra, t3)
  call ddsqrt (t3, t4)
  call dddiv (t1, t4, t2)
  call ddmul (t2, sum, t3)
endif

120 continue

call ddeq (t3, ss)

return
end subroutine ddbesselir

subroutine ddbesseljnr (nu, rr, ss)

!   This evaluates the modified Bessel function BesselJ (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.2.2 for modest RR,
!   and DLMF 10.17.3 for large RR, relative to precision.

implicit none
integer, intent(in):: nu
real (ddknd), intent(in):: rr(1:2)
real (ddknd), intent(out):: ss(1:2)
integer, parameter:: itrmax = 1000000
real (ddknd), parameter:: dfrac = 1.5d0, pi = 3.1415926535897932385d0
integer ic1, ic2, k, nua, n1
real (ddknd) d1, d2
real (ddknd) f1(1:2), f2(1:2), sum1(1:2), &
  sum2(1:2), td1(1:2), td2(1:2), tn1(1:2), &
  tn2(1:2), t1(1:2), t2(1:2), t3(1:2), &
  t41(1:2), t42(1:2), t5(1:2), rra(1:2), &
  rr2(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

!   Check for RR = 0.

if (ddsgn (rr) == 0) then
  write (ddldb, 2)
2 format ('*** DDBESSELJNR: Second argument is zero')
  call ddabrt
endif

ddnw1 = min (2 * ddnw + 1, ddnwx)
nua = abs (nu)
call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw*ddnbt, eps)
call ddmdc (rr, d1, n1)
d1 = abs (d1) * 2.d0 ** n1

if (d1 < dfrac * ddnw1 * dddpw) then
  ddnw1 = min (ddnw + nint (d1 / (dfrac * dddpw)), 2 * ddnw + 1, ddnwx)
  call ddabs (rr, rra)
  call dddmc (1.d0, 0, tn1)
  call dddmc (1.d0, 0, f1)
  call dddmc (1.d0, 0, f2)
  call ddmul (rra, rra, t2)
  call ddmuld (t2, 0.25d0, t1)

  do k = 1, nua
    call ddmuld (f2, dble (k), t2)
    call ddeq (t2, f2)
  enddo

  call ddmul (f1, f2, td1)
  call dddiv (tn1, td1, t2)
  call ddeq (t2, sum1)

  do k = 1, itrmax
    call ddmuld (f1, dble (k), t2)
    call ddeq (t2, f1)
    call ddmuld (f2, dble (k + nua), t2)
    call ddeq (t2, f2)
    call ddmul (t1, tn1, t2)
    call ddneg (t2, tn1)
    call ddmul (f1, f2, td1)
    call dddiv (tn1, td1, t2)
    call ddadd (sum1, t2, t3)
    call ddeq (t3, sum1)

    call ddabs (t2, tc1)
    call ddmul (eps, sum1, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 100
  enddo

  write (ddldb, 4)
4 format ('*** DDBESSELJNR: Loop end error 1')
  call ddabrt

100 continue

  call ddmuld (rra, 0.5d0, t1)
  call ddnpwr (t1, nua, t2)
  call ddmul (sum1, t2, t3)
else
! xa2 = xa**2
! t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
! tn1 = mpreal (1.d0, mpnw)
! tn2 = (t1 - 1.d0) / 8.d0
! td1 = mpreal (1.d0, mpnw)
! td2 = xa
! sum1 = tn1 / td1
! sum2 = tn2 / t32

  ddnw1 = min (ddnw + 1, ddnwx)
  call ddabs (rr, rra)
  call ddmul (rra, rra, rr2)
  d1 = 4.d0 * dble (nua)**2
  call dddmc (d1, 0, t1)
  call dddmc (1.d0, 0, tn1)
  call ddsub (t1, tn1, t2)
  call dddivd (t2, 8.d0, tn2)
  call dddmc (1.d0, 0, td1)
  call ddeq (rra, td2)
  call dddiv (tn1, td1, sum1)
  call dddiv (tn2, td2, sum2)

  do k = 1, itrmax
!   tn1 = -tn1 * (t1 - (2.d0*(2.d0*k-1.d0) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k) - 1.d0)**2)

    d1 = (2.d0*(2.d0*k-1.d0) - 1.d0)**2
    d2 = (2.d0*(2.d0*k) - 1.d0)**2
    call dddmc (d1, 0, t2)
    call ddsub (t1, t2, t3)
    call dddmc (d2, 0, t2)
    call ddsub (t1, t2, t5)
    call ddmul (t3, t5, t2)
    call ddmul (tn1, t2, t3)
    call ddneg (t3, tn1)

!   td1 = td1 * dble (2*k-1) * dble (2*k) * 64.d0 * xa2

    d1 = dble (2*k-1) * dble (2*k) * 64.d0
    call ddmuld (td1, d1, t2)
    call ddmul (t2, rr2, td1)

!   t41 = tn1 / td1
!   sum1 = sum1 + t41

    call dddiv (tn1, td1, t41)
    call ddadd (sum1, t41, t2)
    call ddeq (t2, sum1)

!   tn2 = -tn2 * (t1 - (2.d0*(2.d0*k) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k+1.d0) - 1.d0)**2)

    d1 = (2.d0*(2.d0*k) - 1.d0)**2
    d2 = (2.d0*(2.d0*k+1.d0) - 1.d0)**2
    call dddmc (d1, 0, t2)
    call ddsub (t1, t2, t3)
    call dddmc (d2, 0, t2)
    call ddsub (t1, t2, t5)
    call ddmul (t3, t5, t2)
    call ddmul (tn2, t2, t3)
    call ddneg (t3, tn2)

!   td2 = td2 * dble (2*k) * dble (2*k+1) * 64.d0 * xa2

    d1 = dble (2*k) * dble (2*k+1) * 64.d0
    call ddmuld (td2, d1, t2)
    call ddmul (t2, rr2, td2)

!  t42 = tn2 / td2
!  sum2 = sum2 + t42

    call dddiv (tn2, td2, t42)
    call ddadd (sum2, t42, t2)
    call ddeq (t2, sum2)

!  if (abs (t41) / abs (sum1) < eps .and. abs (t42) / abs (sum2) < eps ) goto 110

    call ddabs (t41, tc1)
    call ddmul (eps, sum1, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    call ddabs (t42, tc1)
    call ddmul (eps, sum2, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic2)
    if (ic1 <= 0 .and. ic2 <= 0) goto 110
  enddo

  write (ddldb, 5)
5 format ('*** DDBESSELJNR: Loop end error 2')
  call ddabrt

110 continue

! t1 = xa - 0.5d0 * nua * pi - 0.25d0 * pi
! besselj = sqrt (2.d0 / (pi * xa)) * (cos (t1) * sum1 - sin (t1) * sum2)

  call ddmuld (ddpicon, 0.5d0 * nua, t1)
  call ddsub (rra, t1, t2)
  call ddmuld (ddpicon, 0.25d0, t1)
  call ddsub (t2, t1, t3)
  call ddcssnr (t3, t41, t42)
  call ddmul (t41, sum1, t1)
  call ddmul (t42, sum2, t2)
  call ddsub (t1, t2, t5)
  call ddmul (ddpicon, rra, t1)
  call dddmc (2.d0, 0, t2)
  call dddiv (t2, t1, t3)
  call ddsqrt (t3, t1)
  call ddmul (t1, t5, t3)
endif

if (mod (nu, 2) /= 0) then
!  if (nu < 0 .and. x > 0.d0 .or. nu > 0 .and. x < 0.d0) besselj = - besselj

  ic1 = ddsgn (rr)
  if (nu < 0 .and. ic1 > 0 .or. nu > 0 .and. ic1 < 0) then
    call ddneg (t3, t2)
    call ddeq (t2, t3)
  endif
endif

call ddeq (t3, ss)

return
end subroutine ddbesseljnr

subroutine ddbesseljr (qq, rr, ss)

!   This evaluates the modified Bessel function BesselJ (QQ,RR) for QQ and RR
!   both MPR. The algorithm is DLMF formula 10.2.2 for modest RR,
!   and DLMF 10.17.3 for large RR, relative to precision.

implicit none
real (ddknd), intent(in):: qq(1:2), rr(1:2)
real (ddknd), intent(out):: ss(1:2)
integer ic1, ic2, i1, i2, k, n1
real (ddknd) d1, d2
integer, parameter:: itrmax = 1000000
real (ddknd), parameter:: dfrac = 1.5d0, pi = 3.1415926535897932385d0
real (ddknd) f1(1:2), f2(1:2), sum1(1:2), &
  sum2(1:2), td1(1:2), td2(1:2), tn1(1:2), &
  tn2(1:2), t1(1:2), t2(1:2), t3(1:2), &
  t4(1:2), t41(1:2), t42(1:2), t5(1:2), &
  rra(1:2), rr2(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

ddnw1 = min (2 * ddnw + 1, ddnwx)
call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw*ddnbt, eps)

!   If QQ is integer, call mpbesseljnr; if RR <= 0, then error.

call ddinfr (qq, t1, t2)
i1 = ddsgn (qq)
i2 = ddsgn (rr)
if (ddsgn (t2) == 0) then
  call ddmdc (qq, d1, n1)
  d1 = d1 * 2.d0**n1
  n1 = nint (d1)
  call ddbesseljnr (n1, rr, t3)
  goto 120
elseif (i2 <= 0) then
  write (ddldb, 3)
3 format ('*** DDBESSELJR: Second argument <= 0')
  call ddabrt
endif

call ddmdc (rr, d1, n1)
d1 = abs (d1) * 2.d0 ** n1

if (d1 < dfrac * ddnw1 * dddpw) then
  ddnw1 = min (ddnw + nint (d1 / (dfrac * dddpw)), 2 * ddnw + 1, ddnwx)
  call ddabs (rr, rra)
  call dddmc (1.d0, 0, tn1)
  call dddmc (1.d0, 0, f1)
  call ddadd (qq, f1, t1)
  call ddgammar (t1, f2)
  call ddmul (rra, rra, t2)
  call ddmuld (t2, 0.25d0, t1)

  call ddmul (f1, f2, td1)
  call dddiv (tn1, td1, t2)
  call ddeq (t2, sum1)

  do k = 1, itrmax
    call ddmuld (f1, dble (k), t2)
    call ddeq (t2, f1)
    call dddmc (dble (k), 0, t3)
    call ddadd (qq, t3, t4)
    call ddmul (f2, t4, t3)
    call ddeq (t3, f2)
    call ddmul (t1, tn1, t2)
    call ddneg (t2, tn1)
    call ddmul (f1, f2, td1)
    call dddiv (tn1, td1, t2)
    call ddadd (sum1, t2, t3)
    call ddeq (t3, sum1)

    call ddabs (t2, tc1)
    call ddmul (eps, sum1, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 100
  enddo

  write (ddldb, 4)
4 format ('*** DDBESSELJR: Loop end error 1')
  call ddabrt

100 continue

  call ddmuld (rr, 0.5d0, t1)
  call ddpower (t1, qq, t2)
  call ddmul (sum1, t2, t3)
else
! xa2 = xa**2
! t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
! tn1 = mpreal (1.d0, mpnw)
! tn2 = (t1 - 1.d0) / 8.d0
! td1 = mpreal (1.d0, mpnw)
! td2 = xa
! sum1 = tn1 / td1
! sum2 = tn2 / t32

  ddnw1 = min (ddnw + 1, ddnwx)
  call ddabs (rr, rra)
  call ddmul (rra, rra, rr2)
  call ddmul (qq, qq, t2)
  call ddmuld (t2, 4.d0, t1)
  call dddmc (1.d0, 0, tn1)
  call ddsub (t1, tn1, t2)
  call dddivd (t2, 8.d0, tn2)
  call dddmc (1.d0, 0, td1)
  call ddeq (rra, td2)
  call dddiv (tn1, td1, sum1)
  call dddiv (tn2, td2, sum2)

  do k = 1, itrmax
!   tn1 = -tn1 * (t1 - (2.d0*(2.d0*k-1.d0) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k) - 1.d0)**2)

    d1 = (2.d0*(2.d0*k-1.d0) - 1.d0)**2
    d2 = (2.d0*(2.d0*k) - 1.d0)**2
    call dddmc (d1, 0, t2)
    call ddsub (t1, t2, t3)
    call dddmc (d2, 0, t2)
    call ddsub (t1, t2, t5)
    call ddmul (t3, t5, t2)
    call ddmul (tn1, t2, t3)
    call ddneg (t3, tn1)

!   td1 = td1 * dble (2*k-1) * dble (2*k) * 64.d0 * xa2

    d1 = dble (2*k-1) * dble (2*k) * 64.d0
    call ddmuld (td1, d1, t2)
    call ddmul (t2, rr2, td1)

!   t41 = tn1 / td1
!   sum1 = sum1 + t41

    call dddiv (tn1, td1, t41)
    call ddadd (sum1, t41, t2)
    call ddeq (t2, sum1)

!   tn2 = -tn2 * (t1 - (2.d0*(2.d0*k) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k+1.d0) - 1.d0)**2)

    d1 = (2.d0*(2.d0*k) - 1.d0)**2
    d2 = (2.d0*(2.d0*k+1.d0) - 1.d0)**2
    call dddmc (d1, 0, t2)
    call ddsub (t1, t2, t3)
    call dddmc (d2, 0, t2)
    call ddsub (t1, t2, t5)
    call ddmul (t3, t5, t2)
    call ddmul (tn2, t2, t3)
    call ddneg (t3, tn2)

!   td2 = td2 * dble (2*k) * dble (2*k+1) * 64.d0 * xa2

    d1 = dble (2*k) * dble (2*k+1) * 64.d0
    call ddmuld (td2, d1, t2)
    call ddmul (t2, rr2, td2)

!  t42 = tn2 / td2
!  sum2 = sum2 + t42

    call dddiv (tn2, td2, t42)
    call ddadd (sum2, t42, t2)
    call ddeq (t2, sum2)

!  if (abs (t41) / abs (sum1) < eps .and. abs (t42) / abs (sum2) < eps ) goto 110

    call ddabs (t41, tc1)
    call ddmul (eps, sum1, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    call ddabs (t42, tc1)
    call ddmul (eps, sum2, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic2)
    if (ic1 <= 0 .and. ic2 <= 0) goto 110
  enddo

  write (ddldb, 5)
5 format ('*** DDBESSELJR: Loop end error 2')
  call ddabrt

110 continue

! t1 = xa - 0.5d0 * nua * pi - 0.25d0 * pi
! besselj = sqrt (2.d0 / (pi * xa)) * (cos (t1) * sum1 - sin (t1) * sum2)

!   call mpmuld (mppicon, 0.5d0 * nua, t1, mpnw1)

  call ddmul (ddpicon, qq, t2)
  call ddmuld (t2, 0.5d0, t1)
  call ddsub (rra, t1, t2)
  call ddmuld (ddpicon, 0.25d0, t1)
  call ddsub (t2, t1, t3)
  call ddcssnr (t3, t41, t42)
  call ddmul (t41, sum1, t1)
  call ddmul (t42, sum2, t2)
  call ddsub (t1, t2, t5)
  call ddmul (ddpicon, rra, t1)
  call dddmc (2.d0, 0, t2)
  call dddiv (t2, t1, t3)
  call ddsqrt (t3, t1)
  call ddmul (t1, t5, t3)
endif

! if (mod (nu, 2) /= 0) then
!  if (nu < 0 .and. x > 0.d0 .or. nu > 0 .and. x < 0.d0) besselj = - besselj

!   ic1 = mpsgn (rr)
!   if (nu < 0 .and. ic1 > 0 .or. nu > 0 .and. ic1 < 0) then
!     call mpneg (t3, t2, mpnw1)
!     call mpeq (t2, t3, mpnw1)
!   endif
! endif

120 continue

call ddeq (t3, ss)

return
end subroutine ddbesseljr

subroutine ddbesselknr (nu, rr, ss)

!   This evaluates the modified Bessel function BesselK (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.31.1 for modest RR,
!   and DLMF 10.40.2 for large RR, relative to precision.

implicit none
integer, intent(in):: nu
real (ddknd), intent(in):: rr(1:2)
real (ddknd), intent(out):: ss(1:2)
integer, parameter:: itrmax = 1000000
integer ic1, k, nua, n1
real (ddknd) d1
real (ddknd), parameter:: dfrac = 1.5d0, egam = 0.5772156649015328606d0, &
  pi = 3.1415926535897932385d0
real (ddknd) f1(1:2), f2(1:2), f3(1:2), f4(1:2), &
  f5(1:2), sum1(1:2), sum2(1:2), sum3(1:2), td(1:2), &
  tn(1:2), t1(1:2), t2(1:2), t3(1:2), t4(1:2), &
  rra(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

!   Check for RR = 0.

if (ddsgn (rr) == 0) then
  write (ddldb, 2)
2 format ('*** DDBESSELKNR: Second argument is zero')
  call ddabrt
endif

ddnw1 = min (ddnw + 1, ddnwx)
nua = abs (nu)
ddnw1 = min (ddnw + 1, ddnwx)
call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw1*ddnbt, eps)
call ddabs (rr, rra)
call ddmdc (rra, d1, n1)
d1 = d1 * 2.d0 ** n1

if (d1 < dfrac * ddnw1 * dddpw) then
  call ddmul (rra, rra, t2)
  call ddmuld (t2, 0.25d0, t1)
  call dddmc (1.d0, 0, f1)
  call dddmc (1.d0, 0, f2)
  call dddmc (1.d0, 0, f3)
  call dddmc (0.d0, 0,  sum1)

  do k = 1, nua - 1
    call ddmuld (f1, dble (k), t2)
    call ddeq (t2, f1)
  enddo

  do k = 0, nua - 1
    if (k > 0) then
      call dddivd (f1, dble (nua - k), t2)
      call ddeq (t2, f1)
      call ddmul (t1, f2, t2)
      call ddneg (t2, f2)
      call ddmuld (f3, dble (k), t2)
      call ddeq (t2, f3)
    endif
    call ddmul (f1, f2, t3)
    call dddiv (t3, f3, t2)
    call ddadd (sum1, t2, t3)
    call ddeq (t3, sum1)
  enddo

  call ddmuld (sum1, 0.5d0, t2)
  call ddmuld (rra, 0.5d0, t3)
  call ddnpwr (t3, nua, t4)
  call dddiv (t2, t4, sum1)

  call ddmuld (rra, 0.5d0, t2)
  call ddlog (t2, t3)
  d1 = (-1.d0) ** (nua + 1)
  call ddmuld (t3, d1, t2)
  call ddbesselinr (nua, rra, t3)
  call ddmul (t2, t3, sum2)

  call ddneg (ddegammacon, f1)
  call ddeq (f1, f2)
  call dddmc (1.d0, 0, f3)
  call dddmc (1.d0, 0, f4)
  call dddmc (1.d0, 0, f5)

  do k = 1, nua
    call dddmc (1.d0, 0, t2)
    call dddivd (t2, dble (k), t3)
    call ddadd (f2, t3, t4)
    call ddeq (t4, f2)
    call ddmuld (f5, dble (k), t2)
    call ddeq (t2, f5)
  enddo

  call ddadd (f1, f2, t2)
  call ddmul (t2, f3, t3)
  call ddmul (f4, f5, t4)
  call dddiv (t3, t4, sum3)

  do k = 1, itrmax
    call dddmc (1.d0, 0, t2)
    call dddivd (t2, dble (k), t3)
    call ddadd (f1, t3, t4)
    call ddeq (t4, f1)
    call dddivd (t2, dble (nua + k), t3)
    call ddadd (f2, t3, t4)
    call ddeq (t4, f2)
    call ddmul (t1, f3, t2)
    call ddeq (t2, f3)
    call ddmuld (f4, dble (k), t2)
    call ddeq (t2, f4)
    call ddmuld (f5, dble (nua + k), t2)
    call ddeq (t2, f5)
    call ddadd (f1, f2, t2)
    call ddmul (t2, f3, t3)
    call ddmul (f4, f5, t4)
    call dddiv (t3, t4, t2)
    call ddadd (sum3, t2, t3)
    call ddeq (t3, sum3)

    call ddabs (t2, tc1)
    call ddmul (eps, sum3, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 100
  enddo

  write (ddldb, 5)
5 format ('*** DDBESSELKNR: Loop end error 1')
  call ddabrt

100 continue

  call ddmuld (rra, 0.5d0, t2)
  call ddnpwr (t2, nua, t3)
  d1 = (-1.d0)**nua * 0.5d0
  call ddmuld (t3, d1, t4)
  call ddmul (t4, sum3, t2)
  call ddeq (t2, sum3)
  call ddadd (sum1, sum2, t2)
  call ddadd (t2, sum3, t3)
else
!  sum1 = mpreal (1.d0, mpnw)
!  t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
!  t2 = mpreal (1.d0, mpnw)
!  t3 = mpreal (1.d0, mpnw)

  call dddmc (1.d0, 0, sum1)
  d1 = 4.d0 * dble (nua)**2
  call dddmc (d1, 0, t1)
  call dddmc (1.d0, 0, tn)
  call dddmc (1.d0, 0, td)

  do k = 1, itrmax
!  t2 = t2 * (t1 - (2.d0*k - 1.d0)**2)
!  t3 = t3 * dble (k) * 8.d0 * xa
!  t4 = t2 / t3
!  sum1 = sum1 + t4

    d1 = 2.d0 * k - 1.d0
    call dddmc (d1, 0, t2)
    call ddmul (t2, t2, t3)
    call ddsub (t1, t3, t2)
    call ddmul (tn, t2, t3)
    call ddeq (t3, tn)
    call ddmuld (rra, 8.d0 * k, t2)
    call ddmul (td, t2, t3)
    call ddeq (t3, td)
    call dddiv (tn, td, t4)
    call ddadd (sum1, t4, t3)
    call ddeq (t3, sum1)

!   if (abs (t4) / abs (sum1) < eps) goto 110

    call ddabs (t4, tc1)
    call ddmul (eps, sum1, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 110
  enddo

write (ddldb, 6)
6 format ('*** DDBESSELKNR: Loop end error 2')
call ddabrt

110 continue

! t1 = sqrt (mppi (mpnw) / (2.d0 * xa)) / exp (xa)
! besseli = t1 * sum1

  call ddexp (rra, t1)
  call ddmuld (rra, 2.d0, t2)
  call dddiv (ddpicon, t2, t3)
  call ddsqrt (t3, t4)
  call dddiv (t4, t1, t2)
  call ddmul (t2, sum1, t3)
endif

! if (x < 0.d0 .and. mod (nu, 2) /= 0) besselk = - besselk

if (ddsgn (rr) < 0 .and. mod (nu, 2) /= 0) then
  call ddneg (t3, t4)
  call ddeq (t4, t3)
endif
call ddeq (t3, ss)
return
end subroutine ddbesselknr

subroutine ddbesselkr (qq, rr, ss)

!   This evaluates the Bessel function BesselK (QQ,RR) for QQ and RR
!   both MPR. This uses DLMF formula 10.27.4.

implicit none
real (ddknd), intent(in):: qq(1:2), rr(1:2)
real (ddknd), intent(out):: ss(1:2)
integer i1, i2, n1
real (ddknd) d1
real (ddknd) t1(1:2), t2(1:2), t3(1:2), t4(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

ddnw1 = min (ddnw + 1, ddnwx)

!   If QQ is integer, call mpbesselknr; if qq < 0 and rr <= 0, then error.

call ddinfr (qq, t1, t2)
i1 = ddsgn (qq)
i2 = ddsgn (rr)
if (ddsgn (t2) == 0) then
  call ddmdc (qq, d1, n1)
  d1 = d1 * 2.d0**n1
  n1 = nint (d1)
  call ddbesselknr (n1, rr, t1)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (ddldb, 2)
2 format ('*** DDBESSELKR: First argument < 0 and second argument <= 0')
  call ddabrt
endif

call ddneg (qq, t1)
call ddbesselir (t1, rr, t2)
call ddbesselir (qq, rr, t3)
call ddsub (t2, t3, t4)
call ddmul (qq, ddpicon, t1)
call ddcssnr (t1, t2, t3)
call dddiv (t4, t3, t2)
call ddmul (ddpicon, t2, t3)
call ddmuld (t3, 0.5d0, t1)

120 continue

call ddeq (t1, ss)
return
end subroutine ddbesselkr

subroutine ddbesselynr (nu, rr, ss)

!   This evaluates the modified Bessel function BesselY (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.8.1 for modest RR,
!   and DLMF 10.17.4 for large RR, relative to precision.

implicit none
integer, intent(in):: nu
real (ddknd), intent(in):: rr(1:2)
real (ddknd), intent(out):: ss(1:2)
integer, parameter:: itrmax = 1000000
real (ddknd), parameter:: dfrac = 1.5d0, egam = 0.5772156649015328606d0, &
  pi = 3.1415926535897932385d0
integer ic1, ic2, k, nua, n1
real (ddknd) d1, d2
real (ddknd) f1(1:2), f2(1:2), f3(1:2), f4(1:2), &
  f5(1:2), rra(1:2), rr2(1:2), sum1(1:2), &
  sum2(1:2), sum3(1:2), td1(1:2), td2(1:2), &
  tn1(1:2), tn2(1:2), t1(1:2), t2(1:2), &
  t3(1:2), t4(1:2), t41(1:2), t42(1:2), &
  t5(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

!   Check for RR = 0.

if (ddsgn (rr) == 0) then
  write (ddldb, 2)
2 format ('*** DDBESSELYNR: argument is negative or too large')
  call ddabrt
endif

ddnw1 = min (2 * ddnw + 1, ddnwx)
nua = abs (nu)
call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw*ddnbt, eps)
call ddmdc (rr, d1, n1)
d1 = abs (d1) * 2.d0 ** n1

if (d1 < dfrac * ddnw1 * dddpw) then
  ddnw1 = min (ddnw + nint (d1 / (dfrac * dddpw)), 2 * ddnw + 1, ddnwx)
  call ddabs (rr, rra)
  call ddmul (rra, rra, t2)
  call ddmuld (t2, 0.25d0, t1)
  call dddmc (1.d0, 0, f1)
  call dddmc (1.d0, 0, f2)
  call dddmc (1.d0, 0, f3)
  call dddmc (0.d0, 0,  sum1)

  do k = 1, nua - 1
    call ddmuld (f1, dble (k), t2)
    call ddeq (t2, f1)
  enddo

  do k = 0, nua - 1
    if (k > 0) then
      call dddivd (f1, dble (nua - k), t2)
      call ddeq (t2, f1)
      call ddmul (t1, f2, t2)
      call ddeq (t2, f2)
      call ddmuld (f3, dble (k), t2)
      call ddeq (t2, f3)
    endif
    call ddmul (f1, f2, t3)
    call dddiv (t3, f3, t2)
    call ddadd (sum1, t2, t3)
    call ddeq (t3, sum1)
  enddo

  call ddmuld (rra, 0.5d0, t3)
  call ddnpwr (t3, nua, t4)
  call dddiv (sum1, t4, t3)
  call ddneg (t3, sum1)

  call ddmuld (rra, 0.5d0, t2)
  call ddlog (t2, t3)
  call ddmuld (t3, 2.d0, t2)
  call ddbesseljnr (nua, rra, t3)
  call ddmul (t2, t3, sum2)

  call ddneg (ddegammacon, f1)
  call ddeq (f1, f2)
  call dddmc (1.d0, 0, f3)
  call dddmc (1.d0, 0, f4)
  call dddmc (1.d0, 0, f5)

  do k = 1, nua
    call dddmc (1.d0, 0, t2)
    call dddivd (t2, dble (k), t3)
    call ddadd (f2, t3, t4)
    call ddeq (t4, f2)
    call ddmuld (f5, dble (k), t2)
    call ddeq (t2, f5)
  enddo

  call ddadd (f1, f2, t2)
  call ddmul (t2, f3, t3)
  call ddmul (f4, f5, t4)
  call dddiv (t3, t4, sum3)

  do k = 1, itrmax
    call dddmc (1.d0, 0, t2)
    call dddivd (t2, dble (k), t3)
    call ddadd (f1, t3, t4)
    call ddeq (t4, f1)
    call dddivd (t2, dble (nua + k), t3)
    call ddadd (f2, t3, t4)
    call ddeq (t4, f2)
    call ddmul (t1, f3, t2)
    call ddneg (t2, f3)
    call ddmuld (f4, dble (k), t2)
    call ddeq (t2, f4)
    call ddmuld (f5, dble (nua + k), t2)
    call ddeq (t2, f5)
    call ddadd (f1, f2, t2)
    call ddmul (t2, f3, t3)
    call ddmul (f4, f5, t4)
    call dddiv (t3, t4, t2)
    call ddadd (sum3, t2, t3)
    call ddeq (t3, sum3)

    call ddabs (t2, tc1)
    call ddmul (eps, sum3, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 100
  enddo

  write (ddldb, 6)
6 format ('*** DDBESSELYNR: Loop end error 1')
  call ddabrt

100 continue

  call ddmuld (rra, 0.5d0, t2)
  call ddnpwr (t2, nua, t3)
  call ddmul (t3, sum3, t2)
  call ddneg (t2, sum3)

  call ddadd (sum1, sum2, t2)
  call ddadd (t2, sum3, t4)
  call ddeq (ddpicon, t2)
  call dddiv (t4, t2, t3)
else

! xa2 = xa**2
! t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
! tn1 = mpreal (1.d0, mpnw)
! tn2 = (t1 - 1.d0) / 8.d0
! td1 = mpreal (1.d0, mpnw)
! td2 = xa
! sum1 = tn1 / td1
! sum2 = tn2 / t32

  ddnw1 = min (ddnw + 1, ddnwx)
  call ddabs (rr, rra)
  call ddmul (rra, rra, rr2)
  d1 = 4.d0 * dble (nua)**2
  call dddmc (d1, 0, t1)
  call dddmc (1.d0, 0, tn1)
  call ddsub (t1, tn1, t2)
  call dddivd (t2, 8.d0, tn2)
  call dddmc (1.d0, 0, td1)
  call ddeq (rra, td2)
  call dddiv (tn1, td1, sum1)
  call dddiv (tn2, td2, sum2)

  do k = 1, itrmax
!   tn1 = -tn1 * (t1 - (2.d0*(2.d0*k-1.d0) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k) - 1.d0)**2)

    d1 = (2.d0*(2.d0*k-1.d0) - 1.d0)**2
    d2 = (2.d0*(2.d0*k) - 1.d0)**2
    call dddmc (d1, 0, t2)
    call ddsub (t1, t2, t3)
    call dddmc (d2, 0, t2)
    call ddsub (t1, t2, t5)
    call ddmul (t3, t5, t2)
    call ddmul (tn1, t2, t3)
    call ddneg (t3, tn1)

!   td1 = td1 * dble (2*k-1) * dble (2*k) * 64.d0 * xa2

    d1 = dble (2*k-1) * dble (2*k) * 64.d0
    call ddmuld (td1, d1, t2)
    call ddmul (t2, rr2, td1)

!   t41 = tn1 / td1
!   sum1 = sum1 + t41

    call dddiv (tn1, td1, t41)
    call ddadd (sum1, t41, t2)
    call ddeq (t2, sum1)

!   tn2 = -tn2 * (t1 - (2.d0*(2.d0*k) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k+1.d0) - 1.d0)**2)

    d1 = (2.d0*(2.d0*k) - 1.d0)**2
    d2 = (2.d0*(2.d0*k+1.d0) - 1.d0)**2
    call dddmc (d1, 0, t2)
    call ddsub (t1, t2, t3)
    call dddmc (d2, 0, t2)
    call ddsub (t1, t2, t5)
    call ddmul (t3, t5, t2)
    call ddmul (tn2, t2, t3)
    call ddneg (t3, tn2)

!   td2 = td2 * dble (2*k) * dble (2*k+1) * 64.d0 * xa2

    d1 = dble (2*k) * dble (2*k+1) * 64.d0
    call ddmuld (td2, d1, t2)
    call ddmul (t2, rr2, td2)

!  t42 = tn2 / td2
!  sum2 = sum2 + t42

    call dddiv (tn2, td2, t42)
    call ddadd (sum2, t42, t2)
    call ddeq (t2, sum2)

!  if (abs (t41) / abs (sum1) < eps .and. abs (t42) / abs (sum2) < eps ) goto 110

    call ddabs (t41, tc1)
    call ddmul (eps, sum1, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    call ddabs (t42, tc1)
    call ddmul (eps, sum2, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic2)
    if (ic1 <= 0 .and. ic2 <= 0) goto 110
  enddo

  write (ddldb, 5)
5 format ('*** DDBESSELYNR: Loop end error 2')
  call ddabrt

110 continue

! t1 = xa - 0.5d0 * nua * pi - 0.25d0 * pi
! besselj = sqrt (2.d0 / (pi * xa)) * (cos (t1) * sum1 - sin (t1) * sum2)

  call ddmuld (ddpicon, 0.5d0 * nua, t1)
  call ddsub (rra, t1, t2)
  call ddmuld (ddpicon, 0.25d0, t1)
  call ddsub (t2, t1, t3)
  call ddcssnr (t3, t41, t42)
  call ddmul (t42, sum1, t1)
  call ddmul (t41, sum2, t2)
  call ddadd (t1, t2, t5)
  call ddmul (ddpicon, rra, t1)
  call dddmc (2.d0, 0, t2)
  call dddiv (t2, t1, t3)
  call ddsqrt (t3, t1)
  call ddmul (t1, t5, t3)
endif

if (mod (nu, 2) /= 0) then
!   if (nu < 0 .and. x > 0.d0 .or. nu > 0 .and. x < 0.d0) bessely = - bessely

  ic1 = ddsgn (rr)
  if (nu < 0 .and. ic1 > 0 .or. nu > 0 .and. ic1 < 0) then
    call ddneg (t3, t4)
    call ddeq (t4, t3)
  endif
endif

call ddeq (t3, ss)
return
end subroutine ddbesselynr

subroutine ddbesselyr (qq, rr, ss)

!   This evaluates the modified Bessel function BesselY (QQ,RR).
!   NU is an integer. The algorithm is DLMF formula 10.2.2.

implicit none
real (ddknd), intent(in):: qq(1:2), rr(1:2)
real (ddknd), intent(out):: ss(1:2)
integer i1, i2, n1
real (ddknd) d1
real (ddknd) t1(1:2), t2(1:2), t3(1:2), t4(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

ddnw1 = min (ddnw + 1, ddnwx)

!   If QQ is integer, call mpbesselynr; if qq < 0 and rr <= 0, then error.

call ddinfr (qq, t1, t2)
i1 = ddsgn (qq)
i2 = ddsgn (rr)
if (ddsgn (t2) == 0) then
  call ddmdc (qq, d1, n1)
  d1 = d1 * 2.d0**n1
  n1 = nint (d1)
  call ddbesselynr (n1, rr, t1)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (ddldb, 2)
2 format ('*** DDBESSELYR: First argument < 0 and second argument <= 0')
  call ddabrt
endif

call ddmul (qq, ddpicon, t1)
call ddcssnr (t1, t2, t3)
call ddbesseljr (qq, rr, t4)
call ddmul (t4, t2, t1)
call ddneg (qq, t2)
call ddbesseljr (t2, rr, t4)
call ddsub (t1, t4, t2)
call dddiv (t2, t3, t1)

120 continue

call ddeq (t1, ss)
return
end subroutine ddbesselyr

subroutine dddigammabe (nb2, berne, x, y)

!  This evaluates the digamma function, using asymptotic formula DLMF 5.11.2:
!  dig(x) ~ log(x) - 1/(2*x) - Sum_{k=1}^inf B[2k] / (2*k*x^(2*k)).
!  Before using this formula, the recursion dig(x+1) = dig(x) + 1/x is used
!  to shift the argument up by IQQ, where IQQ is set based on MPNW below.
!  The array berne contains precomputed even Bernoulli numbers (see MPBERNER
!  above). Its dimensions must be as shown below. NB2 must be greater than
!  1.4 x precision in decimal digits.

implicit none
integer, intent (in):: nb2
real (ddknd), intent(in):: berne(1:2,nb2), x(1:2)
real (ddknd), intent(out):: y(1:2)
real (ddknd), parameter:: dber = 1.4d0, dfrac = 0.4d0
integer k, i1, i2, ic1, iqq, n1
real (ddknd) d1
real (ddknd) f1(1:2), sum1(1:2), sum2(1:2), &
  t1(1:2), t2(1:2), t3(1:2), t4(1:2), &
  t5(1:2), xq(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

ddnw1 = min (ddnw + 1, ddnwx)
iqq = dfrac * ddnw1 * dddpw
call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw*ddnbt, eps)
call dddmc (1.d0, 0, f1)

!   Check if argument is less than or equal to 0 -- undefined.

if (ddsgn (x) <= 0) then
  write (ddldb, 2)
2 format ('*** DDDIGAMMABE: Argument is less than or equal to 0')
  call ddabrt
endif

!   Check if berne array has been initialized.

call ddmdc (berne(1:2,1), d1, n1)
d1 = d1 * 2.d0 ** n1
if (ddwprecr (berne(1:2,1)) < ddnw .or. &
  abs (d1 - 1.d0 / 6.d0) > ddrdfz .or. nb2 < int (dber * dddpw * ddnw)) then
  write (ddldb, 3) int (dber * dddpw * ddnw)
3 format ('*** DDDIGAMMABE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries using DDBERNE or DDBERNER.')
  call ddabrt
endif

! sum1 = mpreal (0.d0, nwds)
! sum2 = mpreal (0.d0, nwds)
! xq = x + dble (iqq)

call dddmc (0.d0, 0, sum1)
call dddmc (0.d0, 0, sum2)
call dddmc (dble (iqq), 0, t1)
call ddadd (x, t1, t2)
call ddeq (t2, xq)

do k = 0, iqq - 1
!   sum1 = sum1 + f1 / (x + dble (k))

  call dddmc (dble (k), 0, t1)
  call ddadd (x, t1, t2)
  call dddiv (f1, t2, t3)
  call ddadd (sum1, t3, t1)
  call ddeq (t1, sum1)
enddo

! t1 = mpreal (1.d0, nwds)
! t2 = xq ** 2

call dddmc (1.d0, 0, t1)
call ddmul (xq, xq, t2)

do k = 1, nb2
!  t1 = t1 * t2
!  t3 = bb(k) / (2.d0 * dble (k) * t1)
!  sum2 = sum2 + t3

  call ddmul (t1, t2, t3)
  call ddeq (t3, t1)
  call ddmuld (t1, 2.d0 * dble (k), t4)
  call dddiv (berne(1:2,k), t4, t3)
  call ddadd (sum2, t3, t4)
  call ddeq (t4, sum2)

!  if (abs (t3 / sum2) < eps) goto 100

  call ddabs (t3, tc1)
  call ddmul (eps, sum2, tc3)
  call ddabs (tc3, tc2)
  call ddcpr (tc1, tc2, ic1)
  if (ic1 <= 0) goto 110
enddo

write (ddldb, 4)
4 format ('*** DDDIGAMMABE: Loop end error: Increase NB2')
call ddabrt

110 continue

! digammax = -sum1 + log (xq) - 1.d0 / (2.d0 * xq) - sum2

call ddneg (sum1, t1)
call ddlog (xq, t2)
call ddadd (t1, t2, t3)
call ddmuld (xq, 2.d0, t4)
call dddiv (f1, t4, t5)
call ddsub (t3, t5, t2)
call ddsub (t2, sum2, t1)
call ddeq (t1, y)
return
end subroutine dddigammabe

subroutine dderfr (z, terf)

!   This evaluates the erf function, using a combination of two series.
!   In particular, the algorithm is (where B = (mpnw + 1) * mpnbt, and
!   dcon is a constant defined below):

!   if (z == 0) then
!     erf = 0
!   elseif (z > sqrt(B*log(2))) then
!     erf = 1
!   elseif (z < -sqrt(B*log(2))) then
!     erf = -1
!   elseif (abs(z) < B/dcon + 8) then
!     erf = 2 / (sqrt(pi)*exp(z^2)) * Sum_{k>=0} 2^k * z^(2*k+1)
!             / (1.3....(2*k+1))
!   else
!     erf = 1 - 1 / (sqrt(pi)*exp(z^2))
!             * Sum_{k>=0} (-1)^k * (1.3...(2*k-1)) / (2^k * z^(2*k+1))
!   endif

implicit none
real (ddknd), intent(in):: z(1:2)
real (ddknd), intent(out):: terf(1:2)
integer, parameter:: itrmx = 100000
real (ddknd), parameter:: dcon = 100.d0, pi = 3.1415926535897932385d0
integer ic1, ic2, ic3, k, nbt, n1
real (ddknd) d1, d2

real (ddknd) t1(1:2), t2(1:2), t3(1:2), t4(1:2), &
  t5(1:2), t6(1:2), t7(1:2), z2(1:2), tc1(1:2), &
  tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

ddnw1 = min (ddnw + 1, ddnwx)
call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw1*ddnbt, eps)

nbt = ddnw * ddnbt
d1 = aint (1.d0 + sqrt (nbt * log (2.d0)))
d2 = aint (nbt / dcon + 8.d0)
call dddmc (d1, 0, t1)
call dddmc (d2, 0, t2)
call ddcpr (z, t1, ic1)
! t1(2) = - t1(2)
call ddneg (t1, t3)
call ddeq (t3, t1)
call ddcpr (z, t1, ic2)
call ddcpr (z, t2, ic3)

if (ddsgn (z) == 0) then
  call dddmc (0.d0, 0, terf)
elseif (ic1 > 0) then
  call dddmc (1.d0, 0, terf)
elseif (ic2 < 0) then
  call dddmc (-1.d0, 0, terf)
elseif (ic3 < 0) then
  call ddmul (z, z, z2)
  call dddmc (0.d0, 0, t1)
  call ddeq (z, t2)
  call dddmc (1.d0, 0, t3)
  call dddmc (1.d10, 0, t5)

  do k = 0, itrmx
    if (k > 0) then
      call ddmuld (z2, 2.d0, t6)
      call ddmul (t6, t2, t7)
      call ddeq (t7, t2)
      d1 = 2.d0 * k + 1.d0
      call ddmuld (t3, d1, t6)
      call ddeq (t6, t3)
    endif

    call dddiv (t2, t3, t4)
    call ddadd (t1, t4, t6)
    call ddeq (t6, t1)
    call dddiv (t4, t1, t6)
    call ddcpr (t6, eps, ic1)
    call ddcpr (t6, t5, ic2)
    if (ic1 <= 0 .or. ic2 >= 0) goto 120
    call ddeq (t6, t5)
  enddo

write (ddldb, 3) 1, itrmx
3 format ('*** DDERFR: iteration limit exceeded',2i10)
call ddabrt

120 continue

  call ddmuld (t1, 2.d0, t3)
  call ddsqrt (ddpicon, t4)
  call ddexp (z2, t5)
  call ddmul (t4, t5, t6)
  call dddiv (t3, t6, t7)
  call ddeq (t7, terf)
else
  call ddmul (z, z, z2)
  call dddmc (0.d0, 0, t1)
  call dddmc (1.d0, 0, t2)
  call ddabs (z, t3)
  call dddmc (1.d10, 0, t5)

  do k = 0, itrmx
    if (k > 0) then
      d1 = -(2.d0 * k - 1.d0)
      call ddmuld (t2, d1, t6)
      call ddeq (t6, t2)
      call ddmul (t2, t3, t6)
      call ddeq (t6, t3)
    endif

    call dddiv (t2, t3, t4)
    call ddadd (t1, t4, t6)
    call ddeq (t6, t1)
    call dddiv (t4, t1, t6)
    call ddcpr (t6, eps, ic1)
    call ddcpr (t6, t5, ic2)
    if (ic1 <= 0 .or. ic2 >= 0) goto 130
    call ddeq (t6, t5)
  enddo

write (ddldb, 3) 2, itrmx
call ddabrt

130 continue

  call dddmc (1.d0, 0, t2)
  call ddsqrt (ddpicon, t3)
  call ddexp (z2, t4)
  call ddmul (t3, t4, t5)
  call dddiv (t1, t5, t6)
  call ddsub (t2, t6, t7)
  call ddeq (t7, terf)
  if (ddsgn (z) < 0) then
    call ddneg (terf, t6)
    call ddeq (t6, terf)
  endif
endif

return
end subroutine dderfr

subroutine dderfcr (z, terfc)

!   This evaluates the erfc function, using a combination of two series.
!   In particular, the algorithm is (where B = (mpnw + 1) * mpnbt, and
!   dcon is a constant defined below):

!   if (z == 0) then
!     erfc = 1
!   elseif (z > sqrt(B*log(2))) then
!     erfc = 0
!   elseif (z < -sqrt(B*log(2))) then
!     erfc = 2
!   elseif (abs(z) < B/dcon + 8) then
!     erfc = 1 - 2 / (sqrt(pi)*exp(z^2)) * Sum_{k>=0} 2^k * z^(2*k+1)
!               / (1.3....(2*k+1))
!   else
!     erfc = 1 / (sqrt(pi)*exp(z^2))
!             * Sum_{k>=0} (-1)^k * (1.3...(2*k-1)) / (2^k * z^(2*k+1))
!   endif

implicit none
real (ddknd), intent(in):: z(1:2)
real (ddknd), intent(out):: terfc(1:2)
integer, parameter:: itrmx = 100000
real (ddknd), parameter:: dcon = 100.d0, pi = 3.1415926535897932385d0
integer ic1, ic2, ic3, k, nbt, n1
real (ddknd) d1, d2
real (ddknd) t1(1:2), t2(1:2), t3(1:2), t4(1:2), &
  t5(1:2), t6(1:2), t7(1:2), z2(1:2), tc1(1:2), &
  tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

ddnw1 = min (ddnw + 1, ddnwx)
call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw1*ddnbt, eps)

nbt = ddnw * ddnbt
d1 = aint (1.d0 + sqrt (nbt * log (2.d0)))
d2 = aint (nbt / dcon + 8.d0)
call dddmc (d1, 0, t1)
call dddmc (d2, 0, t2)
call ddcpr (z, t1, ic1)
call ddneg (t1, t3)
call ddeq (t3, t1)
call ddcpr (z, t1, ic2)
call ddcpr (z, t2, ic3)

if (ddsgn (z) == 0) then
  call dddmc (1.d0, 0, terfc)
elseif (ic1 > 0) then
  call dddmc (0.d0, 0, terfc)
elseif (ic2 < 0) then
  call dddmc (2.d0, 0, terfc)
elseif (ic3 < 0) then
  call ddmul (z, z, z2)
  call dddmc (0.d0, 0, t1)
  call ddeq (z, t2)
  call dddmc (1.d0, 0, t3)
  call dddmc (1.d10, 0, t5)

  do k = 0, itrmx
    if (k > 0) then
      call ddmuld (z2, 2.d0, t6)
      call ddmul (t6, t2, t7)
      call ddeq (t7, t2)
      d1 = 2.d0 * k + 1.d0
      call ddmuld (t3, d1, t6)
      call ddeq (t6, t3)
    endif

    call dddiv (t2, t3, t4)
    call ddadd (t1, t4, t6)
    call ddeq (t6, t1)
    call dddiv (t4, t1, t6)
    call ddcpr (t6, eps, ic1)
    call ddcpr (t6, t5, ic2)
    if (ic1 <= 0 .or. ic2 >= 0) goto 120
    call ddeq (t6, t5)
  enddo

write (ddldb, 3) 1, itrmx
3 format ('*** DDERFCR: iteration limit exceeded',2i10)
call ddabrt

120 continue

  call dddmc (1.d0, 0, t2)
  call ddmuld (t1, 2.d0, t3)
  call ddsqrt (ddpicon, t4)
  call ddexp (z2, t5)
  call ddmul (t4, t5, t6)
  call dddiv (t3, t6, t7)
  call ddsub (t2, t7, t6)
  call ddeq (t6, terfc)
else
  call ddmul (z, z, z2)
  call dddmc (0.d0, 0, t1)
  call dddmc (1.d0, 0, t2)
  call ddabs (z, t3)
  call dddmc (1.d10, 0, t5)

  do k = 0, itrmx
    if (k > 0) then
      d1 = -(2.d0 * k - 1.d0)
      call ddmuld (t2, d1, t6)
      call ddeq (t6, t2)
      call ddmul (t2, t3, t6)
      call ddeq (t6, t3)
    endif

    call dddiv (t2, t3, t4)
    call ddadd (t1, t4, t6)
    call ddeq (t6, t1)
    call dddiv (t4, t1, t6)
    call ddcpr (t6, eps, ic1)
    call ddcpr (t6, t5, ic2)
    if (ic1 <= 0 .or. ic2 >= 0) goto 130
    call ddeq (t6, t5)
  enddo

write (ddldb, 3) 2, itrmx
call ddabrt

130 continue

  call ddsqrt (ddpicon, t3)
  call ddexp (z2, t4)
  call ddmul (t3, t4, t5)
  call dddiv (t1, t5, t6)
  if (ddsgn (z) < 0) then
    call dddmc (2.d0, 0, t2)
    call ddsub (t2, t6, t7)
    call ddeq (t7, t6)
  endif

  call ddeq (t6, terfc)
endif

return
end subroutine dderfcr

subroutine ddexpint (x, y)

!   This evaluates the exponential integral function Ei(x):
!   Ei(x) = - incgamma (0, -x)

implicit none
real (ddknd), intent(in):: x(1:2)
real (ddknd), intent(out):: y(1:2)
real (ddknd) t1(1:2), t2(1:2), t3(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

if (ddsgn (x) == 0) then
  write (ddldb, 2)
2 format ('*** DDEXPINT: argument is zero')
  call ddabrt
endif

call dddmc (0.d0, 0, t1)
call ddneg (x, t2)
call ddincgammar (t1, t2, t3)
call ddneg (t3, y)
return
end subroutine ddexpint

subroutine ddgammar (t, z)

!   This evaluates the gamma function, using an algorithm of R. W. Potter.
!   The argument t must not exceed 10^8 in size (this limit is set below),
!   must not be zero, and if negative must not be integer.

!   In the parameter statement below:
!     itrmx = limit of number of iterations in series; default = 100000.
!     con1 = 1/2 * log (10) to DP accuracy.
!     dmax = maximum size of input argument.

implicit none
real (ddknd), intent(in):: t(1:2)
real (ddknd), intent(out):: z(1:2)
integer, parameter:: itrmx = 100000
real (ddknd), parameter:: al2 = 0.69314718055994530942d0, dmax = 1d8, &
  pi = 3.1415926535897932385d0
integer i, i1, ic1, j, nt, n1, n2, n3
real (ddknd) alpha, d1, d2, d3
real (ddknd) f1(1:2), sum1(1:2), sum2(1:2), tn(1:2), &
  t1(1:2), t2(1:2), t3(1:2), t4(1:2), t5(1:2), &
  t6(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

ddnw1 = min (ddnw + 1, ddnwx)
call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw1*ddnbt, eps)

call ddmdc (t, d1, n1)
d1 = d1 * 2.d0**n1
call ddnint (t, t1)
call ddcpr (t, t1, ic1)
i1 = ddsgn (t)
if (i1 == 0 .or. d1 > dmax .or. (i1 < 0 .and. ic1 == 0)) then
  write (ddldb, 2) dmax
2 format ('*** DDGAMMAR: input argument must have absolute value <=',f10.0,','/ &
  'must not be zero, and if negative must not be an integer.')
  call ddabrt
endif

call dddmc (1.d0, 0, f1)

!   Find the integer and fractional parts of t.

call ddinfr (t, t2, t3)

if (ddsgn (t3) == 0) then

!   If t is a positive integer, then apply the usual factorial recursion.

  call ddmdc (t2, d2, n2)
  nt = d2 * 2.d0 ** n2
  call ddeq (f1, t1)

  do i = 2, nt - 1
    call ddmuld (t1, dble (i), t2)
    call ddeq (t2, t1)
  enddo

  call ddeq (t1, z)
  goto 120
elseif (ddsgn (t) > 0) then

!   Apply the identity Gamma[t+1] = t * Gamma[t] to reduce the input argument
!   to the unit interval.

  call ddmdc (t2, d2, n2)
  nt = d2 * 2.d0 ** n2
  call ddeq (f1, t1)
  call ddeq (t3, tn)

  do i = 1, nt
    call dddmc (dble (i), 0, t4)
    call ddsub (t, t4, t5)
    call ddmul (t1, t5, t6)
    call ddeq (t6, t1)
  enddo
else

!   Apply the gamma identity to reduce a negative argument to the unit interval.

  call ddsub (f1, t, t4)
  call ddinfr (t4, t3, t5)
  call ddmdc (t3, d3, n3)
  nt = d3 * 2.d0 ** n3

  call ddeq (f1, t1)
  call ddsub (f1, t5, t2)
  call ddeq (t2, tn)

  do i = 0, nt - 1
    call dddmc (dble (i), 0, t4)
    call ddadd (t, t4, t5)
    call dddiv (t1, t5, t6)
    call ddeq (t6, t1)
  enddo
endif

!   Calculate alpha = bits of precision * log(2) / 2, then take the next even
!   integer value, so that alpha/2 and alpha^2/4 can be calculated exactly in DP.

alpha = 2.d0 * aint (0.25d0 * (ddnw1 + 1) * ddnbt * al2 + 1.d0)
d2 = 0.25d0 * alpha**2
call ddeq (tn, t2)
call dddiv (f1, t2, t3)
call ddeq (t3, sum1)

!   Evaluate the series with t.

do j = 1, itrmx
  call dddmc (dble (j), 0, t6)
  call ddadd (t2, t6, t4)
  call ddmuld (t4, dble (j), t5)
  call dddiv (t3, t5, t6)
  call ddmuld (t6, d2, t3)
  call ddadd (sum1, t3, t4)
  call ddeq (t4, sum1)

  call ddabs (t3, tc1)
  call ddmul (eps, sum1, tc3)
  call ddabs (tc3, tc2)
  call ddcpr (tc1, tc2, ic1)
  if (ic1 <= 0) goto 100
enddo

write (ddldb, 3) 1, itrmx
3 format ('*** DDGAMMAR: iteration limit exceeded',2i10)
call ddabrt

100 continue

call ddneg (tn, t2)
call dddiv (f1, t2, t3)
call ddeq (t3, sum2)

!   Evaluate the same series with -t.

do j = 1, itrmx
  call dddmc (dble (j), 0, t6)
  call ddadd (t2, t6, t4)
  call ddmuld (t4, dble (j), t5)
  call dddiv (t3, t5, t6)
  call ddmuld (t6, d2, t3)
  call ddadd (sum2, t3, t4)
  call ddeq (t4, sum2)

  call ddabs (t3, tc1)
  call ddmul (eps, sum2, tc3)
  call ddabs (tc3, tc2)
  call ddcpr (tc1, tc2, ic1)
  if (ic1 <= 0) goto 110
enddo

write (ddldb, 3) 2, itrmx
call ddabrt

110 continue

!   Compute sqrt (pi * sum1 / (tn * sin (pi * tn) * sum2))
!   and (alpha/2)^tn terms. Also, multiply by the factor t1, from the
!   If block above.

call ddeq (ddpicon, t2)
call ddmul (t2, tn, t3)
call ddcssnr (t3, t4, t5)
call ddmul (t5, sum2, t6)
call ddmul (tn, t6, t5)
call ddmul (t2, sum1, t3)
call dddiv (t3, t5, t6)
call ddneg (t6, t4)
call ddeq (t4, t6)
call ddsqrt (t6, t2)
call dddmc (0.5d0 * alpha, 0, t3)
call ddlog (t3, t4)
call ddmul (tn, t4, t5)
call ddexp (t5, t6)
call ddmul (t2, t6, t3)
call ddmul (t1, t3, t4)

!   Round to mpnw words precision.

call ddeq (t4, z)

120 continue

return
end subroutine ddgammar

subroutine ddhurwitzzetan (is, aa, zz)

!   This returns the Hurwitz zeta function of IS and AA, using an algorithm from:
!   David H. Bailey and Jonathan M. Borwein, "Crandall's computation of the
!   incomplete gamma function and the Hurwitz zeta function with applications to
!   Dirichlet L-series," Applied Mathematics and Computation, vol. 268C (Oct 2015),
!   pg. 462-477, preprint at:
!   https://www.davidhbailey.com/dhbpapers/lerch.pdf
!   This is limited to IS >= 2 and 0 < AA < 1.

implicit none
integer, intent(in):: is
real (ddknd), intent(in):: aa(1:2)
real (ddknd), intent(out):: zz(1:2)
integer, parameter:: itrmax = 1000000
real (ddknd), parameter:: pi = 3.1415926535897932385d0
integer i1, ic1, ic2, ic3, k, n1
real (ddknd) d1, dk
real (ddknd) gs1(1:2), gs2(1:2), ss(1:2), &
  sum1(1:2), sum2(1:2), sum3(1:2), ss1(1:2), ss2(1:2), &
  ss3(1:2), ss4(1:2), s1(1:2), s2(1:2), s3(1:2), &
  t1(1:2), t2(1:2), t3(1:2), t4(1:2), t5(1:2), &
  t6(1:2), t7(1:2), t8(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

ddnw1 = min (ddnw + 1, ddnwx)
call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw1*ddnbt, eps)

if (is <= 0) then
  write (ddldb, 3)
3 format ('*** DDHURWITZZETAN: IS less than or equal to 0:')
  call ddabrt
endif

call dddmc (0.d0, 0, t1)
call dddmc (1.d0, 0, t2)
call ddcpr (aa, t1, ic1)
call ddcpr (aa, t2, ic2)
if (ic1 <= 0 .or. ic2 >= 0) then
  write (ddldb, 4)
4 format ('*** DDHURWITZZETAN: AA must be in the range (0, 1)')
  call ddabrt
endif

! ss = mpreal (dble (is), mpnw)
! ss1 = 0.5d0 * ss
! ss2 = 0.5d0 * (ss + 1.d0)
! ss3 = 0.5d0 * (1.d0 - ss)
! ss4 = 1.d0 - 0.5d0 * ss

call dddmc (dble(is), 0, ss)
call ddmuld (ss, 0.5d0, ss1)
call dddmc (1.d0, 0, t1)
call ddadd (t1, ss, t2)
call ddmuld (t2, 0.5d0, ss2)
call ddsub (t1, ss, t2)
call ddmuld (t2, 0.5d0, ss3)
call ddsub (t1, ss1, ss4)

! gs1 = gamma (ss1)
! gs2 = gamma (ss2)
! t1 = pi * aa ** 2

call ddgammar (ss1, gs1)
call ddgammar (ss2, gs2)
call ddmul (aa, aa, t2)
call ddmul (ddpicon, t2, t1)

! sum1 = (incgamma (ss1, t1) / gs1 + incgamma (ss2, t1) / gs2) / abs (aa)**is
! sum2 = mpreal (0.d0, mpnw)
! sum3 = mpreal (0.d0, mpnw)

call ddincgammar (ss1, t1, t2)
call dddiv (t2, gs1, t3)
call ddincgammar (ss2, t1, t2)
call dddiv (t2, gs2, t4)
call ddadd (t3, t4, t2)
call ddabs (aa, t3)
call ddnpwr (t3, is, t4)
call dddiv (t2, t4, sum1)
call dddmc (0.d0, 0, sum2)
call dddmc (0.d0, 0, sum3)

do k = 1, itrmax
  dk = dble (k)

!  t1 = pi * (dk + aa)**2
!  t2 = pi * (-dk + aa)**2
!  t3 = dk**2 * pi
!  t4 = 2.d0 * pi * dk * aa

  call dddmc (dk, 0, t5)
  call ddadd (t5, aa, t6)
  call ddmul (t6, t6, t5)
  call ddmul (ddpicon, t5, t1)
  call dddmc (-dk, 0, t5)
  call ddadd (t5, aa, t6)
  call ddmul (t6, t6, t7)
  call ddmul (ddpicon, t7, t2)
  call ddmul (t5, t5, t6)
  call ddmul (ddpicon, t6, t3)
  call ddmuld (ddpicon, 2.d0 * dk, t5)
  call ddmul (t5, aa, t4)

!  s1 = (incgamma (ss1, t1) / gs1 + incgamma (ss2, t1) / gs2) / abs (dk + aa)**is

  call ddincgammar (ss1, t1, t5)
  call dddiv (t5, gs1, t6)
  call ddincgammar (ss2, t1, t5)
  call dddiv (t5, gs2, t7)
  call ddadd (t6, t7, t5)
  call dddmc (dk, 0, t6)
  call ddadd (t6, aa, t7)
  call ddabs (t7, t6)
  call ddnpwr (t6, is, t7)
  call dddiv (t5, t7, s1)

!  s2 = (incgamma (ss1, t2) / gs1 - incgamma (ss2, t2) / gs2) / abs (-dk + aa)**is

  call ddincgammar (ss1, t2, t5)
  call dddiv (t5, gs1, t6)
  call ddincgammar (ss2, t2, t5)
  call dddiv (t5, gs2, t7)
  call ddsub (t6, t7, t5)
  call dddmc (-dk, 0, t6)
  call ddadd (t6, aa, t7)
  call ddabs (t7, t6)
  call ddnpwr (t6, is, t7)
  call dddiv (t5, t7, s2)

!  sum1 = sum1 + s1
!  sum2 = sum2 + s2

  call ddadd (sum1, s1, t5)
  call ddeq (t5, sum1)
  call ddadd (sum2, s2, t5)
  call ddeq (t5, sum2)

!  s3 = (incgamma (ss3, t3) * cos (t4)/ gs1 + incgamma (ss4, t3) * sin (t4) / gs2) &
!    / mpreal (dk, mpnw)**(1-is)

  call ddincgammar (ss3, t3, t5)
  call ddcssnr (t4, t6, t7)
  call ddmul (t5, t6, t8)
  call dddiv (t8, gs1, t5)
  call ddincgammar (ss4, t3, t6)
  call ddmul (t6, t7, t8)
  call dddiv (t8, gs2, t6)
  call ddadd (t5, t6, t7)
  call dddmc (dk, 0, t5)
  i1 = 1 - is
  call ddnpwr (t5, i1, t6)
  call dddiv (t7, t6, s3)

!  sum3 = sum3 + s3

  call ddadd (sum3, s3, t5)
  call ddeq (t5, sum3)

!  if (abs (s1) < eps * abs (sum1) .and. abs (s2) < eps * abs (sum2) .and. &
!    abs (s3) < eps * abs (sum3)) goto 100

  call ddabs (s1, tc1)
  call ddmul (eps, sum1, tc3)
  call ddabs (tc3, tc2)
  call ddcpr (tc1, tc2, ic1)
  call ddabs (s2, tc1)
  call ddmul (eps, sum2, tc3)
  call ddabs (tc3, tc2)
  call ddcpr (tc1, tc2, ic2)
  call ddabs (s3, tc1)
  call ddmul (eps, sum3, tc3)
  call ddabs (tc3, tc2)
  call ddcpr (tc1, tc2, ic3)
  if (ic1 <= 0 .and. ic2 <= 0 .and. ic3 <= 0) goto 100
enddo

write (ddldb, 5)
5 format ('*** DDHURWITZZETAN: Loop end error')
call ddabrt

100 continue

if (mod (is, 2) == 0) then
!  t1 = pi ** (is / 2) / ((ss - 1.d0) * gamma (ss1))

  i1 = is / 2
  call ddnpwr (ddpicon, i1, t2)
  call dddmc (1.d0, 0, t3)
  call ddsub (ss, t3, t4)
  call ddgammar (ss1, t5)
  call ddmul (t4, t5, t6)
  call dddiv (t2, t6, t1)
else
!  t1 = sqrt (pi) * pi ** ((is - 1) / 2) / ((ss - 1.d0) * gamma (ss1))

  i1 = (is - 1) / 2
  call ddnpwr (ddpicon, i1, t2)
  call ddsqrt (ddpicon, t3)
  call ddmul (t2, t3, t4)
  call dddmc (1.d0, 0, t2)
  call ddsub (ss, t2, t3)
  call ddgammar (ss1, t5)
  call ddmul (t3, t5, t6)
  call dddiv (t4, t6, t1)
endif

!t2 = pi ** is / sqrt (pi)

call ddnpwr (ddpicon, is, t3)
call ddsqrt (ddpicon, t4)
call dddiv (t3, t4, t2)

! hurwitzzetan = t1 + 0.5d0 * sum1 + 0.5d0 * sum2 + t2 * sum3

call ddmuld (sum1, 0.5d0, t3)
call ddmuld (sum2, 0.5d0, t4)
call ddmul (sum3, t2, t5)
call ddadd (t1, t3, t6)
call ddadd (t6, t4, t7)
call ddadd (t7, t5, t1)
call ddeq (t1, zz)

return
end subroutine ddhurwitzzetan

subroutine ddhurwitzzetanbe (nb2, berne, iss, aa, zz)

!  This evaluates the Hurwitz zeta function, using the combination of
!  the definition formula (for large iss), and an Euler-Maclaurin scheme
!  (see formula 25.2.9 of the DLMF). The array berne contains precomputed
!  even Bernoulli numbers (see MPBERNER above). Its dimensions must be as
!  shown below. NB2 must be greater than 1.4 x precision in decimal digits.

implicit none
integer, intent(in):: nb2, iss     
real (ddknd), intent(in):: berne(1:2,nb2), aa(1:2)
real (ddknd), intent(out):: zz(1:2)
integer, parameter:: itrmax = 1000000
real (ddknd), parameter:: dber = 1.4d0, dcon = 0.6d0
integer i1, i2, ic1, iqq, k, n1
real (ddknd) d1, dp
real (ddknd) aq(1:2), aq2(1:2), s1(1:2), s2(1:2), &
  s3(1:2), s4(1:2), t1(1:2), t2(1:2), t3(1:2), &
  t4(1:2), t5(1:2), t6(1:2), eps(1:2), f1(1:2), tc1(1:2), &
  tc2(1:2), tc3(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

!   Check if berne array has been initialized.

call ddmdc (berne(1:2,1), d1, n1)
d1 = d1 * 2.d0 ** n1
if (ddwprecr (berne(1:2,1)) < ddnw .or. &
  abs (d1 - 1.d0 / 6.d0) > ddrdfz .or. nb2 < int (dber * dddpw * ddnw)) then
  write (ddldb, 3) int (dber * dddpw * ddnw)
3 format ('*** DDHURWITZZETANBE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries by calling DDBERNE or DDBERNER.')
  call ddabrt
endif

ddnw1 = min (ddnw + 1, ddnwx)
call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw1*ddnbt, eps)
call dddmc (1.d0, 0, f1)
call dddmc (0.d0, 0, s1)
call dddmc (0.d0, 0, s2)
call dddmc (0.d0, 0, s3)
call dddmc (0.d0, 0, s4)

if (iss <= 0) then
  write (ddldb, 4)
4 format ('*** DDHURWITZZETANBE: ISS <= 0')
  call ddabrt
endif

if (ddsgn (aa) < 0) then
  write (ddldb, 5)
5 format ('*** DDHURWITZZETANBE: AA < 0')
  call ddabrt
endif

dp = anint (ddnw1 * dddpw)

!   If iss > a certain value, then use definition formula.

if (iss > 2.303d0 * dp / log (2.515d0 * dp)) then
  do k = 0, itrmax
!    t1 = 1.d0 / (aa + dble (k))**iss
!    s1 = s1 + t1

    call dddmc (dble (k), 0, t1)
    call ddadd (aa, t1, t2)
    call ddnpwr (t2, iss, t3)
    call dddiv (f1, t3, t1)
    call ddadd (s1, t1, t2)
    call ddeq (t2, s1)

!    if (abs (t1 / s1) < eps) goto 110

    call ddabs (t1, tc1)
    call ddmul (eps, s1, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 110
  enddo

  write (6, 6)
6 format ('*** DDHURWITZZETANBE: Loop end error 1')
  call ddabrt
endif

call ddmdc (aa, d1, n1)
d1 = d1 * 2.d0**n1
iqq = max (dcon * ddnw1 * dddpw - d1, 0.d0)

do k = 0, iqq - 1
!  s1 = s1 + 1.d0 / (aa + dble (k))**iss

  call dddmc (dble (k), 0, t1)
  call ddadd (aa, t1, t2)
  call ddnpwr (t2, iss, t3)
  call dddiv (f1, t3, t1)
  call ddadd (s1, t1, t2)
  call ddeq (t2, s1)
enddo

! aq = aa + dble (iqq)

call dddmc (dble (iqq), 0, t1)
call ddadd (aa, t1, aq)

! s2 = 1.d0 / (dble (iss - 1) * aq**(iss -  1))

call dddmc (dble (iss - 1), 0, t1)
call ddnpwr (aq, iss - 1, t2)
call ddmul (t1, t2, t3)
call dddiv (f1, t3, s2)

! s3 = 1.d0 / (2.d0 * aq**iss)

call ddnpwr (aq, iss, t1)
call ddmuld (t1, 2.d0, t2)
call dddiv (f1, t2, s3)

! t1 = mpreal (dble (iss), nwds)
! t2 = mpreal (1.d0, nwds)
! t3 = aq**(iss - 1)
! aq2 = aq**2

call dddmc (dble (iss), 0, t1)
call dddmc (1.d0, 0, t2)
call ddnpwr (aq, iss - 1, t3)
call ddmul (aq, aq, aq2)

do k = 1, nb2
!  if (k > 1) t1 = t1 * dble (iss + 2*k - 3) * dble (iss + 2*k - 2)

  if (k > 1) then
    call ddmuld (t1, dble (iss + 2*k - 3), t5)
    call ddmuld (t5, dble (iss + 2*k - 2), t1)
  endif

!  t2 = t2 * dble (2 * k - 1) * dble (2 * k)

  call ddmuld (t2, dble (2 * k - 1), t5)
  call ddmuld (t5, dble (2 * k), t2)

!  t3 = t3 * aq2
!  t4 = rb(k) * t1 / (t2 * t3)
!  s4 = s4 + t4

  call ddmul (t3, aq2, t5)
  call ddeq (t5, t3)
  call ddmul (berne(1:2,k), t1, t5)
  call ddmul (t2, t3, t6)
  call dddiv (t5, t6, t4)
  call ddadd (s4, t4, t5)
  call ddeq (t5, s4)

!  if (abs (t4) < eps) goto 110

  call ddabs (t4, tc1)
  call ddmul (eps, s4, tc3)
  call ddabs (tc3, tc2)
  call ddcpr (tc1, tc2, ic1)
  if (ic1 <= 0) goto 110
enddo

write (6, 7)
7 format ('*** DDHURWITZZETANBE: End loop error 2; call DDBERNE with larger NB.')
call ddabrt

110 continue

! hurwitz_be = s1 + s2 + s3 + s4

call ddadd (s1, s2, t5)
call ddadd (t5, s3, t6)
call ddadd (t6, s4, s1)
call ddeq (s1, zz)

return
end subroutine ddhurwitzzetanbe

subroutine ddhypergeompfq (np, nq, aa, bb, xx, yy)

!  This returns the HypergeometricPFQ function, namely the sum of the infinite series

!  Sum_0^infinity poch(aa(1),n)*poch(aa(2),n)*...*poch(aa(np),n) /
!      poch(bb(1),n)*poch(bb(2),n)*...*poch(bb(nq),n) * xx^n / n!

!  This subroutine evaluates the HypergeometricPFQ function directly according to
!  this definitional formula. The arrays aa and bb must be dimensioned as shown below.
!  NP and NQ are limited to [1,10].

implicit none
integer, intent(in):: np, nq
real (ddknd), intent(in):: aa(1:2,np), bb(1:2,nq), xx(1:2)
real (ddknd), intent(out):: yy(1:2)
integer, parameter:: itrmax = 1000000, npq = 10
integer i1, i2, ic1, j, k
real (ddknd) sum(1:2), td(1:2), tn(1:2), t1(1:2), &
  t2(1:2), t3(1:2), t4(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

ddnw1 = min (ddnw + 1, ddnwx)
call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw1*ddnbt, eps)

if (np < 1 .or. np > npq .or. nq < 1 .or. nq > npq) then
  write (ddldb, 2) npq
2 format ('*** DDHYPERGEOMPFQ: NP and NQ must be between 1 and',i4)
  call ddabrt
endif

call dddmc (1.d0, 0, sum)
call dddmc (1.d0, 0, td)
call dddmc (1.d0, 0, tn)

do k = 1, itrmax
  call dddmc (dble (k - 1), 0, t1)

  do j = 1, np
    call ddadd (t1, aa(1:2,j), t2)
    call ddmul (tn, t2, t3)
    call ddeq (t3, tn)
  enddo

  do j = 1, nq
    call ddadd (t1, bb(1:2,j), t2)
    call ddmul (td, t2, t3)
    call ddeq (t3, td)
  enddo

  call ddmul (tn, xx, t2)
  call ddeq (t2, tn)
  call ddmuld (td, dble (k), t3)
  call ddeq (t3, td)
  call dddiv (tn, td, t1)
  call ddadd (sum, t1, t2)
  call ddeq (t2, sum)

  call ddabs (t1, tc1)
  call ddmul (eps, sum, tc3)
  call ddabs (tc3, tc2)
  call ddcpr (tc1, tc2, ic1)
  if (ic1 <= 0) goto 100
enddo

    write (ddldb, 3) itrmax
3   format ('*** DDHYPERGEOMPFQ: Loop end error',i10)
    call ddabrt

100  continue

call ddeq (sum, yy)
return
end subroutine ddhypergeompfq

subroutine ddincgammar (s, z, g)

!  This returns the incomplete gamma function, using a combination of formula
!  8.7.3 of the DLMF (for modest-sized z), formula 8.11.2 (for large z),
!  a formula from the Wikipedia page for the case S = 0, and another formula
!  from the Wikipedia page for the case S = negative integer. The formula
!  for the case S = 0 requires increased working precision, up to 2.5X normal,
!  depending on the size of Z.

implicit none
real (ddknd), intent(in):: s(1:2), z(1:2)
real (ddknd), intent(out):: g(1:2)
integer, parameter:: itrmax = 1000000
real (ddknd), parameter:: dmax = 0.833d0, egam = 0.5772156649015328606d0
integer ic1, k, nn, n1, n2
real (ddknd) d1, d2, bits
real (ddknd) t0(1:2), t1(1:2), t2(1:2), &
  t3(1:2), t4(1:2), t5(1:2), f1(1:2), &
  tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

n1 = ddsgn (s)
n2 = ddsgn (z)
if (n2 == 0 .or. n1 /= 0 .and. n2 < 0) then
  write (ddldb, 2)
2 format ('*** DDINCGAMMAR: The second argument must not be zero,'/ &
    'and must not be negative unless the first is zero.')
  call ddabrt
endif

ddnw1 = min (ddnw + 1, ddnwx)
call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw1*ddnbt, eps)

call dddmc (1.d0, 0, f1)
call ddmdc (z, d1, n1)
d1 = d1 * 2.d0 ** n1
bits = ddnw1 * ddnbt

if (abs (d1) < dmax * bits) then

!   This is for modest-sized z.

  call ddinfr (s, t1, t2)
  call ddcpr (s, t1, ic1)
  call ddmdc (s, d2, n2)
  nn = d2 * 2.d0**n2

  if (ic1 == 0 .and. nn == 1) then

!   S = 1; result is exp (-z).

    call ddneg (z, t0)
    call ddexp (t0, t1)
    goto 200
  elseif (ic1 == 0 .and. nn <= 0) then

!    S is zero or a negative integer -- use a different algorithm. In
!    either event, first compute incgamma for S = 0. For large Z, the
!    working precision must be increased, up to 2.5X times normal.

!    mpnw2 = min (mpnw1 + 1.5d0 * d1 / (dmax * bits) * mpnw, 5*mpnw/2+1.d0)
    ddnw2 = ddnw1
    call ddeq (z, t0)
    call ddeq (z, t1)
    call dddmc (1.d0, 0, t2)

    do k = 2, itrmax
      if (mod (k, 2) == 1) then
        d1 = dble (k)
        call dddivd (f1, d1, t3)
        call ddadd (t2, t3, t4)
        call ddeq (t4, t2)
      endif
      call ddmul (z, t1, t3)
      d1 = 2.d0 * dble (k)
      call dddivd (t3, d1, t1)
      call ddmul (t1, t2, t3)
      call ddadd (t0, t3, t4)
      call ddeq (t4, t0)

      call ddabs (t3, tc1)
      call ddmul (eps, t0, tc3)
      call ddabs (tc3, tc2)
      call ddcpr (tc1, tc2, ic1)
      if (ic1 <= 0) goto 100
    enddo

    write (ddldb, 4)
4   format ('*** DDINCGAMMAR: Loop end error 1')
    call ddabrt

100  continue

    call ddneg (ddegammacon, t1)
    call ddabs (z, t3)
    call ddlog (t3, t2)
    call ddsub (t1, t2, t3)
    call ddmuld (z, -0.5d0, t4)
    call ddexp (t4, t5)
    call ddmul (t5, t0, t4)
    call ddadd (t3, t4, t1)
    if (nn == 0) goto 200

!   S is negative integer (not zero).

    nn = abs (nn)
    call dddmc (1.d0, 0, t0)
    call ddeq (t0, t2)

    do k = 1, nn - 1
      call ddmuld (t0, dble (k), t2)
      call ddeq (t2, t0)
    enddo

    call ddmuld (t0, dble (nn), t5)

    do k = 1, nn - 1
      call ddmul (t2, z, t3)
      call dddivd (t3, dble (nn - k), t4)
      call ddneg (t4, t2)
      call ddadd (t0, t2, t3)
      call ddeq (t3, t0)
    enddo

    call ddexp (z, t2)
    call dddiv (t0, t2, t3)
    call ddnpwr (z, nn, t4)
    call dddiv (t3, t4, t2)

    if (mod (nn, 2) == 0) then
      call ddadd (t2, t1, t3)
    else
      call ddsub (t2, t1, t3)
    endif
    call dddiv (t3, t5, t1)
    goto 200
  endif

  call ddgammar (s, t1)
  call ddmul (s, t1, t3)
  call dddiv (f1, t3, t2)
  call ddeq (t2, t0)

  do k = 1, itrmax
    call ddmul (t2, z, t5)
    call dddmc (dble (k), 0, t3)
    call ddadd (s, t3, t4)
    call dddiv (t5, t4, t2)
    call ddadd (t0, t2, t3)
    call ddeq (t3, t0)

    call ddabs (t2, tc1)
    call ddmul (eps, t0, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 110
  enddo

  write (ddldb, 5) itrmax
5   format ('*** DDINCGAMMAR: Loop end error 1')
  call ddabrt

110 continue

  call ddpower (z, s, t2)
  call ddexp (z, t3)
  call dddiv (t2, t3, t4)
  call ddmul (t4, t0, t5)
  call ddsub (f1, t5, t2)
  call ddmul (t1, t2, t3)
  call ddeq (t3, t1)
  goto 200
else

!   This is for large z. Note that if S is a positive integer, this loop
!   is finite.

  call dddmc (1.d0, 0, t0)
  call dddmc (1.d0, 0, t1)

  do k = 1, itrmax
    call dddmc (dble (k), 0, t2)
    call ddsub (s, t2, t3)
    call ddmul (t1, t3, t4)
    call dddiv (t4, z, t1)
    call ddadd (t0, t1, t2)
    call ddeq (t2, t0)

    call ddabs (t1, tc1)
    call ddmul (eps, t0, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 120
  enddo

  write (ddldb, 6)
6 format ('*** DDINCGAMMAR: Loop end error 3')
  call ddabrt

120 continue

   call ddsub (s, f1, t2)
   call ddpower (z, t2, t3)
   call ddexp (z, t4)
   call dddiv (t3, t4, t2)
   call ddmul (t2, t0, t1)
   goto 200
endif

200 continue

call ddeq (t1, g)

return
end subroutine ddincgammar

subroutine ddpolygamma (nn, x, y)

!   This returns polygamma (nn, x) for nn >= 0 and 0 < x < 1, by calling
!   mphurwitzzetan.

implicit none
integer, intent(in):: nn
real (ddknd), intent(in):: x(1:2)
real (ddknd), intent(out):: y(1:2)
integer ic1, ic2, k
real (ddknd) t1(1:2), t2(1:2), t3(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

ddnw1 = min (ddnw, ddnwx)

if (nn <= 0) then
  write (ddldb, 2)
2 format ('*** DDPOLYGAMMA: NN <= 0')
  call ddabrt
endif

call dddmc (0.d0, 0, t1)
call dddmc (1.d0, 0, t2)
call ddcpr (x, t1, ic1)
call ddcpr (x, t2, ic2)
if (ic1 <= 0 .or. ic2 >= 0) then
  write (ddldb, 3)
3 format ('*** DDPOLYGAMMA: X must be in the range (0, 1)')
  call ddabrt
endif

call dddmc (1.d0, 0, t1)

do k = 1, nn
  call ddmuld (t1, dble(k), t2)
  call ddeq (t2, t1)
enddo

if (mod (nn + 1, 2) == 1) then
  call ddneg (t1, t2)
  call ddeq (t2, t1)
endif
call ddhurwitzzetan (nn + 1, x, t2)
call ddmul (t1, t2, t3)
call ddeq (t3, y)

return
end subroutine ddpolygamma

subroutine ddpolygammabe (nb2, berne, nn, x, y)

!  This returns polygamma (nn, x) for nn >= 0, by calling mphurwitzzetanbe.
!  The array berne contains precomputed even Bernoulli numbers (see MPBERNER
!  above). Its dimensions must be as shown below. NB2 must be greater than
!  1.4 x precision in decimal digits.

implicit none
integer, intent(in):: nb2, nn
real (ddknd), intent(in):: berne(1:2,nb2), x(1:2)
real (ddknd), intent(out):: y(1:2)
real (ddknd), parameter:: dber = 1.4d0
integer i1, i2, k, n1
real (ddknd) d1
real (ddknd) t1(1:2), t2(1:2), t3(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

!   Check if berne array has sufficient entries.

call ddmdc (berne(1:2,1), d1, n1)
d1 = d1 * 2.d0 ** n1
if (ddwprecr (berne(1:2,1)) < ddnw .or. &
  abs (d1 - 1.d0 / 6.d0) > ddrdfz .or. nb2 < int (dber * dddpw * ddnw)) then
  write (ddldb, 3) int (dber * dddpw * ddnw)
3 format ('*** DDPOLYGAMMABE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries by calling DDBERNE or DDBERER.')
  call ddabrt
endif

ddnw1 = min (ddnw, ddnwx)

if (nn <= 0) then
  write (ddldb, 4)
4 format ('*** DDPOLYGAMMABE: NN <= 0')
  call ddabrt
endif

if (ddsgn (x) < 0) then
  write (ddldb, 5)
5 format ('*** DDPOLYGAMMABE: X < 0')
  call ddabrt
endif

call dddmc (1.d0, 0, t1)

do k = 1, nn
  call ddmuld (t1, dble(k), t2)
  call ddeq (t2, t1)
enddo

if (mod (nn + 1, 2) == 1) then
  call ddneg (t1, t2)
  call ddeq (t2, t1)
endif
call ddhurwitzzetanbe (nb2, berne, nn + 1, x, t2)
call ddmul (t1, t2, t3)
call ddeq (t3, y)

return
end subroutine ddpolygammabe

subroutine ddpolylogini (nn, arr)

!   Initializes the MP array arr with data for mppolylogneg.
!   NN must be in the range (-nmax, -1).

implicit none
integer, intent(in):: nn
real (ddknd), intent(out):: arr(1:2,1:abs(nn))
integer, parameter:: nmax = 1000
integer i1, i2, k, n, nna
real (ddknd) aa(1:2,2,abs(nn)), t1(1:2), t2(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

nna = abs (nn)
ddnw1 = min (ddnw + 1, ddnwx)
i1 = 2
i2 = 1
call dddmc (1.d0, 0, aa(1:2,1,1))
call dddmc (1.d0, 0, aa(1:2,2,1))

do k = 2, nna
  call dddmc (0.d0, 0, aa(1:2,1,k))
  call dddmc (0.d0, 0, aa(1:2,2,k))
enddo

do n = 2, nna
  i1 = 3 - i1
  i2 = 3 - i1

  do k = 2, n
    call ddmuld (aa(1:2,i1,k-1), dble (n + 1 - k), t1)
    call ddmuld (aa(1:2,i1,k), dble (k), t2)
    call ddadd (t1, t2, aa(1:2,i2,k))
  enddo
enddo

do k = 1, nna
  call ddeq (aa(1:2,i2,k), arr(1:2,k))
enddo

return
end subroutine ddpolylogini

subroutine ddpolylogneg (nn, arr, x, y)

!   This returns polylog (nn, x) for the case nn < 0. Before calling this,
!   one must call mppolylognini to initialize the array arr for this NN.
!   The dimensions of arr must be as shown below.
!   NN must be in the range (-nmax, -1).
!   The parameter nmxa is the maximum number of additional words of
!   precision needed to overcome cancelation errors when x is negative,
!   for nmax = 1000.

implicit none
integer, intent(in):: nn
real (ddknd), intent(in):: arr(1:2,1:abs(nn)), x(1:2)
real (ddknd), intent(out):: y(1:2)
integer, parameter:: nmax = 1000, nmxa = 8525 / ddnbt + 1
integer i1, i2, k, n1, n2, nna
real (ddknd) d1, d2
real (ddknd) t1(1:2+nmxa), t2(1:2+nmxa), t3(1:2+nmxa), &
  t4(1:2+nmxa)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

nna = abs (nn)
call ddmdc (arr(1:2,1), d1, n1)
d1 = d1 * 2.d0 ** n1
call ddmdc (arr(1:2,nna), d2, n2)
d2 = d2 * 2.d0 ** n2

if (d1 /= 1.d0 .or. d2 /= 1.d0) then
  write (ddldb, 2)
2 format ('*** DDPOLYLOGNEG: Uninitialized or inadequately sized arrays'/ &
  'Call ddpolylogini or polylog_ini to initialize array. See documentation.')
  call ddabrt
endif

ddnw1 = min (ddnw + 1, ddnwx)

if (ddsgn (x) < 0) then
  i1 = (nna + 1) / 2
  call ddmdc (arr(1:2,i1), d1, n1)
  ddnw1 = min (ddnw1 + (n1 + 1) / ddnbt + 1, ddnwx)
endif

call ddeq (x, t1)
call ddeq (t1, t2)

do k = 2, nna
  call ddmul (x, t1, t3)
  call ddeq (t3, t1)
  call ddmul (arr(1:2,k), t1, t3)
  call ddadd (t2, t3, t4)
  call ddeq (t4, t2)
enddo

call dddmc (1.d0, 0, t3)
call ddsub (t3, x, t4)
call ddnpwr (t4, nna + 1, t3)
call dddiv (t2, t3, t4)
call ddeq (t4, y)

return
end subroutine ddpolylogneg

subroutine ddpolylogpos (nn, x, y)

!   This returns polylog (nn, x) for the case nn >= 0.

implicit none
integer, intent(in):: nn
real (ddknd), intent(in):: x(1:2)
real (ddknd), intent(out):: y(1:2)
integer, parameter:: itrmax = 1000000
integer ic1, k
real (ddknd) t1(1:2), t2(1:2), t3(1:2), t4(1:2), &
  t5(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps (1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

if (nn < 0) then
  write (ddldb, 1)
1 format ('*** DDPOLYLOGPOS: N is less than zero.'/ &
  'For negative n, call ddpolylogneg or polylog_neg. See documentation.')
  call ddabrt
endif

ddnw1 = min (ddnw + 1, ddnwx)
call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw1*ddnbt, eps)

call ddabs (x, t1)
call dddmc (1.d0, 0, t2)
call ddcpr (t1, t2, ic1)
if (ic1 >= 0) then
  write (ddldb, 3)
3 format ('*** DDPOLYLOGPOS: |X| must be less than one.')
  call ddabrt
endif

if (nn == 0) then
  call dddmc (1.d0, 0, t1)
  call ddsub (t1, x, t2)
  call dddiv (x, t2, t3)
  call ddeq (t3, y)
else
  call ddeq (x, t1)
  call ddeq (x, t2)

  do k = 2, itrmax
    call ddmul (x, t2, t3)
    call ddeq (t3, t2)
    call dddmc (dble (k), 0, t3)
    call ddnpwr (t3, nn, t4)
    call dddiv (t2, t4, t3)
    call ddadd (t1, t3, t4)
    call ddeq (t4, t1)

    call ddabs (t3, tc1)
    call ddmul (eps, t1, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 100
  enddo

  write (ddldb, 4)
4 format ('*** DDPOLYLOGPOS: Loop end error')
  call ddabrt

100 continue

  call ddeq (t1, y)
endif

return
end subroutine ddpolylogpos

subroutine ddstruvehn (nu, ss, zz)

!   This returns the StruveH function with integer arg NU and MPFR argument SS.

implicit none
integer, intent(in):: nu
real (ddknd), intent(in):: ss(1:2)
real (ddknd), intent(out):: zz(1:2)
integer, parameter:: itrmax = 1000000
real (ddknd), parameter:: dmax = 1000.d0, pi = 3.1415926535897932385d0
integer ic1, k, n1
real (ddknd) d1
real (ddknd) sum(1:2), td1(1:2), td2(1:2), tn1(1:2), &
  tnm1(1:2), t1(1:2), t2(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

if (nu < 0) then
  write (ddldb, 2)
2 format ('*** DDSTRUVEHN: NU < 0')
  call ddabrt
endif

call ddmdc (ss, d1, n1)
d1 = abs (d1) * 2.d0**n1
if (d1 > dmax) then
  write (ddldb, 3)
3 format ('*** DDSTRUVEHN: ABS(SS) >',f8.2)
  call ddabrt
endif

ddnw1 = min (ddnw * (1.d0 + d1 / dmax), dble (ddnwx))

call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw*ddnbt, eps)

! tn1 = mpreal (1.d0, nwds1)
! tnm1 = -0.25d0 * mpreal (ss, nwds1)**2
! td1 = 0.5d0 * sqrt (mppi (nwds1))
! td2 = td1

call dddmc (1.d0, 0, tn1)
call ddmul (ss, ss, t1)
call ddmuld (t1, -0.25d0, tnm1)
call ddsqrt (ddpicon, t1)
call ddmuld (t1, 0.5d0, td1)
call ddeq (td1, td2)

! do k = 1, nu
!  td2 = (k + 0.5d0) * td2
! enddo

do k = 1, nu
  d1 = k + 0.5d0
  call ddmuld (td2, d1, t1)
  call ddeq (t1, td2)
enddo

! sum = tn1 / (td1 * td2)

call ddmul (td1, td2, t2)
call dddiv (tn1, t2, sum)

do k = 1, itrmax

!  tn1 = tnm1 * tn1
!  td1 = (k + 0.5d0) * td1
!  td2 = (nu + k + 0.5d0) * td2
!  t1 = tn1 / (td1 * td2)
!  sum = sum + t1

  call ddmul (tnm1, tn1, t1)
  call ddeq (t1, tn1)
  d1 = k + 0.5d0
  call ddmuld (td1, d1, t1)
  call ddeq (t1, td1)
  d1 = nu + k + 0.5d0
  call ddmuld (td2, d1, t1)
  call ddeq (t1, td2)
  call ddmul (td1, td2, t2)
  call dddiv (tn1, t2, t1)
  call ddadd (sum, t1, t2)
  call ddeq (t2, sum)

!  if (abs (t1) < eps) goto 100

  call ddabs (t1, tc1)
  call ddmul (eps, sum, tc3)
  call ddabs (tc3, tc2)
  call ddcpr (tc1, tc2, ic1)
  if (ic1 <= 0) goto 100
enddo

write (ddldb, 5)
5 format ('*** DDSTRUVEHN: Loop end error')
call ddabrt

100 continue

! struvehn = (0.5d0 * ss)**(nu + 1) * sum

call ddmuld (ss, 0.5d0, t1)
n1 = nu + 1
call ddnpwr (t1, n1, t2)
call ddmul (t2, sum, t1)
call ddeq (t1, zz)
return
end subroutine ddstruvehn

subroutine ddzetar (ss, zz)

!   This returns the zeta function of an MPR argument SS using an algorithm
!   due to Peter Borwein.

implicit none
real (ddknd), intent(in):: ss(1:2)
real (ddknd), intent(out):: zz(1:2)
integer, parameter:: itrmax = 1000000
real (ddknd), parameter:: dfrac = 1.d0+dddpw, pi = 3.1415926535897932385d0
integer i, ic1, iss, j, n, n1, n2
real (ddknd) d1, d2
real (ddknd) f1(1:2), s(1:2), t1(1:2), t2(1:2), &
  t3(1:2), t4(1:2), t5(1:2), tn(1:2), tt(1:2), &
  tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)
real (ddknd) sgn

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

ddnw1 = min (ddnw + 1, ddnwx)
call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw1*ddnbt, eps)

call dddmc (1.d0, 0, f1)
call ddcpr (ss, f1, ic1)
call ddinfr (ss, t1, t2)

if (ic1 == 0) then
  write (ddldb, 2)
2 format ('*** DDZETAR: argument is 1')
  call ddabrt
elseif (ddsgn (t2) == 0) then

!   The argument is an integer value. Call mpzetaintr instead.

  call ddmdc (ss, d1, n1)
  iss = d1 * 2.d0**n1
  call ddzetaintr (iss, t1)
  goto 200
elseif (ddsgn (ss) < 0) then

!   If arg < 0, compute zeta(1-ss), and later apply Riemann's formula.

  call ddsub (f1, ss, tt)
else
  call ddeq (ss, tt)
endif

!  Check if argument is large enough that computing with definition is faster.

d1 = ddnbt * ddnw * log (2.d0) / log (2.d0 * ddnbt * ddnw / 3.d0)
call ddmdc (tt, d2, n2)
d2 = d2 * 2.d0 ** n2

if (d2 > d1) then

!   Evaluate the infinite series.

  call dddmc (1.d0, 0, t1)

  do i = 2, itrmax
    call dddmc (dble (i), 0, t4)
    call ddpower (t4, tt, t2)
    call dddiv (f1, t2, t3)
    call ddadd (t1, t3, t2)
    call ddeq (t2, t1)

    call ddabs (t3, tc1)
    call ddmul (eps, t1, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 200
  enddo

  write (ddldb, 3) 1, itrmax
3 format ('*** DDZETAR: iteration limit exceeded',2i10)
  call ddabrt
endif

n = dfrac * ddnw1
call dddmc (2.d0, 0, t1)
call ddnpwr (t1, n, tn)
call ddneg (tn, t1)
call dddmc (0.d0, 0, t2)
call dddmc (0.d0, 0, s)

sgn = 1.d0

do j = 0, 2 * n - 1
  call dddmc (dble (j + 1), 0, t4)
  call ddpower (t4, tt, t3)
  call dddiv (t1, t3, t4)
  call ddmuld (t4, sgn, t5)
  call ddadd (s, t5, t4)
  call ddeq (t4, s)
  sgn = - sgn

  if (j < n - 1) then
    call dddmc (0.d0, 0, t2)
  elseif (j == n - 1) then
    call dddmc (1.d0, 0, t2)
  else
    call ddmuld (t2, dble (2 * n - j), t3)
    call dddivd (t3, dble (j + 1 - n), t2)
  endif
  call ddadd (t1, t2, t3)
  call ddeq (t3, t1)
enddo

call ddsub (f1, tt, t3)
call dddmc (2.d0, 0, t2)
call ddpower (t2, t3, t4)
call ddsub (f1, t4, t2)
call ddmul (tn, t2, t3)
call dddiv (s, t3, t1)
call ddneg (t1, t2)
call ddeq (t2, t1)

!   If original argument was negative, apply Riemann's formula.

if (ddsgn (ss) < 0) then
  call ddgammar (tt, t3)
  call ddmul (t1, t3, t2)
  call ddmul (ddpicon, tt, t1)
  call ddmuld (t1, 0.5d0, t3)
  call ddcssnr (t3, t4, t5)
  call ddmul (t2, t4, t1)
  call ddmuld (ddpicon, 2.d0, t2)
  call ddpower (t2, tt, t3)
  call dddiv (t1, t3, t2)
  call ddmuld (t2, 2.d0, t1)
endif

200 continue

call ddeq (t1, zz)
return
end subroutine ddzetar

subroutine ddzetaintr (iss, zz)

!   This returns the zeta function of the integer argument ISS using an algorithm
!   due to Peter Borwein.

implicit none
integer, intent(in):: iss
real (ddknd), intent(out):: zz(1:2)
integer, parameter:: itrmax = 1000000
real (ddknd), parameter:: dfrac = 1.d0+dddpw, pi = 3.1415926535897932385d0
integer i, ic1, j, n, n1, itt
real (ddknd) d1, sgn
real (ddknd) f1(1:2), s(1:2), t1(1:2), t2(1:2), &
  t3(1:2), t4(1:2), t5(1:2), tn(1:2), &
  tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

ddnw1 = min (ddnw + 1, ddnwx)
call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw1*ddnbt, eps)

call dddmc (1.d0, 0, f1)

if (iss == 1) then
  write (ddldb, 2)
2 format ('*** DDZETAINTR: argument is 1')
  call ddabrt
elseif (iss == 0) then

!   Argument is zero -- result is -1/2.

  call dddmc (-0.5d0, 0, t1)
  goto 200
elseif (iss < 0) then

!   If argument is a negative even integer, the result is zero.

  if (mod (iss, 2) == 0) then
    call dddmc (0.d0, 0, t1)
    goto 200
  endif

!   Otherwise if arg < 0, compute zeta(1-is), and later apply Riemann's formula.

  itt = 1 - iss
else
  itt = iss
endif

!  Check if argument is large enough that computing with definition is faster.

d1 = ddnbt * ddnw * log (2.d0) / log (2.d0 * ddnbt * ddnw / 3.d0)

if (itt > d1) then

!   Evaluate the infinite series.

  call dddmc (1.d0, 0, t1)

  do i = 2, itrmax
    call dddmc (dble (i), 0, t4)
    call ddnpwr (t4, itt, t2)
    call dddiv (f1, t2, t3)
    call ddadd (t1, t3, t2)
    call ddeq (t2, t1)

    call ddabs (t3, tc1)
    call ddmul (eps, t1, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 200
  enddo

  write (ddldb, 3) 1, itrmax
3 format ('*** DDZETAINTR: iteration limit exceeded',2i10)
  call ddabrt
endif

n = dfrac * ddnw1
call dddmc (2.d0, 0, t1)
call ddnpwr (t1, n, tn)
call ddneg (tn, t1)
call dddmc (0.d0, 0, t2)
call dddmc (0.d0, 0, s)

sgn = 1.d0

do j = 0, 2 * n - 1
  call dddmc (dble (j + 1), 0, t4)
  call ddnpwr (t4, itt, t3)
  call dddiv (t1, t3, t4)
  call ddmuld (t4, sgn, t5)
  call ddadd (s, t5, t4)
  call ddeq (t4, s)
  sgn = - sgn

  if (j < n - 1) then
    call dddmc (0.d0, 0, t2)
  elseif (j == n - 1) then
    call dddmc (1.d0, 0, t2)
  else
    call ddmuld (t2, dble (2 * n - j), t3)
    call dddivd (t3, dble (j + 1 - n), t2)
  endif

  call ddadd (t1, t2, t3)
  call ddeq (t3, t1)
enddo

call dddmc (2.d0, 0, t2)
call ddnpwr (t2, 1 - itt, t4)
call ddsub (f1, t4, t2)
call ddmul (tn, t2, t3)
call dddiv (s, t3, t1)
call ddneg (t1, t2)
call ddeq (t2, t1)

!   If original argument was negative, apply Riemann's formula.

if (iss < 0) then
  call dddmc (1.d0, 0, t3)
  do i = 1, itt - 1
    call ddmuld (t3, dble (i), t4)
    call ddeq (t4, t3)
  enddo

  call ddmul (t1, t3, t2)
  call ddmuld (ddpicon, dble (itt), t1)
  call ddmuld (t1, 0.5d0, t3)
  call ddcssnr (t3, t4, t5)
  call ddmul (t2, t4, t1)
  call ddmuld (ddpicon, 2.d0, t2)
  call ddnpwr (t2, itt, t3)
  call dddiv (t1, t3, t2)
  call ddmuld (t2, 2.d0, t1)
endif

200 continue

call ddeq (t1, zz)
return
end subroutine ddzetaintr

subroutine ddzetabe (nb2, berne, s, z)

!  This evaluates the Riemann zeta function, using the combination of
!  the definition formula (for large s), and an Euler-Maclaurin scheme
!  (see formula 25.2.9 of the DLMF). The array berne contains precomputed
!  even Bernoulli numbers (see MPBERNER above). Its dimensions must be as
!  shown below. NB2 must be greater than 1.4 x precision in decimal digits.

implicit none
integer, intent(in):: nb2
real (ddknd), intent(in):: berne(1:2,nb2), s(1:2)
real (ddknd), intent(out):: z(1:2)
integer, parameter:: itrmax = 1000000
real (ddknd), parameter:: dber = 1.5d0, dfrac = 0.6d0, &
  pi = 3.1415926535897932385d0
integer i, i1, i2, ic1, k, n1, n2, nn
real (ddknd) d1, d2
real (ddknd) t0(1:2), t1(1:2), t2(1:2), t3(1:2), &
  t4(1:2), t5(1:2), t6(1:2), t7(1:2), t8(1:2), &
  t9(1:2), tt(1:2), f1(1:2), tc1(1:2), tc2(1:2), tc3(1:2), eps(1:2)

!  End of declaration
integer ddnw, ddnw1, ddnw2



ddnw = ddnwx

!   Check if berne array has been initialized.

call ddmdc (berne(1:2,1), d1, n1)
d1 = d1 * 2.d0 ** n1
if (ddwprecr (berne(1:2,1)) < ddnw .or. &
  abs (d1 - 1.d0 / 6.d0) > ddrdfz .or. nb2 < int (dber * dddpw * ddnw)) then
  write (ddldb, 3) int (dber * dddpw * ddnw)
3 format ('*** DDZETABE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries.')
  call ddabrt
endif

i = 0
k = 0
ddnw1 = min (ddnw + 1, ddnwx)

call dddmc (2.d0, 0, tc1)
call ddnpwr (tc1, -ddnw1*ddnbt, eps)

!   Check if argument is 1 -- undefined.

call dddmc (1.d0, 0, t0)
call ddcpr (s, t0, ic1)
if (ic1 == 0) then
  write (ddldb, 2)
2 format ('*** DDZETABE: argument is 1')
  call ddabrt
endif

call dddmc (1.d0, 0, f1)

!   Check if argument is zero. If so, result is - 1/2.

if (ddsgn (s) == 0) then
  call dddmc (-0.5d0, 0, t1)
  goto 200
endif

!   Check if argument is negative.

if (ddsgn (s) < 0) then

!   Check if argument is a negative even integer. If so, the result is zero.

  call ddmuld (s, 0.5d0, t1)
  call ddinfr (t1, t2, t3)
  if (ddsgn (t3) == 0) then
    call dddmc (0.d0, 0, t1)
    goto 200
  endif

!   Otherwise compute zeta(1-s), and later apply the reflection formula.

  call ddsub (f1, s, tt)
else
  call ddeq (s, tt)
endif

!  Check if argument is large enough that computing with definition is faster.

d1 = ddlogb * ddnw1 / log (32.d0 * ddnw1)
call ddmdc (tt, d2, n2)
d2 = d2 * 2.d0**n2
if (d2 > d1) then
  call dddmc (1.d0, 0, t1)

  do i = 2, itrmax
    call dddmc (dble (i), 0, t4)
    call ddpower (t4, tt, t2)
    call dddiv (f1, t2, t3)
    call ddadd (t1, t3, t2)
    call ddeq (t2, t1)

    call ddabs (t3, tc1)
    call ddmul (eps, t1, tc3)
    call ddabs (tc3, tc2)
    call ddcpr (tc1, tc2, ic1)
    if (ic1 <= 0) goto 200
  enddo

  write (ddldb, 4) 1, itrmax
4 format ('*** DDZETABE: iteration limit exceeded',2i10)
  call ddabrt
endif

call dddmc (1.d0, 0, t0)
nn = dfrac * dddpw * ddnw1

do k = 2, nn
  call dddmc (dble (k), 0, t2)
  call ddpower (t2, tt, t1)
  call dddiv (f1, t1, t2)
  call ddadd (t0, t2, t3)
  call ddeq (t3, t0)
enddo

call dddmc (dble (nn), 0, t2)
call ddsub (tt, f1, t3)
call ddmul (t1, t3, t4)
call dddiv (t2, t4, t3)
call ddadd (t0, t3, t2)
call dddmc (0.5d0, 0, t3)
call dddiv (t3, t1, t4)
call ddsub (t2, t4, t0)

call ddeq (tt, t3)
d1 = 12.d0 * dble (nn)
call ddmuld (t1, d1, t4)
call dddiv (t3, t4, t2)
call ddmuld (t1, dble (nn), t5)
call dddmc (dble (nn), 0, t6)
call ddmul (t6, t6, t9)

do k = 2, min (nb2, itrmax)
  call dddmc (dble (2 * k - 2), 0, t4)
  call ddadd (tt, t4, t6)
  call dddmc (dble (2 * k - 3), 0, t7)
  call ddadd (tt, t7, t8)
  call ddmul (t6, t8, t7)
  call ddmul (t3, t7, t4)
  call dddmc (dble (2 * k - 1), 0, t6)
  call dddmc (dble (2 * k - 2), 0, t7)
  call ddmul (t6, t7, t8)
  call dddiv (t4, t8, t3)
  call ddmul (t5, t9, t6)
  call ddeq (t6, t5)
  call ddmul (t3, berne(1:2,k), t4)
  call ddmuld (t5, dble (2 * k), t6)
  call dddiv (t4, t6, t7)
  call ddadd (t2, t7, t4)
  call ddeq (t4, t2)

  call ddabs (t7, tc1)
  call ddmul (eps, t2, tc3)
  call ddabs (tc3, tc2)
  call ddcpr (tc1, tc2, ic1)
  if (ic1 <= 0) goto 110
enddo

write (ddldb, 4) 2, min (nb2, itrmax)
call ddabrt

110 continue

call ddadd (t0, t2, t1)

!   If original argument was negative, apply the reflection formula.

if (ddsgn (s) < 0) then
  call ddgammar (tt, t3)
  call ddmul (t1, t3, t2)
  call ddmul (ddpicon, tt, t1)
  call ddmuld (t1, 0.5d0, t3)
  call ddcssnr (t3, t4, t5)
  call ddmul (t2, t4, t1)
  call ddmuld (ddpicon, 2.d0, t2)
  call ddpower (t2, tt, t3)
  call dddiv (t1, t3, t2)
  call ddmuld (t2, 2.d0, t1)
endif

200 continue

call ddeq (t1, z)

return
end subroutine ddzetabe

end module ddfune
