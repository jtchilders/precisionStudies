program tpslqm1dd

!   This performs the one-level multipair PSLQ algorithm on an algebraic input.
!   This version uses the DDFUN double-double precision package.

!   David H Bailey   7 Jan 2023

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2023 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!   This program demonstrates pslqm1, which performs the one-level multipair
!   PSLQ algorithm on an input vector. A variety of sample input vectors can
!   be generated as inputs to pslqm1, as given in the parameters below. The
!   pslqm1 routine is suitable for relations up to degree 20 or so; above this
!   level pslqm2 or pslqm3 should be used (which require higher precision).
!   This version uses the DDFUN double-double precision package.

!   For additional details, see:

!   David H. Bailey and David J. Broadhurst, "Parallel integer relation
!   detection: Techniques and applications," Mathematics of Computation,
!   vol. 70, no. 236 (Oct 2000), pg. 1719-1736, preprint available at
!   http://www.davidhbailey.com/dhbpapers/ppslq.pdf.

!   The pslqm1 routine is 100% THREAD SAFE -- all requisite parameters and
!   arrays are passed through subroutine arguments.

!   These parameters are set in the parameter statement below.

!     idb   Debug level (0 - 3); default = 2.
!     itr   Number of trials when kq = 1, 2 or 3; default = 1.
!     n     Integer relation vector length, = 1 + polynomial degree.
!           When kq = 0, set n = kr * ks + 1; this is the default.
!     kq    0: for the algebraic case [1, al, al^2, ... al^(n-1)], where
!             al = 3^(1/kr) - 2^(1/ks); this is the default.
!           1: for testing algebraic relations of a number read from a file.
!           2: for testing additive relations of numbers read from a file.
!           3: for testing multiplicative relations of numbers read from a file.
!           4: for custom input.
!     kr    Degree of root of 3 when kq = 0; default = 3.
!     ks    Degree of root of 2 when kq = 0; default = 3.
!     ndp   Full precision level in digits; default = 32.
!     nep   Log10 of full precision epsilon for detections; default = 3 - ndp.
!           ***Must not be smaller than the accuracy of input data.
!     rb    Log10 of max size (Euclidean norm) of acceptable relation; DP;
!           default = 20. Run will abort if this is exceeded.

use ddmodule
implicit none
integer, parameter:: idb = 2, itr = 1, kq = 0, kr = 3, ks = 3, &
  m1 = 40, m2 = 2, n = kr * ks + 1, ndp = 32, nep = 3 - ndp
real (ddknd), parameter:: rb = 20.d0
real (ddknd), external:: second
real (ddknd) d1, d2, rm, tm, tm0, tm1
integer i, i1, iq, j, j1, k
integer lnm(n)
character(1) chr1(m1)
character(64) nam(n), namx
type (dd_real) al, eps, t1, t2, r(n), x(n)

eps = ddreal (10.d0) ** nep
write (6, 1) itr, n, kq, kr, ks, rb, ndp, nep
1 format ('PSLQM1 Test Program'/ &
  'No. trials = ',i3,3x,'n =',i4,3x,'kq =',i2/ &
  'kr =',i2,3x,'ks =',i2,3x,'rb =',1p,d12.4/ &
  'Full precision level ndp =',i6,' digits'/ &
  'Full precision epsilon level nep = ',i6)

if (kq .eq. 1 .or. kq .eq. 2 .or. kq .eq. 3) then
  open (11, file = 'pslq.inp')
  rewind 11
endif

do k = 1, itr
  write (6, 2) k
2 format (/'Start trial', i3)

  if (kq .eq. 0) then

!   This code generates al = 3^(1/kr) - 2^(1/ks).  al is algebraic of degree
!   kr * ks.  Set n = kr * ks + 1 to recover the polynomial satisfied by al.

    al = ddnrt (ddreal (3.d0), kr) - ddnrt (ddreal (2.d0), ks)
  elseif (kq .eq. 1) then

!   Read an algebraic constant from a file.

    call ddread (11, al)
  elseif (kq .eq. 2) then

!   Read constants from a file for additive test.

    do i = 1, n
      call ddread (11, al)
      x(i) = al
      write (namx, '(''con'',i3.3)') i
      nam(i) = namx(1:6)
      lnm(i) = 6
    enddo
  elseif (kq .eq. 3) then

!   Read constants from a file for multiplicative test.

    do i = 1, n
      call ddread (11, al)
      x(i) = log (al)
      write (namx, '(''log(con'',i3.3,'')'')') i
      nam(i) = namx(1:11)
      lnm(i) = 11
    enddo
  elseif (kq .eq. 4) then

!   Produce X vector by a custom scheme.

  endif

!   If kq is 0 or 1, generate x = [1, al, al^2, ..., al^(n-1)].

  if (kq .eq. 0 .or. kq .eq. 1) then
    x(1) = ddreal (1.d0)
    nam(1) = '1'
    lnm(1) = 1

    do i = 2, n
      x(i) = al * x(i-1)
      write (namx, '(''al^'',i3)') i - 1
      nam(i) = namx(1:6)
      lnm(i) = 6
    enddo
  endif

!   Perform relation search.

  tm0 = second ()
  call pslqm1 (idb, n, rb, eps, x, iq, r)
  tm1 = second ()

!   Output relation, if one was found.

  if (iq == 1) then
    write (6, 3)
3   format (/'Recovered relation:')

    do i = 1, n
      call ddfform (r(i), m1, m2, chr1)
      write (6, '(80a1)') (chr1(i1), i1 = 1, m1)
    enddo
  endif

  write (6, 6) tm1 - tm0
6 format ('CPU Time =',f12.4)
enddo

if (iq == 1) then
  write (6, '(a)') 'TEST PASSED'
else
  write (6, '(a)') 'TEST FAILED'
endif
stop
end program tpslqm1dd

!------------------------------

!   The following code performs the one-level, multi-pair PSLQ algorithm.
!   David H. Bailey     7 Jan 2023

subroutine pslqm1 (idb, n, rb, eps, x, iq, r)

!   Arguments are as follows:
!     Name  Type    Description
!     idb   int     Debug flag (0-3); increasing idb produces more output.
!     n     int     Length of input vector x and output relation r.
!     rb    DP      Log10 of max size (Euclidean norm) of acceptable relation.
!                   NOTE: Run is aborted when this is exceeded.
!     eps   dd_real Tolerance for full precision relation detection. 
!     x     dd_real Input dd_real vector.
!     iq    int     Output flag: 0 (unsuccessful) or 1 (successful).
!     r     dd_real Output integer relation vector, if successful.

!   The following parameters are set in this routine:
!     ipi   int     Iteration print interval when idb >= 2; default = 25.
!     ipm   int     Iteration  check interval for iterations; default = 100.
!     itm   int     Maximum iteration count; default = 10^5.
!                   NOTE: Run is aborted when this is exceeded. If itm >= 10^7,
!                     change all "i7" in formats to "i8" or as high as needed.
!     nsq   int     Size of tables used in iterdp and itermp; default = 8.

use ddmodule
implicit none
integer, intent(in):: idb, n
real (ddknd), intent(in):: rb
type (dd_real), intent(in):: eps, x(n)
integer, intent(out):: iq
type (dd_real), intent(out):: r(n)
integer, parameter:: ipi = 25, ipm = 100, itm = 100000, nsq = 8
real (ddknd), external:: dplog10, second
type (dd_real), external:: bound
integer i, imq, it, izm, j, j1, n1, n2, n3, n4
real (ddknd) d1, d2, d3, d4, tm0, tm1, times(2)
type (dd_real) b(n,n), h(n,n), s(n), syq(n,nsq), y(n), rn, t1, t2, t3, t4

!   Initialize.

if (idb >= 2) write (6, 1) n
1 format ('PSLQM1 integer relation detection: n =',i5)
iq = 0
it = 0
imq = 0
rn = ddreal (0.d0)

do i = 1, 2
  times(i) = 0.d0
enddo

if (idb >= 2) write (6, 2) it
2 format ('Iteration',i7,3x,'Initialization')
tm0 = second ()
call initmp (idb, n, nsq, b, h, syq, x, y)
tm1 = second ()
times(1) = tm1 - tm0

!   Main iterations.

if (idb >= 2) write (6, 3) it
3 format ('Iteration',i7,3x,'Start iterations')

100 continue

!   Perform one iteration.

it = it + 1
if (idb == 3 .or. idb >= 2 .and. mod (it, ipi) == 0) write (6, 4) it
4 format ('Iteration',i7)
tm0 = second ()
call itermp (idb, it, n, nsq, eps, b, h, syq, y, imq, izm)
tm1 = second ()
times(2) = times(2) + (tm1 - tm0)

!   Test conditions on itermp output flag izm:
!   0: Update was uneventful; periodically output min, max, norm; continue.
!   1: A small value was found in update; go output relation.
!   2: Precision is exhausted; quit.

if (izm == 0) then

!   Periodically output min, max and norm.

  if (mod (it, ipm) == 0) then

!   Find min and max absolute value of y vector.

    call minmax (n, y, t1, t2)
    if (idb >= 2) then
      call decmd (t1, d1, n1)
      call decmd (t2, d2, n2)
      write (6, 5) it, d1, n1, d2, n2
5     format ('Iteration',i7,3x,'Min, max of y =',f11.6,'e',i6,f11.6,'e',i6)
    endif

!   Compute norm bound.

    t1 = bound (n, h)
    rn = max (rn, t1)
    if (idb >= 2) then
      call decmd (t1, d1, n1)
      call decmd (rn, d2, n2)
      write (6, 6) it, d1, n1, d2, n2
6     format ('Iteration',i7,3x,'Norm bound =',f11.6,'e',i5,4x,'Max. bound =', &
        f11.6,'e',i5)
    endif

!   Check if iteration limit or norm bound limit is exceeded; if so, quit.

    if (it > itm) then
      if (idb >= 1) write (6, 7) itm
7     format ('Iteration limit exceeded',i7)
      goto 120
    endif
    if (dplog10 (rn) > rb) then
      if (idb >= 1) write (6, 8) rb
8     format ('Norm bound limit exceeded.',1p,d15.6)
      goto 120
    endif
  endif
  goto 100
elseif (izm == 1) then
  goto 110
elseif (izm == 2) then
  goto 120
endif

110 continue

!   A relation has been detected.

tm0 = second ()
t1 = ddreal (1.d300)
t2 = ddreal (0.d0)
t3 = ddreal (0.d0)

!   Select the relation corresponding to the smallest y entry and compute norm.

do j = 1, n
  if (abs (y(j)) < t1) then
    j1 = j
    t1 = abs (y(j))
  endif
  t2 = max (t2, abs (y(j)))
enddo

do i = 1, n
  r(i) = b(j1,i)
  t3 = t3 + r(i) ** 2
enddo

t3 = sqrt (t3)
call decmd (t3, d3, n3)
d3 = d3 * 10.d0 ** n3

!   Output the final norm bound and other info.

if (idb >= 1) then
  t4 = bound (n, h)
  call decmd (t1, d1, n1)
  call decmd (t2, d2, n2)
  call decmd (t3, d3, n3)
  call decmd (t4, d4, n4)
  write (6, 9) it, d1, n1, d2, n2, d4, n4
9 format ('Iteration',i7,3x,'Relation detected'/ &
  'Min, max of y =',0p,f11.6,'e',i5,f11.6,'e',i5/'Max. bound =',f11.6,'e',i5)
  write (6, 10) j1, d3, n3, d1, n1
10 format ('Index of relation =',i4,3x,'Norm =',f11.6,'e',i5,3x, &
  'Residual =',f11.6,'e',i5)
endif

!   If run was successful, set iq = 1.

if (dplog10 (t3) <= rb) then
  iq = 1
else
  if (idb >= 2) write (6, 11)
11 format ('Relation is too large.')
endif

120 continue

!   Output CPU run times and return.

if (idb >= 2) write (6, 12) times
12 format ('CPU times:'/(5f12.2))

return
end subroutine pslqm1

!------------------------------

!   First-level subroutines.

subroutine minmax (n, y, y1, y2)

!   This returns min|y_k| and max|y_k| using full precision.
!   Input: n, y.
!   Output: y1, y2.

use ddmodule
implicit none
integer, intent(in):: n
type (dd_real), intent(in):: y(n)
type (dd_real), intent(out):: y1, y2
integer i
type (dd_real) t1, t2, t3

t1 = ddreal (1.d300)
t2 = ddreal (0.d0)

!   Find the min and max absolute value in the y vector.

do i = 1, n
  t3 = abs (y(i))
  t1 = min (t1, t3)
  t2 = max (t2, t3)
enddo

y1 = t1
y2 = t2
return
end subroutine minmax

subroutine initmp (idb, n, nsq, b, h, syq, x, y)

!   This initializes arrays at the beginning.
!   This is performed in full precision.
!   Input: idb, n, nsq, x.
!   Output: b, h, syq, y.

use ddmodule
implicit none
integer, intent(in):: idb, n, nsq
type (dd_real), intent(in):: x(n)
type (dd_real), intent(out):: b(n,n), h(n,n), syq(n,nsq), y(n)
integer i, i1, j
real (ddknd) d1
type (dd_real) s(n), t1, t2

if (idb >= 3) then
  write (6, 1)
1 format ('initmp: Input x vector:')
  call matoutmd (1, n, x)
endif

!   Set b to the identity matrix.

do j = 1, n
  do i = 1, n
    b(i,j) = ddreal (0.d0)
  enddo

  b(j,j) = ddreal (1.d0)
enddo

t1 = ddreal (0.d0)

!   Compute the s vector, the square root of the partial sum of squares of x,
!   and the y vector, which is the normalized x vector.

do i = n, 1, -1
  t1 = t1 + x(i) ** 2
  s(i) = sqrt (t1)
enddo

t1 = 1.d0 / s(1)

do i = 1, n
  y(i) = t1 * x(i)
  s(i) = t1 * s(i)
enddo

!   Compute the initial h matrix.

do j = 1, n - 1
  do i = 1, j - 1
    h(i,j) = ddreal (0.d0)
  enddo

  h(j,j) = s(j+1) / s(j)
  t1 = y(j) / (s(j) * s(j+1))

  do i = j + 1, n
    h(i,j) = - y(i) * t1
  enddo
enddo

!   Zero the syq array.

do j = 1, nsq
  do i = 1, n
    syq(i,j) = ddreal (0.d0)
  enddo
enddo

if (idb >= 3) then
  write (6, 2)
2 format ('initmp: Initial y vector:')
  call matoutmd (1, n, y)
  write (6, 3)
3 format ('initmp: Initial h matrix:')
  call matoutmd (n, n - 1, h)
endif

return
end subroutine initmp

subroutine itermp (idb, it, n, nsq, eps, b, h, syq, y, imq, izm)

!   This performs one iteration of the PSLQM algorithm using DD arithmetic.
!   Input: idb, it, n, nsq, eps, b, h, syq, imq.
!   Output: b, h, syq, y, imq, izm.

use ddmodule
implicit none
integer, intent(in):: idb, it, n, nsq
integer, intent(inout):: imq
type (dd_real), intent(in):: eps
type (dd_real), intent(inout):: b(n,n), h(n,n), syq(n,nsq)
type (dd_real), intent(out):: y(n)
integer, intent(out):: izm
integer, parameter:: ntl = 4
integer i, ii, ij, im, im1, j, j1, j2, k, mpr, mq, n1
real (ddknd) d1, d2
integer ip(n), ir(n), is(n)
type (dd_real) q(n), t(n,n)
type (dd_real) gam, t1, t2, t3, t4, teps

teps = 2.d0 ** ntl * eps
izm = 0
mpr = nint (0.4d0 * n)
gam = sqrt (ddreal (4.d0) / ddreal (3.d0))

!   Compute q vector = {gam^i * |h(i,i)|}, then sort in ascending order.

do i = 1, n - 1
  q(i) = gam ** i * abs (h(i,i))
enddo

call qsortmp (n - 1, q, ip)

!   Select up to mpr disjoint pairs of indices (m,m+1), where m is an index
!   from the list of the largest q(i).

do i = 1, n
  is(i) = 0
enddo

if (imq == 0) then
  mq = mpr
else
  mq = 1
  imq = 0
endif
ii = n

do i = 1, mq
100 continue
  ii = ii - 1
  if (ii == 0) then
    mq = i - 1
    goto 110
  endif
  j1 = ip(ii)
  j2 = j1 + 1
  if (is(j1) .ne. 0 .or. is(j2) .ne. 0) goto 100
  ir(i) = j1
  is(j1) = 1
  is(j2) = 1
enddo

110 continue

!   Exchange the pairs of entries of y, and rows of b and h.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  t1 = y(im)
  y(im) = y(im1)
  y(im1) = t1

  do i = 1, n
    t1 = b(im,i)
    b(im,i) = b(im1,i)
    b(im1,i) = t1
  enddo

  do i = 1, n - 1
    t1 = h(im,i)
    h(im,i) = h(im1,i)
    h(im1,i) = t1
  enddo
enddo

!   Eliminate the "corners" produced by the above permutation in h.

do j = 1, mq
  im = ir(j)
  im1 = im + 1
  if (im <= n - 2) then
    t1 = h(im,im)
    t2 = h(im,im1)
    t3 = sqrt (t1 ** 2 + t2 ** 2)
    t1 = t1 / t3
    t2 = t2 / t3

    do i = im, n
      t3 = h(i,im)
      t4 = h(i,im1)
      h(i,im) = t1 * t3 + t2 * t4
      h(i,im1) = - t2 * t3 + t1 * t4
    enddo
  endif
enddo

!   Perform reduction on h, using the diagonal scheme.  Multipliers are
!   saved in the t array.

do i = 2, n
  do j = 1, n - i + 1
    ij = i + j - 1

    do k = j + 1, ij - 1
      h(ij,j) = h(ij,j) - t(ij,k) * h(k,j)
    enddo

    t(ij,j) = anint (h(ij,j) / h(j,j))
    h(ij,j) = h(ij,j) - t(ij,j) * h(j,j)
  enddo
enddo

!   Update y, using the t array.  Find min absolute value of y.

t1 = abs (y(n))
j1 = n

do j = 1, n - 1
  do i = j + 1, n
    y(j) = y(j) + t(i,j) * y(i)
  enddo

  if (abs (y(j)) < t1) then
    j1 = j
    t1 = abs (y(j))
  endif
enddo

!   Update b, using the t array.

do k = 1, n
  do j = 1, n - 1
    do i = j + 1, n
      b(j,k) = b(j,k) + t(i,j) * b(i,k)
    enddo
  enddo
enddo

!  Find the largest entry of b in the same row as the smallest y.

t2 = ddreal (0.d0)

do i = 1, n
  t2 = max (t2, abs (b(j1,i)))
enddo

if (t1 <= t2 * teps) then
  if (idb >= 2) then
    call decmd (t1, d1, n1) 
    write (6, 1) it, d1, n1
1   format ('Iteration',i7,3x,'itermp: Small value in y =',f11.6,'e',i5)
  endif
  if (t1 <= t2 * eps) then
    izm = 1
  else
    if (idb >= 1) write (6, 2) it
2   format ('Iteration',i7,3x,'itermp: Precision exhausted.')
    izm = 2
  endif
endif

!   Compare the y vector with those of recent iterations.  If a duplicate is
!   found, then the next iteration must be performed with mq = 1.

do j = 1, nsq
  t1 = ddreal (0.d0)

  do i = 1, n
    t1 = max (t1, abs (y(i) - syq(i,j)))
  enddo

  if (t1 <= t2 * teps) then
    if (idb >= 2) write (6, 3) it, j
 3  format ('Iteration',i7,3x,'itermp: Duplicate found, j =',i6)
    imq = 1
    goto 120
  endif
enddo

!   Place the vector y in the table syq.

120 continue

k = 1 + mod (it, nsq)

do i = 1, n
  syq(i,k) = y(i)
enddo

if (idb >= 3) then
  write (6, 4)
4 format ('itermp: Updated y:')
!  call matoutmd (1, n, y)
  call matoutmp (1, n, y)
  write (6, 5)
5 format ('itermp: Updated b matrix:')
  call matoutmd (n, n, b)
  write (6, 6)
6 format ('itermp: Updated h matrix:')
  call matoutmd (n, n - 1, h)
endif

return
end subroutine itermp

!------------------------------

!   Second- and third-level subroutines.

type (dd_real) function bound (n, h)

!   This computes the norm bound using dd arithmetic.

use ddmodule
implicit none
integer, intent(in):: n
type (dd_real), intent(in):: h(n,n)
integer i
type (dd_real) t1, t2

t1 = ddreal (0.d0)

do i = 1, n - 1
  t1 = max (t1, abs (h(i,i)))
enddo

bound = 1.d0 / t1

return
end function bound

real (kind (0.d0)) function dplog10 (a)

!   For input dd value a, this routine returns a DP approximation to log10 (a).

use ddmodule
implicit none
type (dd_real), intent(in):: a
real (ddknd) da, t1

da = dble (a)
if (da == 0.d0) then
  dplog10 = -999999.d0
else
  dplog10 = log10 (abs (da))
endif

100 continue
return
end function dplog10

subroutine decmd (a, b, ib)

!   For input dd value a, this routine returns DP b and integer ib such that 
!   a = b * 10^ib, with 1 <= abs (b) < 10 for nonzero a.

use ddmodule
implicit none
type (dd_real), intent(in):: a
real (ddknd), intent(out):: b
integer, intent(out):: ib
integer ia
real (ddknd) da, t1


da = dble (a)
if (da .ne. 0.d0) then
  t1 = log10 (abs (da))
  ib = t1
  if (t1 < 0.d0) ib = ib - 1
  b = sign (10.d0 ** (t1 - ib), da)
else
  b = 0.d0
  ib = 0
endif

return
end subroutine decmd

subroutine matoutmd (n1, n2, a)

!   This outputs the dd matrix a as a DP matrix.

use ddmodule
implicit none
integer, intent(in):: n1, n2
type (dd_real), intent(in):: a(n1,n2)
integer i, j, ix(n2)
real (ddknd) dx(n2)

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call decmd (a(i,j), dx(j), ix(j))
  enddo

  write (6, 2) (dx(j), ix(j), j = 1, n2)
2 format (4(f13.8,'e',i5))
enddo

return
end subroutine matoutmd

subroutine matoutmp (n1, n2, a)

!   This outputs the dd matrix a.  It may be used in place of calls to matoutmd
!   in the code above if greater accuracy is desired in debug output.

use ddmodule
implicit none
integer, intent(in):: n1, n2
type (dd_real), intent(in):: a(n1,n2)
integer, parameter:: m1 = 45, m2 = 35
integer i, j

do i = 1, n1
  write (6, 1) i
1 format ('Row',i3)

  do j = 1, n2
    call ddwrite (6, m1, m2, a(i,j))
  enddo
enddo

return
end subroutine matoutmp

subroutine qsortmp (n, a, ip)

!   This routine sorts the entries of the N-long dd vector A into ascending
!   order using the quicksort algorithm.  The permutation vector that would
!   sort the vector is returned in IP.
!   Input: n, a.
!   Output: ip.

use ddmodule
implicit none
integer, intent(in):: n
type (dd_real), intent(in):: a(n)
integer, intent(out):: ip(n)

integer i, iq, it, j, jq, jz, k, l
integer ik(50), jk(50)
type (dd_real) s0, s1, s2

do i = 1, n
  ip(i) = i
enddo

if (n == 1) return

k = 1
ik(1) = 1
jk(1) = n

130 continue

i = ik(k)
j = jk(k)
iq = i
jq = j
it = (i + j + 1) / 2
l = ip(j)
ip(j) = ip(it)
ip(it) = l
s0 = a(ip(j))
j = j - 1

140 continue

do l = i, j
  if (s0 < a(ip(l))) goto 160
enddo

i = j
goto 190

160 continue

i = l

do l = j, i, -1
  if (s0 > a(ip(l))) goto 180
enddo

j = i
goto 190

180 continue

j = l
if (i >= j) goto 190
l = ip(i)
ip(i) = ip(j)
ip(j) = l
goto 140

190 continue

if (s0 >= a(ip(i))) goto 200
l = ip(jq)
ip(jq) = ip(i)
ip(i) = l

200 continue

k = k - 1
jz = 0
if (j == iq) goto 210
k = k + 1
jk(k) = j
jz = 1

210 continue

i = i + 1
if (i == jq) goto 220
k = k + 1
ik(k) = i
jk(k) = jq
if (jz == 0) goto 220
if (j - iq >= jq - i) goto 220
ik(k-1) = i
jk(k-1) = jq
ik(k) = iq
jk(k) = j

220 continue

if (k > 0) goto 130

return
end subroutine qsortmp
