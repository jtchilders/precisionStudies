
!  DDFUN: A double-double package with special functions

!  High-level language interface module (DDMODULE).

!  Revision date:  16 Mar 2023

!  AUTHOR:
!     David H. Bailey
!     Lawrence Berkeley National Lab (retired) and University of California, Davis
!     Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2023 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF PACKAGE:
!    This package permits one to perform floating-point computations (real and complex)
!    to double-double precision (approximately 30 digits), by making only relatively
!    minor changes to existing Fortran programs. It is intended as a substitute for
!    IEEE 128-bit quad precision on systems where quad precision is not available. All
!    basic arithmetic operations and transcendental functions are supported. The
!    package should run correctly on any Unix-based system supporting a Fortran-2008
!    compiler and IEEE 64-bit floating-point arithmetic. Note however that results are
!    NOT guaranteed to the last bit.

!    In addition to fast execution times, one key feature of this package is a
!    100% THREAD-SAFE design, which means that user-level applications can be
!    easily converted for parallel execution, say using a threaded parallel
!    environment such as OpenMP. There are no global shared variables (except
!    static compile-time data), and no initialization is necessary.

!  DOCUMENTATION:
!    See the README-ddfun.txt file in the main DDFUN directory.

!  DESCRIPTION OF THIS MODULE (DDMODULE):
!    This module contains all high-level Fortran language interfaces.

module ddmodule
use ddfuna
use ddfune
implicit none
type dd_real
  sequence
  real (ddknd) ddr(2)
end type
type dd_complex
  sequence
  real (ddknd) ddc(4)
end type

private &
  dd_eqqq, dd_eqqz, dd_eqdq, dd_eqqd, dd_eqxq, dd_eqqx, &
  dd_addqq, dd_addqz, dd_adddq, dd_addqd, dd_addxq, dd_addqx, &
  dd_subqq, dd_subqz, dd_subdq, dd_subqd, dd_subxq, dd_subqx, dd_negq, &
  dd_mulqq, dd_mulqz, dd_muldq, dd_mulqd, dd_mulxq, dd_mulqx, &
  dd_divqq, dd_divqz, dd_divdq, dd_divqd, dd_divxq, dd_divqx, &
  dd_expqq, dd_expqi, dd_expdq, dd_expqd, &
  dd_eqtqq, dd_eqtqz, dd_eqtdq, dd_eqtqd, dd_eqtxq, dd_eqtqx, &
  dd_netqq, dd_netqz, dd_netdq, dd_netqd, dd_netxq, dd_netqx, &
  dd_letqq, dd_letdq, dd_letqd, dd_getqq, dd_getdq, dd_getqd, &
  dd_lttqq, dd_lttdq, dd_lttqd, dd_gttqq, dd_gttdq, dd_gttqd

private &
  dd_eqzq, dd_eqzz, dd_eqdz, dd_eqzd, dd_eqxz, dd_eqzx, &
  dd_addzq, dd_addzz, dd_adddz, dd_addzd, dd_addxz, dd_addzx, &
  dd_subzq, dd_subzz, dd_subdz, dd_subzd, dd_subxz, dd_subzx, dd_negz, &
  dd_mulzq, dd_mulzz, dd_muldz, dd_mulzd, dd_mulxz, dd_mulzx, &
  dd_divzq, dd_divzz, dd_divdz, dd_divzd, dd_divxz, dd_divzx, &
  dd_expzi, dd_expzz, dd_exprz, dd_expzr, &
  dd_eqtzq, dd_eqtzz, dd_eqtdz, dd_eqtzd, dd_eqtxz, dd_eqtzx, &
  dd_netzq, dd_netzz, dd_netdz, dd_netzd, dd_netxz, dd_netzx

private &
  dd_absq, dd_absz, dd_acos, dd_acosh, dd_agm, dd_aimag, dd_aint, &
  dd_anint, dd_asin, dd_asinh, dd_atan, dd_atan2, dd_atanh, dd_atoq, &
  dd_bessel_i, dd_bessel_in, dd_bessel_j, dd_bessel_j0, dd_bessel_j1, &
  dd_bessel_jn, dd_bessel_k, dd_bessel_kn, dd_bessel_y, dd_bessel_yn, &
  dd_conjg, dd_cos, dd_cosz, dd_cosh, dd_cssh, dd_cssn, dd_digamma_be, &
  dd_dtoq, dd_dtoz, dd_ddtoz, dd_eform, dd_egamma, dd_erf, &
  dd_erfc, dd_exp, dd_expint, dd_expz, dd_fform, dd_gamma, &
  dd_hurwitz_zetan, dd_hurwitz_zetan_be, dd_hypergeom_pfq, dd_hypot, &
  dd_incgamma, dd_inpq, dd_inpz, dd_itoq, dd_log, dd_logz, dd_log2, &
  dd_log10, dd_maxq, dd_maxq3, dd_minq, dd_minq3, dd_modq, dd_nrt, &
  dd_outq, dd_outz, dd_pi, dd_poly, dd_polygamma, dd_polygamma_be, &
  dd_polylog_ini, dd_polylog_neg, dd_polylog_pos, &
  dd_qtod, dd_qqtoz, dd_qtoq, dd_qtox, dd_qtoz, dd_signq, &
  dd_sin, dd_sinz, dd_sinh, dd_sqrtq, dd_sqrtz, dd_struve_hn, dd_tan, &
  dd_tanh, dd_xtoq, dd_xtoz, dd_zeta, &
  dd_zeta_be, dd_zeta_int, dd_ztod, dd_ztoq, dd_ztox, dd_ztoz

!   These are routine names in module DDFUNA (ddfuna.f90):

private &
  ddabrt, ddabs, ddadd, ddacosh, ddang, ddasinh, ddatanh, ddcadd, &
  ddcdiv, ddceq, ddcmul, ddcpr, ddcpwr, ddcsqrt, ddcsshr, ddcssnr, &
  ddcsub, dddiv, dddivd, dddmc, dddpddc, dddddpc, ddeq, dddigin, &
  dddigout, ddeformat, ddfformat, ddexp, ddinfr, ddinp, ddinpc, ddlog, &
  ddlog2c, ddmdc, ddmul, ddmuld, ddmuldd, ddneg, ddnint, ddpower, &
  ddnpwr, ddnrtf, ddout, ddpic, ddpoly, ddqqc, ddsqrt, ddsub, ddxzc, &
  ddmzc

!   These are routine names in module DDFUNE (ddfune.f90):

private &
  ddwprecr, ddspacer, ddberner, ddpolyadd, ddpolysub, ddpolymul, &
  ddbesselir, ddbesselinr, ddbesseljr, ddbesseljnr, ddbesselkr, &
  ddbesselknr,  ddbesselyr, ddbesselynr, dddigammabe, dderfr, &
  dderfcr, ddexpint, ddgammar, ddhurwitzzetan, ddhurwitzzetanbe, &
  ddhypergeompfq, ddincgammar, ddpolygamma, ddpolygammabe, &
  ddpolylogini, ddpolylogneg, ddpolylogpos, ddstruvehn, ddzetar, &
  ddzetaintr, ddzetabe

  

!  Operator extension interface blocks.

interface assignment (=)
  module procedure dd_eqqq
  module procedure dd_eqqz
  module procedure dd_eqdq
  module procedure dd_eqqd
  module procedure dd_eqxq
  module procedure dd_eqqx
  module procedure dd_eqqa
  module procedure dd_eqzq
  module procedure dd_eqzz
  module procedure dd_eqdz
  module procedure dd_eqzd
  module procedure dd_eqxz
  module procedure dd_eqzx
end interface

interface operator (+)
  module procedure dd_addqq
  module procedure dd_addqz
  module procedure dd_adddq
  module procedure dd_addqd
  module procedure dd_addxq
  module procedure dd_addqx
  module procedure dd_addzq
  module procedure dd_addzz
  module procedure dd_adddz
  module procedure dd_addzd
  module procedure dd_addxz
  module procedure dd_addzx
end interface

interface operator (-)
  module procedure dd_subqq
  module procedure dd_subqz
  module procedure dd_subdq
  module procedure dd_subqd
  module procedure dd_subxq
  module procedure dd_subqx
  module procedure dd_subzq
  module procedure dd_subzz
  module procedure dd_subdz
  module procedure dd_subzd
  module procedure dd_subxz
  module procedure dd_subzx

  module procedure dd_negq
  module procedure dd_negz
end interface

interface operator (*)
  module procedure dd_mulqq
  module procedure dd_mulqz
  module procedure dd_muldq
  module procedure dd_mulqd
  module procedure dd_mulxq
  module procedure dd_mulqx
  module procedure dd_mulzq
  module procedure dd_mulzz
  module procedure dd_muldz
  module procedure dd_mulzd
  module procedure dd_mulxz
  module procedure dd_mulzx
end interface

interface operator (/)
  module procedure dd_divqq
  module procedure dd_divqz
  module procedure dd_divdq
  module procedure dd_divqd
  module procedure dd_divxq
  module procedure dd_divqx
  module procedure dd_divzq
  module procedure dd_divzz
  module procedure dd_divdz
  module procedure dd_divzd
  module procedure dd_divxz
  module procedure dd_divzx
end interface

interface operator (**)
  module procedure dd_expqq
  module procedure dd_expqi
  module procedure dd_expdq
  module procedure dd_expqd
  module procedure dd_expzi
  module procedure dd_expzz
  module procedure dd_exprz
  module procedure dd_expzr
end interface

interface operator (==)
  module procedure dd_eqtqq
  module procedure dd_eqtqz
  module procedure dd_eqtdq
  module procedure dd_eqtqd
  module procedure dd_eqtxq
  module procedure dd_eqtqx
  module procedure dd_eqtzq
  module procedure dd_eqtzz
  module procedure dd_eqtdz
  module procedure dd_eqtzd
  module procedure dd_eqtxz
  module procedure dd_eqtzx
end interface

interface operator (/=)
  module procedure dd_netqq
  module procedure dd_netqz
  module procedure dd_netdq
  module procedure dd_netqd
  module procedure dd_netxq
  module procedure dd_netqx
  module procedure dd_netzq
  module procedure dd_netzz
  module procedure dd_netdz
  module procedure dd_netzd
  module procedure dd_netxz
  module procedure dd_netzx
end interface

interface operator (<=)
  module procedure dd_letqq
  module procedure dd_letdq
  module procedure dd_letqd
end interface

interface operator (>=)
  module procedure dd_getqq
  module procedure dd_getdq
  module procedure dd_getqd
end interface

interface operator (<)
  module procedure dd_lttqq
  module procedure dd_lttdq
  module procedure dd_lttqd
end interface

interface operator (>)
  module procedure dd_gttqq
  module procedure dd_gttdq
  module procedure dd_gttqd
end interface

interface abs
  module procedure dd_absq
  module procedure dd_absz
end interface

interface acos
  module procedure dd_acos
end interface

interface acosh
  module procedure dd_acosh
end interface

interface agm
  module procedure dd_agm
end interface

interface aimag
  module procedure dd_aimag
end interface

interface aint
  module procedure dd_aint
end interface

interface anint
  module procedure dd_anint
end interface

interface asin
  module procedure dd_asin
end interface

interface asinh
  module procedure dd_asinh
end interface

interface atan
  module procedure dd_atan
end interface

interface atan2
  module procedure dd_atan2
end interface

interface atanh
  module procedure dd_atanh
end interface

interface ddberne
  module procedure dd_berne
end interface

interface bessel_i
  module procedure dd_bessel_i
end interface

interface bessel_in
  module procedure dd_bessel_in
end interface

interface bessel_j
  module procedure dd_bessel_j
end interface

interface bessel_jn
  module procedure dd_bessel_jn
end interface

interface bessel_j0
  module procedure dd_bessel_j0
end interface

interface bessel_j1
  module procedure dd_bessel_j1
end interface

interface bessel_k
  module procedure dd_bessel_k
end interface

interface bessel_kn
  module procedure dd_bessel_kn
end interface

interface bessel_y
  module procedure dd_bessel_y
end interface

interface bessel_yn
  module procedure dd_bessel_yn
end interface

interface bessel_y0
  module procedure dd_bessel_y0
end interface

interface bessel_y1
  module procedure dd_bessel_y1
end interface

interface conjg
  module procedure dd_conjg
end interface

interface cos
  module procedure dd_cos
  module procedure dd_cosz
end interface

interface cosh
  module procedure dd_cosh
end interface

interface dble
  module procedure dd_qtod
  module procedure dd_ztod
end interface

interface dcmplx
  module procedure dd_qtox
  module procedure dd_ztox
end interface

interface ddcmplx
  module procedure dd_ztoz
  module procedure dd_qtoz
  module procedure dd_dtoz
  module procedure dd_xtoz
  module procedure dd_qqtoz
  module procedure dd_ddtoz
end interface

interface ddcssh
  module procedure dd_cssh
end interface

interface ddcssn
  module procedure dd_cssn
end interface

interface ddeform
  module procedure dd_eform
end interface

interface ddfform
  module procedure dd_fform
end interface

interface ddnrt
  module procedure dd_nrt
end interface

interface ddlog2
  module procedure dd_log2
end interface

interface ddpi
  module procedure dd_pi
end interface

interface ddpoly
  module procedure dd_poly
end interface

interface ddread
  module procedure dd_inpq
  module procedure dd_inpz
end interface

interface ddreal
  module procedure dd_qtoq
  module procedure dd_ztoq
  module procedure dd_dtoq
  module procedure dd_xtoq
  module procedure dd_atoq
  module procedure dd_itoq
end interface

interface ddwrite
  module procedure dd_outq
  module procedure dd_outz
end interface

interface digamma_be
  module procedure dd_digamma_be
end interface

interface erf
  module procedure dd_erf
end interface

interface erfc
  module procedure dd_erfc
end interface

interface exp
  module procedure dd_exp
  module procedure dd_expz
end interface

interface expint
  module procedure dd_expint
end interface

interface gamma
  module procedure dd_gamma
end interface

interface hurwitz_zetan
  module procedure dd_hurwitz_zetan
end interface

interface hurwitz_zetan_be
  module procedure dd_hurwitz_zetan_be
end interface

interface hypergeom_pfq
  module procedure dd_hypergeom_pfq
end interface

interface hypot
  module procedure dd_hypot
end interface

interface incgamma
  module procedure dd_incgamma
end interface

interface log
  module procedure dd_log
  module procedure dd_logz
end interface

interface log10
  module procedure dd_log10
end interface

interface max
  module procedure dd_maxq
  module procedure dd_maxq3
end interface

interface min
  module procedure dd_minq
  module procedure dd_minq3
end interface

interface mod
  module procedure dd_modq
end interface

interface polygamma
  module procedure dd_polygamma
end interface

interface polygamma_be
  module procedure dd_polygamma_be
end interface

interface polylog_ini
  module procedure dd_polylog_ini
end interface

interface polylog_neg
  module procedure dd_polylog_neg
end interface

interface polylog_pos
  module procedure dd_polylog_pos
end interface

interface sign
  module procedure dd_signq
end interface

interface sin
  module procedure dd_sin
  module procedure dd_sinz
end interface

interface sinh
  module procedure dd_sinh
end interface

interface sqrt
  module procedure dd_sqrtq
  module procedure dd_sqrtz
end interface

interface struve_hn
  module procedure dd_struve_hn
end interface

interface tan
  module procedure dd_tan
end interface

interface tanh
  module procedure dd_tanh
end interface

interface zeta
  module procedure dd_zeta
end interface

interface zeta_be
  module procedure dd_zeta_be
end interface

interface zeta_int
  module procedure dd_zeta_int
end interface

contains

!  Assignment routines.

  subroutine dd_eqqq (qa, qb)
    implicit none
    type (dd_real), intent (out):: qa
    type (dd_real), intent (in):: qb
    qa%ddr(1) = qb%ddr(1)
    qa%ddr(2) = qb%ddr(2)
    return
  end subroutine

  subroutine dd_eqqz (qa, zb)
    implicit none
    type (dd_real), intent (out):: qa
    type (dd_complex), intent (in):: zb
    qa%ddr(1) = zb%ddc(1)
    qa%ddr(2) = zb%ddc(2)
    return
  end subroutine

  subroutine dd_eqdq (da, qb)
    implicit none
    real (ddknd), intent (out):: da
    type (dd_real), intent (in):: qb
    da = qb%ddr(1)
    return
  end subroutine

  subroutine dd_eqqd (qa, db)
    implicit none
    type (dd_real), intent (out):: qa
    real (ddknd), intent (in):: db
    qa%ddr(1) = db
    qa%ddr(2) = 0.d0
    return
  end subroutine

  subroutine dd_eqxq (xa, qb)
    implicit none
    complex (ddknd), intent (out):: xa
    type (dd_real), intent (in):: qb
    xa = qb%ddr(1)
    return
  end subroutine

  subroutine dd_eqqx (qa, xb)
    implicit none
    type (dd_real), intent (out):: qa
    complex (ddknd), intent (in):: xb
    qa%ddr(1) = xb
    qa%ddr(2) = 0.d0
    return
  end subroutine

  subroutine dd_eqqa (qa, ab)
    implicit none
    type (dd_real), intent (out):: qa
    character(*), intent (in):: ab
    character(120) cht
    cht = ab
    call ddinpc (cht, qa%ddr)
    return
  end subroutine

  subroutine dd_eqzq (za, qb)
    implicit none
    type (dd_complex), intent (out):: za
    type (dd_real), intent (in):: qb
    call ddmzc (qb%ddr, za%ddc)
    return
  end subroutine

  subroutine dd_eqzz (za, zb)
    implicit none
    type (dd_complex), intent (out):: za
    type (dd_complex), intent (in):: zb
    call ddceq (zb%ddc, za%ddc)
    return
  end subroutine

  subroutine dd_eqdz (da, zb)
    implicit none
    real (ddknd), intent (out):: da
    type (dd_complex), intent (in):: zb
    da = zb%ddc(1)
    return
  end subroutine

  subroutine dd_eqzd (za, db)
    implicit none
    type (dd_complex), intent (out):: za
    real (ddknd), intent (in):: db
    complex (ddknd) xb
    xb = db
    call ddxzc (xb, za%ddc)
    return
  end subroutine

  subroutine dd_eqxz (xa, zb)
    implicit none
    complex (ddknd), intent (out):: xa
    type (dd_complex), intent (in):: zb
    real (ddknd) db, dc
    db = zb%ddc(1)
    dc = zb%ddc(3)
    xa = cmplx (db, dc, ddknd)
    return
  end subroutine

  subroutine dd_eqzx (za, xb)
    implicit none
    type (dd_complex), intent (out):: za
    complex (ddknd), intent (in):: xb
    call ddxzc (xb, za%ddc)
    return
  end subroutine

!  Addition routines.

  function dd_addqq (qa, qb)
    implicit none
    type (dd_real) dd_addqq
    type (dd_real), intent (in):: qa, qb
    call ddadd (qa%ddr, qb%ddr, dd_addqq%ddr)
    return
  end function

  function dd_addqz (qa, zb)
    implicit none
    type (dd_complex) dd_addqz
    type (dd_real), intent (in):: qa
    type (dd_complex), intent (in):: zb
    real (ddknd) qt1(4)
    call ddmzc (qa%ddr, qt1)
    call ddcadd (qt1, zb%ddc, dd_addqz%ddc)
    return
  end function

  function dd_adddq (da, qb)
    implicit none
    type (dd_real):: dd_adddq
    real (ddknd), intent (in):: da
    type (dd_real), intent (in):: qb
    real (ddknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.d0
    call ddadd (qt1, qb%ddr, dd_adddq%ddr)
    return
  end function

  function dd_addqd (qa, db)
    implicit none
    type (dd_real):: dd_addqd
    type (dd_real), intent (in):: qa
    real (ddknd), intent (in):: db
    real (ddknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.d0
    call ddadd (qa%ddr, qt1, dd_addqd%ddr)
    return
  end function

  function dd_addxq (xa, qb)
    implicit none
    type (dd_complex):: dd_addxq
    complex (ddknd), intent (in):: xa
    type (dd_real), intent (in):: qb
    real (ddknd) qt1(4), qt2(4)
    call ddxzc (xa, qt1)
    call ddmzc (qb%ddr, qt2)
    call ddcadd (qt1, qt2, dd_addxq%ddc)
    return
  end function

  function dd_addqx (qa, xb)
    implicit none
    type (dd_complex):: dd_addqx
    type (dd_real), intent (in):: qa
    complex (ddknd), intent (in):: xb
    real (ddknd) qt1(4), qt2(4)
    call ddmzc (qa%ddr, qt1)
    call ddxzc (xb, qt2)
    call ddcadd (qt1, qt2, dd_addqx%ddc)
    return
  end function

  function dd_addzq (za, qb)
    implicit none
    type (dd_complex):: dd_addzq
    type (dd_complex), intent (in):: za
    type (dd_real), intent (in)::qb
    real (ddknd) qt1(4)
    call ddmzc (qb%ddr, qt1)
    call ddcadd (za%ddc, qt1, dd_addzq%ddc)
    return
  end function

  function dd_addzz (za, zb)
    implicit none
    type (dd_complex):: dd_addzz
    type (dd_complex), intent (in):: za, zb
    call ddcadd (za%ddc, zb%ddc, dd_addzz%ddc)
    return
  end function

  function dd_adddz (da, zb)
    implicit none
    type (dd_complex):: dd_adddz
    real (ddknd), intent (in):: da
    type (dd_complex), intent (in):: zb
    real (ddknd) qt1(4)
    complex (ddknd) xa
    xa = da
    call ddxzc (xa, qt1)
    call ddcadd (qt1, zb%ddc, dd_adddz%ddc)
    return
  end function

  function dd_addzd (za, db)
    implicit none
    type (dd_complex):: dd_addzd
    type (dd_complex), intent (in):: za
    real (ddknd), intent (in):: db
    real (ddknd) qt1(4)
    complex (ddknd) xb
    xb = db
    call ddxzc (xb, qt1)
    call ddcadd (za%ddc, qt1, dd_addzd%ddc)
    return
  end function

  function dd_addxz (xa, zb)
    implicit none
    type (dd_complex):: dd_addxz
    complex (ddknd), intent (in):: xa
    type (dd_complex), intent (in):: zb
    real (ddknd) qt1(4)
    call ddxzc (xa, qt1)
    call ddcadd (qt1, zb%ddc, dd_addxz%ddc)
    return
  end function

  function dd_addzx (za, xb)
    implicit none
    type (dd_complex):: dd_addzx
    type (dd_complex), intent (in):: za
    complex (ddknd), intent (in):: xb
    real (ddknd) qt1(4)
    call ddxzc (xb, qt1)
    call ddcadd (za%ddc, qt1, dd_addzx%ddc)
    return
  end function

!  Subtraction routines.

  function dd_subqq (qa, qb)
    implicit none
    type (dd_real):: dd_subqq
    type (dd_real), intent (in):: qa, qb
    call ddsub (qa%ddr, qb%ddr, dd_subqq%ddr)
    return
  end function

  function dd_subqz (qa, zb)
    implicit none
    type (dd_complex):: dd_subqz
    type (dd_real), intent (in):: qa
    type (dd_complex), intent (in):: zb
    real (ddknd) qt1(4)
    call ddmzc (qa%ddr, qt1)
    call ddcsub (qt1, zb%ddc, dd_subqz%ddc)
    return
  end function

  function dd_subdq (da, qb)
    implicit none
    type (dd_real):: dd_subdq
    real (ddknd), intent (in):: da
    type (dd_real), intent (in):: qb
    real (ddknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.d0
    call ddsub (qt1, qb%ddr, dd_subdq%ddr)
    return
  end function

  function dd_subqd (qa, db)
    implicit none
    type (dd_real):: dd_subqd
    type (dd_real), intent (in):: qa
    real (ddknd), intent (in):: db
    real (ddknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.d0
    call ddsub (qa%ddr, qt1, dd_subqd%ddr)
    return
  end function

  function dd_subxq (xa, qb)
    implicit none
    type (dd_complex):: dd_subxq
    complex (ddknd), intent (in):: xa
    type (dd_real), intent (in):: qb
    real (ddknd) qt1(4), qt2(4)
    call ddxzc (xa, qt1)
    call ddmzc (qb%ddr, qt2)
    call ddcsub (qt1, qt2, dd_subxq%ddc)
    return
  end function

  function dd_subqx (qa, xb)
    implicit none
    type (dd_complex):: dd_subqx
    type (dd_real), intent (in):: qa
    complex (ddknd), intent (in):: xb
    real (ddknd) qt1(4), qt2(4)
    call ddmzc (qa%ddr, qt1)
    call ddxzc (xb, qt2)
    call ddcsub (qt1, qt2, dd_subqx%ddc)
    return
  end function

  function dd_subzq (za, qb)
    implicit none
    type (dd_complex):: dd_subzq
    type (dd_complex), intent (in):: za
    type (dd_real), intent (in):: qb
    real (ddknd) qt1(4)
    call ddmzc (qb%ddr, qt1)
    call ddcsub (za%ddc, qt1, dd_subzq%ddc)
    return
  end function

  function dd_subzz (za, zb)
    implicit none
    type (dd_complex):: dd_subzz
    type (dd_complex), intent (in):: za, zb
    call ddcsub (za%ddc, zb%ddc, dd_subzz%ddc)
    return
  end function

  function dd_subdz (da, zb)
    implicit none
    type (dd_complex):: dd_subdz
    real (ddknd), intent (in):: da
    type (dd_complex), intent (in):: zb
    real (ddknd) qt1(4)
    complex (ddknd) xa
    xa = da
    call ddxzc (xa, qt1)
    call ddcsub (qt1, zb%ddc, dd_subdz%ddc)
    return
  end function

  function dd_subzd (za, db)
    implicit none
    type (dd_complex):: dd_subzd
    type (dd_complex), intent (in):: za
    real (ddknd), intent (in):: db
    real (ddknd) qt1(4)
    complex (ddknd) xb
    xb = db
    call ddxzc (xb, qt1)
    call ddcsub (za%ddc, qt1, dd_subzd%ddc)
    return
  end function

  function dd_subxz (xa, zb)
    implicit none
    type (dd_complex):: dd_subxz
    complex (ddknd), intent (in):: xa
    type (dd_complex), intent (in):: zb
    real (ddknd) qt1(4)
    call ddxzc (xa, qt1)
    call ddcsub (qt1, zb%ddc, dd_subxz%ddc)
    return
  end function

  function dd_subzx (za, xb)
    implicit none
    type (dd_complex):: dd_subzx
    type (dd_complex), intent (in):: za
    complex (ddknd), intent (in):: xb
    real (ddknd) qt1(4)
    call ddxzc (xb, qt1)
    call ddcsub (za%ddc, qt1, dd_subzx%ddc)
    return
  end function
  
!  ddr negation routines.

  function dd_negq (qa)
    implicit none
    type (dd_real):: dd_negq
    type (dd_real), intent (in):: qa
    call ddeq (qa%ddr, dd_negq%ddr)
    dd_negq%ddr(1) = - qa%ddr(1)
    dd_negq%ddr(2) = - qa%ddr(2)
    return
  end function

  function dd_negz (za)
    implicit none
    type (dd_complex):: dd_negz
    type (dd_complex), intent (in):: za
    call ddceq (za%ddc, dd_negz%ddc)
    dd_negz%ddc(1) = - za%ddc(1)
    dd_negz%ddc(2) = - za%ddc(2)
    dd_negz%ddc(3) = - za%ddc(3)
    dd_negz%ddc(4) = - za%ddc(4)
    return
  end function

!  ddr multiply routines.

  function dd_mulqq (qa, qb)
    implicit none
    type (dd_real):: dd_mulqq
    type (dd_real), intent (in):: qa, qb
    call ddmul (qa%ddr, qb%ddr, dd_mulqq%ddr)
    return
  end function

  function dd_mulqz (qa, zb)
    implicit none
    type (dd_complex):: dd_mulqz
    type (dd_real), intent (in):: qa
    type (dd_complex), intent (in):: zb
    real (ddknd) qt1(4)
    call ddmzc (qa%ddr, qt1)
    call ddcmul (qt1, zb%ddc, dd_mulqz%ddc)
    return
  end function

  function dd_muldq (da, qb)
    implicit none
    type (dd_real):: dd_muldq
    real (ddknd), intent (in):: da
    type (dd_real), intent (in):: qb
    call ddmuld (qb%ddr, da, dd_muldq%ddr)
    return
  end function

  function dd_mulqd (qa, db)
    implicit none
    type (dd_real):: dd_mulqd
    type (dd_real), intent (in):: qa
    real (ddknd), intent (in):: db
    call ddmuld (qa%ddr, db, dd_mulqd%ddr)
    return
  end function

  function dd_mulxq (xa, qb)
    implicit none
    type (dd_complex):: dd_mulxq
    complex (ddknd), intent (in):: xa
    type (dd_real), intent (in):: qb
    real (ddknd) qt1(4), qt2(4)
    call ddxzc (xa, qt1)
    call ddmzc (qb%ddr, qt2)
    call ddcmul (qt1, qt2, dd_mulxq%ddc)
    return
  end function

  function dd_mulqx (qa, xb)
    implicit none
    type (dd_complex):: dd_mulqx
    type (dd_real), intent (in):: qa
    complex (ddknd), intent (in):: xb
    real (ddknd) qt1(4), qt2(4)
    call ddmzc (qa%ddr, qt1)
    call ddxzc (xb, qt2)
    call ddcmul (qt1, qt2, dd_mulqx%ddc)
    return
  end function

  function dd_mulzq (za, qb)
    implicit none
    type (dd_complex):: dd_mulzq
    type (dd_complex), intent (in):: za
    type (dd_real), intent (in):: qb
    real (ddknd) qt1(4)
    call ddmzc (qb%ddr, qt1)
    call ddcmul (za%ddc, qt1, dd_mulzq%ddc)
    return
  end function

  function dd_mulzz (za, zb)
    implicit none
    type (dd_complex):: dd_mulzz
    type (dd_complex), intent (in):: za, zb
    call ddcmul (za%ddc, zb%ddc, dd_mulzz%ddc)
    return
  end function

  function dd_muldz (da, zb)
    implicit none
    type (dd_complex):: dd_muldz
    real (ddknd), intent (in):: da
    type (dd_complex), intent (in):: zb
    real (ddknd) qt1(4)
    complex (ddknd) xa
    xa = da
    call ddxzc (xa, qt1)
    call ddcmul (qt1, zb%ddc, dd_muldz%ddc)
    return
  end function

  function dd_mulzd (za, db)
    implicit none
    type (dd_complex):: dd_mulzd
    type (dd_complex), intent (in):: za
    real (ddknd), intent (in):: db
    real (ddknd) qt1(4)
    complex (ddknd) xb
    xb = db
    call ddxzc (xb, qt1)
    call ddcmul (za%ddc, qt1, dd_mulzd%ddc)
    return
  end function

  function dd_mulxz (xa, zb)
    implicit none
    type (dd_complex):: dd_mulxz
    complex (ddknd), intent (in):: xa
    type (dd_complex), intent (in):: zb
    real (ddknd) qt1(4)
    call ddxzc (xa, qt1)
    call ddcmul (qt1, zb%ddc, dd_mulxz%ddc)
    return
  end function

  function dd_mulzx (za, xb)
    implicit none
    type (dd_complex):: dd_mulzx
    type (dd_complex), intent (in):: za
    complex (ddknd), intent (in):: xb
    real (ddknd) qt1(4)
    call ddxzc (xb, qt1)
    call ddcmul (za%ddc, qt1, dd_mulzx%ddc)
    return
  end function

!  ddr divide routines.

  function dd_divqq (qa, qb)
    implicit none
    type (dd_real):: dd_divqq
    type (dd_real), intent (in):: qa, qb
    call dddiv (qa%ddr, qb%ddr, dd_divqq%ddr)
    return
  end function

  function dd_divqz (qa, zb)
    implicit none
    type (dd_complex):: dd_divqz
    type (dd_real), intent (in):: qa
    type (dd_complex), intent (in):: zb
    real (ddknd) qt1(4)
    call ddmzc (qa%ddr, qt1)
    call ddcdiv (qt1, zb%ddc, dd_divqz%ddc)
    return
  end function

  function dd_divdq (da, qb)
    implicit none
    type (dd_real):: dd_divdq
    real (ddknd), intent (in):: da
    type (dd_real), intent (in):: qb
    real (ddknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.d0
    call dddiv (qt1, qb%ddr, dd_divdq%ddr)
    return
  end function

  function dd_divqd (qa, db)
    implicit none
    type (dd_real):: dd_divqd
    type (dd_real), intent (in):: qa
    real (ddknd), intent (in):: db
    call dddivd (qa%ddr, db, dd_divqd%ddr)
    return
  end function

  function dd_divxq (xa, qb)
    implicit none
    type (dd_complex):: dd_divxq
    complex (ddknd), intent (in):: xa
    type (dd_real), intent (in):: qb
    real (ddknd) qt1(4), qt2(4)
    call ddxzc (xa, qt1)
    call ddmzc (qb%ddr, qt2)
    call ddcdiv (qt1, qt2, dd_divxq%ddc)
    return
  end function

  function dd_divqx (qa, xb)
    implicit none
    type (dd_complex):: dd_divqx
    type (dd_real), intent (in):: qa
    complex (ddknd), intent (in):: xb
    real (ddknd) qt1(4), qt2(4)
    call ddmzc (qa%ddr, qt1)
    call ddxzc (xb, qt2)
    call ddcdiv (qt1, qt2, dd_divqx%ddc)
    return
  end function

  function dd_divzq (za, qb)
    implicit none
    type (dd_complex):: dd_divzq
    type (dd_complex), intent (in):: za
    type (dd_real), intent (in):: qb
    real (ddknd) qt1(4)
    call ddmzc (qb%ddr, qt1)
    call ddcdiv (za%ddc, qt1, dd_divzq%ddc)
    return
  end function

  function dd_divzz (za, zb)
    implicit none
    type (dd_complex):: dd_divzz
    type (dd_complex), intent (in):: za, zb
    call ddcdiv (za%ddc, zb%ddc, dd_divzz%ddc)
    return
  end function

  function dd_divdz (da, zb)
    implicit none
    type (dd_complex):: dd_divdz
    real (ddknd), intent (in):: da
    type (dd_complex), intent (in):: zb
    real (ddknd) qt1(4)
    complex (ddknd) xa
    xa = da
    call ddxzc (xa, qt1)
    call ddcdiv (qt1, zb%ddc, dd_divdz%ddc)
    return
  end function

  function dd_divzd (za, db)
    implicit none
    type (dd_complex):: dd_divzd
    type (dd_complex), intent (in):: za
    real (ddknd), intent (in):: db
    real (ddknd) qt1(4)
    complex (ddknd) xb
    xb = db
    call ddxzc (xb, qt1)
    call ddcdiv (za%ddc, qt1, dd_divzd%ddc)
    return
  end function

  function dd_divxz (xa, zb)
    implicit none
    type (dd_complex):: dd_divxz
    complex (ddknd), intent (in):: xa
    type (dd_complex), intent (in):: zb
    real (ddknd) qt1(4)
    call ddxzc (xa, qt1)
    call ddcdiv (qt1, zb%ddc, dd_divxz%ddc)
    return
  end function

  function dd_divzx (za, xb)
    implicit none
    type (dd_complex):: dd_divzx
    type (dd_complex), intent (in):: za
    complex (ddknd), intent (in):: xb
    real (ddknd) qt1(4)
    call ddxzc (xb, qt1)
    call ddcdiv (za%ddc, qt1, dd_divzx%ddc)
    return
  end function

!  ddr exponentiation routines.

  function dd_expqq (qa, qb)
    implicit none
    type (dd_real):: dd_expqq
    type (dd_real), intent (in):: qa, qb
    real (ddknd) qt1(4), qt2(4)
    call ddlog (qa%ddr, qt1)
    call ddmul (qt1, qb%ddr, qt2)
    call ddexp (qt2, dd_expqq%ddr)
    return
  end function

  function dd_expqi (qa, ib)
    implicit none
    type (dd_real):: dd_expqi
    type (dd_real), intent (in):: qa
    integer, intent (in):: ib
    call ddnpwr (qa%ddr, ib, dd_expqi%ddr)
    return
  end function

  function dd_expdq (da, qb)
    implicit none
    type (dd_real):: dd_expdq
    real (ddknd), intent (in):: da
    type (dd_real), intent (in):: qb
    real (ddknd) qt1(4), qt2(4), qt3(4)
    qt1(1) = da
    qt1(2) = 0.d0
    call ddlog (qt1, qt2)
    call ddmul (qt2, qb%ddr, qt3)
    call ddexp (qt3, dd_expdq%ddr)
    return
    end function

  function dd_expqd (qa, db)
    implicit none
    type (dd_real):: dd_expqd
    type (dd_real), intent (in):: qa
    real (ddknd), intent (in):: db
    real (ddknd) qt1(4), qt2(4)
    call ddlog (qa%ddr, qt1)
    call ddmuld (qt1, db, qt2)
    call ddexp (qt2, dd_expqd%ddr)
    return
  end function

  function dd_expzi (za, ib)
    implicit none
    type (dd_complex):: dd_expzi
    type (dd_complex), intent (in):: za
    integer, intent (in):: ib
    call ddcpwr (za%ddc, ib, dd_expzi%ddc)
    return
  end function

  function dd_expzz (za, zb)
    implicit none
    type (dd_complex):: dd_expzz
    type (dd_complex), intent (in):: za, zb
    real (ddknd) r1(4), r2(4), r3(4), r4(4), r5(4), r6(4)
    call ddmul (za%ddc(1), za%ddc(1), r1(1))
    call ddmul (za%ddc(3), za%ddc(3), r2(1))
    call ddadd (r1(1), r2(1), r3(1))
    call ddlog (r3(1), r4(1))
    call ddmuld (r4(1), 0.5d0, r5(1))
    call ddmul (zb%ddc(1), r5(1), r1(1))
    call ddang (za%ddc(1), za%ddc(3), r2(1))
    call ddmul (r2(1), zb%ddc(3), r3(1))
    call ddsub (r1(1), r3(1), r4(1))
    call ddexp (r4(1), r1(1))
    call ddmul (zb%ddc(3), r5(1), r3(1))
    call ddmul (zb%ddc(1), r2(1), r4(1))
    call ddadd (r3(1), r4(1), r6(1))
    call ddcssnr (r6(1), r3(1), r4(1))
    call ddmul (r1(1), r3(1), dd_expzz%ddc(1))
    call ddmul (r1(1), r4(1), dd_expzz%ddc(3))
    return
  end function

  function dd_exprz (qa, zb)
    implicit none
    type (dd_complex):: dd_exprz
    type (dd_real), intent (in):: qa
    type (dd_complex), intent (in):: zb
    real (ddknd) r1(4), r2(4), r3(4), r4(4), r5(4)
    call ddlog (qa%ddr(1), r2(1))
    call ddmul (r2(1), zb%ddc(1), r3(1))
    call ddexp (r3(1), r1(1))
    call ddlog (qa%ddr(1), r2(1))
    call ddmul (r2(1), zb%ddc(3), r3(1))
    call ddcssnr (r3(1), r4(1), r5(1))
    call ddmul (r1(1), r4(1), dd_exprz%ddc(1))
    call ddmul (r1(1), r5(1), dd_exprz%ddc(3))
    return
  end function

  function dd_expzr (za, qb)
    implicit none
    type (dd_complex):: dd_expzr
    type (dd_complex), intent (in):: za
    type (dd_real), intent (in):: qb
    real (ddknd) r1(4), r2(4), r3(4), r4(4), r5(4)
    call ddmul (za%ddc(1), za%ddc(1), r1(1))
    call ddmul (za%ddc(3), za%ddc(3), r2(1))
    call ddadd (r1(1), r2(1), r3(1))
    call ddlog (r3(1), r4(1))
    call ddmuld (r4(1), 0.5d0, r5(1))
    call ddmul (r5(1), qb%ddr(1), r1(1))
    call ddexp (r1(1), r2(1))
    call ddang (za%ddc(1), za%ddc(3), r3(1))
    call ddmul (qb%ddr(1), r3(1), r1(1))
    call ddcssnr (r1(1), r4(1), r5(1)) 
    call ddmul (r2(1), r4(1), dd_expzr%ddc(1))
    call ddmul (r2(1), r5(1), dd_expzr%ddc(3))
    return
  end function

!  Equality test routines.

  function dd_eqtqq (qa, qb)
    implicit none
    logical dd_eqtqq
    type (dd_real), intent (in):: qa, qb
    integer ic
    call ddcpr (qa%ddr, qb%ddr, ic)
    if (ic == 0) then
      dd_eqtqq = .true.
    else
      dd_eqtqq = .false.
    endif
    return
  end function

  function dd_eqtqz (qa, zb)
    implicit none
    logical dd_eqtqz
    type (dd_real), intent (in):: qa
    type (dd_complex), intent (in):: zb
    integer ic1, ic2
    real (ddknd) qt1(4)
    call ddmzc (qa%ddr, qt1)
    call ddcpr (qt1, zb%ddc, ic1)
    call ddcpr (qt1(3:4), zb%ddc(3:4), ic2)
    if (ic1 == 0 .and. ic2 == 0) then
      dd_eqtqz = .true.
    else
      dd_eqtqz = .false.
    endif
    return
  end function

  function dd_eqtdq (da, qb)
    implicit none
    logical dd_eqtdq
    real (ddknd), intent (in):: da
    type (dd_real), intent (in):: qb
    integer ic
    real (ddknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.d0
    call ddcpr (qt1, qb%ddr, ic)
    if (ic == 0) then
      dd_eqtdq = .true.
    else
      dd_eqtdq = .false.
    endif
    return
  end function

  function dd_eqtqd (qa, db)
    implicit none
    logical dd_eqtqd
    type (dd_real), intent (in):: qa
    real (ddknd), intent (in):: db
    integer ic
    real (ddknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.d0
    call ddcpr (qa%ddr, qt1, ic)
    if (ic == 0) then
      dd_eqtqd = .true.
    else
      dd_eqtqd = .false.
    endif
    return
  end function

  function dd_eqtxq (xa, qb)
    implicit none
    logical dd_eqtxq
    complex (ddknd), intent (in):: xa
    type (dd_real), intent (in):: qb
    integer ic1, ic2
    real (ddknd) qt1(4), qt2(4)
    call ddxzc (xa, qt1)
    call ddmzc (qb%ddr, qt2)
    call ddcpr (qt1, qt2, ic1)
    call ddcpr (qt1(3:4), qt2(3:4), ic2)
    if (ic1 == 0 .and. ic2 == 0) then
      dd_eqtxq = .true.
    else
      dd_eqtxq = .false.
    endif
    return
  end function

  function dd_eqtqx (qa, xb)
    implicit none
    logical dd_eqtqx
    type (dd_real), intent (in):: qa
    complex (ddknd), intent (in):: xb
    integer ic1, ic2
    real (ddknd) qt1(4), qt2(4)
    call ddmzc (qa%ddr, qt1)
    call ddxzc (xb, qt2)
    call ddcpr (qt1, qt2, ic1)
    call ddcpr (qt1(3:4), qt2(3:4), ic2)
    if (ic1 == 0 .and. ic2 == 0) then
      dd_eqtqx = .true.
    else
      dd_eqtqx = .false.
    endif
    return
  end function

  function dd_eqtzq (za, qb)
    implicit none
    logical dd_eqtzq
    type (dd_complex), intent (in):: za
    type (dd_real), intent (in):: qb
    integer ic1, ic2
    real (ddknd) qt1(4)
    call ddmzc (qb%ddr, qt1)
    call ddcpr (za%ddc, qt1, ic1)
    call ddcpr (za%ddc(3:4), qt1(3:4), ic2)
    if (ic1 == 0 .and. ic2 == 0) then
      dd_eqtzq = .true.
    else
      dd_eqtzq = .false.
    endif
    return
  end function

  function dd_eqtzz (za, zb)
    implicit none
    logical dd_eqtzz
    type (dd_complex), intent (in):: za, zb
    integer ic1, ic2
    call ddcpr (za%ddc, zb%ddc, ic1)
    call ddcpr (za%ddc(3:4), zb%ddc(3:4), ic2)
    if (ic1 == 0 .and. ic2 == 0) then
      dd_eqtzz = .true.
    else
      dd_eqtzz = .false.
    endif
    return
  end function

  function dd_eqtdz (da, zb)
    implicit none
    logical dd_eqtdz
    real (ddknd), intent (in):: da
    type (dd_complex), intent (in):: zb
    integer ic1, ic2
    real (ddknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.d0
    call ddcpr (qt1, zb%ddc, ic1)
    call ddcpr (qt1(3:4), zb%ddc(3:4), ic2)
    if (ic1 == 0 .and. ic2 == 0) then
      dd_eqtdz = .true.
    else
      dd_eqtdz = .false.
    endif
    return
  end function

  function dd_eqtzd (za, db)
    implicit none
    logical dd_eqtzd
    type (dd_complex), intent (in):: za
    real (ddknd), intent (in):: db
    integer ic1, ic2
    real (ddknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.d0
    call ddcpr (za%ddc, qt1, ic1)
    call ddcpr (za%ddc(3:4), qt1(3:4), ic2)
    if (ic1 == 0 .and. ic2 == 0) then
      dd_eqtzd = .true.
    else
      dd_eqtzd = .false.
    endif
    return
  end function

  function dd_eqtxz (xa, zb)
    implicit none
    logical dd_eqtxz
    complex (ddknd), intent (in):: xa
    type (dd_complex), intent (in):: zb
    integer ic1, ic2
    real (ddknd) qt1(4)
    call ddxzc (xa, qt1)
    call ddcpr (qt1, zb%ddc, ic1)
    call ddcpr (qt1(3:4), zb%ddc(3:4), ic2)
    if (ic1 == 0 .and. ic2 == 0) then
      dd_eqtxz = .true.
    else
      dd_eqtxz = .false.
    endif
    return
  end function

  function dd_eqtzx (za, xb)
    implicit none
    logical dd_eqtzx
    type (dd_complex), intent (in):: za
    complex (ddknd), intent (in):: xb
    integer ic1, ic2
    real (ddknd) qt1(4)
    call ddxzc (xb, qt1)
    call ddcpr (za%ddc, qt1, ic1)
    call ddcpr (za%ddc(3:4), qt1(3:4), ic2)
    if (ic1 == 0 .and. ic2 == 0) then
      dd_eqtzx = .true.
    else
      dd_eqtzx = .false.
    endif
    return
  end function

!  Inequality test routines.

  function dd_netqq (qa, qb)
    implicit none
    logical dd_netqq
    type (dd_real), intent (in):: qa, qb
    integer ic
    call ddcpr (qa%ddr, qb%ddr, ic)
    if (ic /= 0) then
      dd_netqq = .true.
    else
      dd_netqq = .false.
    endif
    return
  end function

  function dd_netqz (qa, zb)
    implicit none
    logical dd_netqz
    type (dd_real), intent (in):: qa
    type (dd_complex), intent (in):: zb
    integer ic1, ic2
    real (ddknd) qt1(4)
    call ddmzc (qa%ddr, qt1)
    call ddcpr (qt1, zb%ddc, ic1)
    call ddcpr (qt1(3:4), zb%ddc(3:4), ic2)
    if (ic1 /= 0 .or. ic2 /= 0) then
      dd_netqz = .true.
    else
      dd_netqz = .false.
    endif
    return
  end function

  function dd_netdq (da, qb)
    implicit none
    logical dd_netdq
    real (ddknd), intent (in):: da
    type (dd_real), intent (in):: qb
    integer ic
    real (ddknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.d0
    call ddcpr (qt1, qb%ddr, ic)
    if (ic /= 0) then
      dd_netdq = .true.
    else
      dd_netdq = .false.
    endif
    return
  end function

  function dd_netqd (qa, db)
    implicit none
    logical dd_netqd
    type (dd_real), intent (in):: qa
    real (ddknd), intent (in):: db
    integer ic
    real (ddknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.d0
    call ddcpr (qa%ddr, qt1, ic)
    if (ic /= 0) then
      dd_netqd = .true.
    else
      dd_netqd = .false.
    endif
    return
  end function

  function dd_netxq (xa, qb)
    implicit none
    logical dd_netxq
    complex (ddknd), intent (in):: xa
    type (dd_real), intent (in):: qb
    integer ic1, ic2
    real (ddknd) qt1(4), qt2(4)
    call ddxzc (xa, qt1)
    call ddmzc (qb%ddr, qt2)
    call ddcpr (qt1, qt2, ic1)
    call ddcpr (qt1(3:4), qt2(3:4), ic2)
    if (ic1 /= 0 .or. ic2 /= 0) then
      dd_netxq = .true.
    else
      dd_netxq = .false.
    endif
    return
  end function

  function dd_netqx (qa, xb)
    implicit none
    logical dd_netqx
    type (dd_real), intent (in):: qa
    complex (ddknd), intent (in):: xb
    integer ic1, ic2
    real (ddknd) qt1(4), qt2(4)
    call ddmzc (qa%ddr, qt1)
    call ddxzc (xb, qt2)
    call ddcpr (qt1, qt2, ic1)
    call ddcpr (qt1(3:4), qt2(3:4), ic2)
    if (ic1 /= 0 .or. ic2 /= 0) then
      dd_netqx = .true.
    else
      dd_netqx = .false.
    endif
    return
  end function

  function dd_netzq (za, qb)
    implicit none
    logical dd_netzq
    type (dd_complex), intent (in):: za
    type (dd_real), intent(in):: qb
    integer ic1, ic2
    real (ddknd) qt1(4)
    call ddmzc (qb%ddr, qt1)
    call ddcpr (za%ddc, qt1, ic1)
    call ddcpr (za%ddc(3:4), qt1(3:4), ic2)
    if (ic1 /= 0 .or. ic2 /= 0) then
      dd_netzq = .true.
    else
      dd_netzq = .false.
    endif
    return
  end function

  function dd_netzz (za, zb)
    implicit none
    logical dd_netzz
    type (dd_complex), intent (in):: za, zb
    integer ic1, ic2
    call ddcpr (za%ddc, zb%ddc, ic1)
    call ddcpr (za%ddc(3:4), zb%ddc(3:4), ic2)
    if (ic1 /= 0 .or. ic2 /= 0) then
      dd_netzz = .true.
    else
      dd_netzz = .false.
    endif
    return
  end function

  function dd_netdz (da, zb)
    implicit none
    logical dd_netdz
    real (ddknd), intent (in):: da
    type (dd_complex), intent (in):: zb
    integer ic1, ic2
    real (ddknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.d0
    call ddcpr (qt1, zb%ddc, ic1)
    call ddcpr (qt1(3:4), zb%ddc(3:4), ic2)
    if (ic1 /= 0 .or. ic2 /= 0) then
      dd_netdz = .true.
    else
      dd_netdz = .false.
    endif
    return
  end function

  function dd_netzd (za, db)
    implicit none
    logical dd_netzd
    type (dd_complex), intent (in):: za
    real (ddknd), intent (in):: db
    integer ic1, ic2
    real (ddknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.d0
    call ddcpr (za%ddc, qt1, ic1)
    call ddcpr (za%ddc(3:4), qt1(3:4), ic2)
    if (ic1 /= 0 .or. ic2 /= 0) then
      dd_netzd = .true.
    else
      dd_netzd = .false.
    endif
    return
  end function

  function dd_netxz (xa, zb)
    implicit none
    logical dd_netxz
    complex (ddknd), intent (in):: xa
    type (dd_complex), intent (in):: zb
    integer ic1, ic2
    real (ddknd) qt1(4)
    call ddxzc (xa, qt1)
    call ddcpr (qt1, zb%ddc, ic1)
    call ddcpr (qt1(3:4), zb%ddc(3:4), ic2)
    if (ic1 /= 0 .or. ic2 /= 0) then
      dd_netxz = .true.
    else
      dd_netxz = .false.
    endif
    return
  end function

  function dd_netzx (za, xb)
    implicit none
    logical dd_netzx
    type (dd_complex), intent (in):: za
    complex (ddknd), intent (in):: xb
    integer ic1, ic2
    real (ddknd) qt1(4)
    call ddxzc (xb, qt1)
    call ddcpr (za%ddc, qt1, ic1)
    call ddcpr (za%ddc(3:4), qt1(3:4), ic2)
    if (ic1 /= 0 .or. ic2 /= 0) then
      dd_netzx = .true.
    else
      dd_netzx = .false.
    endif
    return
  end function

!  Less-than-or-equal test routines.

  function dd_letqq (qa, qb)
    implicit none
    logical dd_letqq
    type (dd_real), intent (in):: qa, qb
    integer ic
    call ddcpr (qa%ddr, qb%ddr, ic)
    if (ic <= 0) then
      dd_letqq = .true.
    else
      dd_letqq = .false.
    endif
    return
  end function

  function dd_letdq (da, qb)
    implicit none
    logical dd_letdq
    real (ddknd), intent (in):: da
    type (dd_real), intent (in):: qb
    integer ic
    real (ddknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.d0
    call ddcpr (qt1, qb%ddr, ic)
    if (ic <= 0) then
      dd_letdq = .true.
    else
      dd_letdq = .false.
    endif
    return
  end function

  function dd_letqd (qa, db)
    implicit none
    logical dd_letqd
    type (dd_real), intent (in):: qa
    real (ddknd), intent (in):: db
    integer ic
    real (ddknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.d0
    call ddcpr (qa%ddr, qt1, ic)
    if (ic <= 0) then
      dd_letqd = .true.
    else
      dd_letqd = .false.
    endif
    return
  end function

!  Greater-than-or-equal test routines.

  function dd_getqq (qa, qb)
    implicit none
    logical dd_getqq
    type (dd_real), intent (in):: qa, qb
    integer ic
    call ddcpr (qa%ddr, qb%ddr, ic)
    if (ic >= 0) then
      dd_getqq = .true.
    else
      dd_getqq = .false.
    endif
    return
  end function

  function dd_getdq (da, qb)
    implicit none
    logical dd_getdq
    real (ddknd), intent (in):: da
    type (dd_real), intent (in):: qb
    integer ic
    real (ddknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.d0
    call ddcpr (qt1, qb%ddr, ic)
    if (ic >= 0) then
      dd_getdq = .true.
    else
      dd_getdq = .false.
    endif
    return
  end function

  function dd_getqd (qa, db)
    implicit none
    logical dd_getqd
    type (dd_real), intent (in):: qa
    real (ddknd), intent (in):: db
    integer ic
    real (ddknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.d0
    call ddcpr (qa%ddr, qt1, ic)
    if (ic >= 0) then
      dd_getqd = .true.
    else
      dd_getqd = .false.
    endif
    return
  end function

!  Less-than test routines.

  function dd_lttqq (qa, qb)
    implicit none
    logical dd_lttqq
    type (dd_real), intent (in):: qa, qb
    integer ic
    call ddcpr (qa%ddr, qb%ddr, ic)
    if (ic < 0) then
      dd_lttqq = .true.
    else
      dd_lttqq = .false.
    endif
    return
  end function

  function dd_lttdq (da, qb)
    implicit none
    logical dd_lttdq
    real (ddknd), intent (in):: da
    type (dd_real), intent (in):: qb
    integer ic
    real (ddknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.d0
    call ddcpr (qt1, qb%ddr, ic)
    if (ic < 0) then
      dd_lttdq = .true.
    else
      dd_lttdq = .false.
    endif
    return
  end function

  function dd_lttqd (qa, db)
    implicit none
    logical dd_lttqd
    type (dd_real), intent (in):: qa
    real (ddknd), intent (in):: db
    integer ic
    real (ddknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.d0
    call ddcpr (qa%ddr, qt1, ic)
    if (ic < 0) then
      dd_lttqd = .true.
    else
      dd_lttqd = .false.
    endif
    return
  end function

!  Greater-than test routines.

  function dd_gttqq (qa, qb)
    implicit none
    logical dd_gttqq
    type (dd_real), intent (in):: qa, qb
    integer ic
    call ddcpr (qa%ddr, qb%ddr, ic)
    if (ic > 0) then
      dd_gttqq = .true.
    else
      dd_gttqq = .false.
    endif
    return
  end function

  function dd_gttdq (da, qb)
    implicit none
    logical dd_gttdq
    real (ddknd), intent (in):: da
    type (dd_real), intent (in):: qb
    integer ic
    real (ddknd) qt1(4)
    qt1(1) = da
    qt1(2) = 0.d0
    call ddcpr (qt1, qb%ddr, ic)
    if (ic > 0) then
      dd_gttdq = .true.
    else
      dd_gttdq = .false.
    endif
    return
  end function

  function dd_gttqd (qa, db)
    implicit none
    logical dd_gttqd
    type (dd_real), intent (in):: qa
    real (ddknd), intent (in):: db
    integer ic
    real (ddknd) qt1(4)
    qt1(1) = db
    qt1(2) = 0.d0
    call ddcpr (qa%ddr, qt1, ic)
    if (ic > 0) then
      dd_gttqd = .true.
    else
      dd_gttqd = .false.
    endif
    return
  end function

!  Other DD subroutines and functions, in alphabetical order.

  function dd_absq (qa)
    implicit none
    type (dd_real):: dd_absq
    type (dd_real), intent (in):: qa
    call ddeq (qa%ddr, dd_absq%ddr)
    if (qa%ddr(1) >= 0.d0) then
      dd_absq%ddr(1) = qa%ddr(1)
      dd_absq%ddr(2) = qa%ddr(2)
    else
      dd_absq%ddr(1) = - qa%ddr(1)
      dd_absq%ddr(2) = - qa%ddr(2)
    endif
    return
  end function

  function dd_absz (za)
    implicit none
    type (dd_real):: dd_absz
    type (dd_complex), intent (in):: za
    real (ddknd) qt1(4), qt2(4), qt3(4)
    call ddmul (za%ddc, za%ddc, qt1)
    call ddmul (za%ddc(3:4), za%ddc(3:4), qt2)
    call ddadd (qt1, qt2, qt3)
    call ddsqrt (qt3, dd_absz%ddr)
    return
  end function

  function dd_acos (qa)
    implicit none
    type (dd_real):: dd_acos
    type (dd_real), intent (in):: qa
    real (ddknd) qt1(4), qt2(4), qt3(4)
    qt1(1) = 1.d0
    qt1(2) = 0.d0
    call ddmul (qa%ddr, qa%ddr, qt2)
    call ddsub (qt1, qt2, qt3)
    call ddsqrt (qt3, qt1)
    call ddang (qa%ddr, qt1, dd_acos%ddr)
    return
  end function

  function dd_acosh (qa)
    implicit none
    type (dd_real):: dd_acosh
    type (dd_real), intent (in):: qa
    call ddacosh (qa%ddr, dd_acosh%ddr)
    return
  end function

  function dd_aimag (za)
    implicit none
    type (dd_real):: dd_aimag
    type (dd_complex), intent (in):: za
    call ddeq (za%ddc(3:4), dd_aimag%ddr)
    return
  end function

  function dd_agm (qa, qb)
    implicit none
    type (dd_real):: dd_agm
    type (dd_real), intent (in):: qa, qb
    call ddagmr (qa%ddr, qb%ddr, dd_agm%ddr)
    return
  end function

  function dd_aint (qa)
    implicit none
    type (dd_real):: dd_aint
    type (dd_real), intent (in):: qa
    real (ddknd) qt1(4)
    call ddinfr (qa%ddr, dd_aint%ddr, qt1)
    return
  end function

  function dd_anint (qa)
    implicit none
    type (dd_real):: dd_anint
    type (dd_real), intent (in):: qa
    call ddnint (qa%ddr, dd_anint%ddr)
    return
  end function

  function dd_asin (qa)
    implicit none
    type (dd_real):: dd_asin
    type (dd_real), intent (in):: qa
    real (ddknd) qt1(4), qt2(4), qt3(4)
    qt1(1) = 1.d0
    qt1(2) = 0.d0
    call ddmul (qa%ddr, qa%ddr, qt2)
    call ddsub (qt1, qt2, qt3)
    call ddsqrt (qt3, qt1)
    call ddang (qt1, qa%ddr, dd_asin%ddr)
    return
  end function

  function dd_asinh (qa)
    implicit none
    type (dd_real):: dd_asinh
    type (dd_real), intent (in):: qa
    call ddasinh (qa%ddr, dd_asinh%ddr)
    return
  end function

  function dd_atan (qa)
    implicit none
    type (dd_real):: dd_atan
    type (dd_real), intent (in):: qa
    real (ddknd) qt1(4)
    qt1(1) = 1.d0
    qt1(2) = 0.d0
    call ddang (qt1, qa%ddr, dd_atan%ddr)
    return
  end function

  function dd_atan2 (qa, qb)
    implicit none
    type (dd_real):: dd_atan2
    type (dd_real), intent (in):: qa, qb
    call ddang (qb%ddr, qa%ddr, dd_atan2%ddr)
    return
  end function

  function dd_atanh (qa)
    implicit none
    type (dd_real):: dd_atanh
    type (dd_real), intent (in):: qa
    call ddatanh (qa%ddr, dd_atanh%ddr)
    return
  end function

  function dd_atoq (aa)
    implicit none
    character(*), intent (in):: aa
    type (dd_real):: dd_atoq
    character(120) t
    t = aa
    call ddinpc (t, dd_atoq%ddr)
    return
  end function

  subroutine dd_berne (nb, rb)
    implicit none
    integer, intent (in):: nb
    type (dd_real), intent (out):: rb(nb)
    call ddberner (nb, rb(1)%ddr)
    return
  end subroutine

  function dd_bessel_i (qa, ra)
    implicit none
    type (dd_real):: dd_bessel_i
    type (dd_real), intent(in):: qa, ra
    call ddbesselir (qa%ddr, ra%ddr, dd_bessel_i%ddr)
    return
  end function

  function dd_bessel_in (nu, ra)
    implicit none
    type (dd_real):: dd_bessel_in
    integer, intent(in):: nu
    type (dd_real), intent(in):: ra
    call ddbesselinr (nu, ra%ddr, dd_bessel_in%ddr)
    return
  end function

  function dd_bessel_j (qa, ra)
    implicit none
    type (dd_real):: dd_bessel_j
    type (dd_real), intent(in):: qa, ra
    call ddbesseljr (qa%ddr, ra%ddr, dd_bessel_j%ddr)
    return
  end function

  function dd_bessel_jn (nu, ra)
    implicit none
    type (dd_real):: dd_bessel_jn
    integer, intent(in):: nu
    type (dd_real), intent(in):: ra
    call ddbesseljnr (nu, ra%ddr, dd_bessel_jn%ddr)
    return
  end function

  function dd_bessel_j0 (ra)
    implicit none
    type (dd_real):: dd_bessel_j0
    integer:: nu
    type (dd_real), intent(in):: ra
    nu = 0
    call ddbesseljnr (nu, ra%ddr, dd_bessel_j0%ddr)
    return
  end function

  function dd_bessel_j1 (ra)
    implicit none
    type (dd_real):: dd_bessel_j1
    integer:: nu
    type (dd_real), intent(in):: ra
    nu = 1
    call ddbesseljnr (nu, ra%ddr, dd_bessel_j1%ddr)
    return
  end function

  function dd_bessel_k (qa, ra)
    implicit none
    type (dd_real):: dd_bessel_k
    type (dd_real), intent(in):: qa, ra
    call ddbesselkr (qa%ddr, ra%ddr, dd_bessel_k%ddr)
    return
  end function

  function dd_bessel_kn (nu, ra)
    implicit none
    type (dd_real):: dd_bessel_kn
    integer, intent(in):: nu 
    type (dd_real), intent(in):: ra
    call ddbesselknr (nu, ra%ddr, dd_bessel_kn%ddr)
    return
  end function

  function dd_bessel_y (qa, ra)
    implicit none
    type (dd_real):: dd_bessel_y
    type (dd_real), intent(in):: qa, ra
    call ddbesselyr (qa%ddr, ra%ddr, dd_bessel_y%ddr)
    return
  end function

  function dd_bessel_yn (nu, ra)
    implicit none
    type (dd_real):: dd_bessel_yn
    integer, intent(in):: nu 
    type (dd_real), intent(in):: ra
    call ddbesselynr (nu, ra%ddr, dd_bessel_yn%ddr)
    return
  end function

  function dd_bessel_y0 (ra)
    implicit none
    type (dd_real):: dd_bessel_y0
    integer:: nu
    type (dd_real), intent(in):: ra
    nu = 0
    call ddbesselynr (nu, ra%ddr, dd_bessel_y0%ddr)
    return
  end function

  function dd_bessel_y1 (ra)
    implicit none
    type (dd_real):: dd_bessel_y1
    integer:: nu
    type (dd_real), intent(in):: ra
    nu = 1
    call ddbesselynr (nu, ra%ddr, dd_bessel_y1%ddr)
    return
  end function

  function dd_conjg (za)
    implicit none
    type (dd_complex):: dd_conjg
    type (dd_complex), intent (in):: za
    call ddceq (za%ddc, dd_conjg%ddc)
    dd_conjg%ddc(3) = - za%ddc(3)
    dd_conjg%ddc(4) = - za%ddc(4)
    return
  end function

  function dd_cos (qa)
    implicit none
    type (dd_real):: dd_cos
    type (dd_real), intent (in):: qa
    real (ddknd) qt1(4)
    call ddcssnr (qa%ddr, dd_cos%ddr, qt1)
    return
  end function

  function dd_cosz (za)
    implicit none
    type (dd_complex):: dd_cosz
    type (dd_complex), intent (in):: za
    real (ddknd) qt1(4), qt2(4), qt3(4), qt4(4), qt5(4), qt6(4)
    call ddeq (za%ddc(3:4), qt2)
    qt2(1) = - qt2(1)
    qt2(2) = - qt2(2)
    call ddexp (qt2, qt1)
    qt3(1) = 1.d0
    qt3(2) = 0.d0
    call dddiv (qt3, qt1, qt2)
    call ddcssnr (za%ddc, qt3, qt4)
    call ddadd (qt1, qt2, qt5)
    call ddmuld (qt5, 0.5d0, qt6)
    call ddmul (qt6, qt3, dd_cosz%ddc)
    call ddsub (qt1, qt2, qt5)
    call ddmuld (qt5, 0.5d0, qt6)
    call ddmul (qt6, qt4, dd_cosz%ddc(3:4))
    return
  end function

  function dd_cosh (qa)
    implicit none
    type (dd_real):: dd_cosh
    type (dd_real), intent (in):: qa
    real (ddknd) qt1(4)
    call ddcsshr (qa%ddr, dd_cosh%ddr, qt1)
    return
  end function

  subroutine dd_cssh (qa, qb, qc)
    implicit none
    type (dd_real), intent (in):: qa
    type (dd_real), intent (out):: qb, qc
    call ddcsshr (qa%ddr, qb%ddr, qc%ddr)
    return
  end subroutine

  subroutine dd_cssn (qa, qb, qc)
    implicit none
    type (dd_real), intent (in):: qa
    type (dd_real), intent (out):: qb, qc
    call ddcssnr (qa%ddr, qb%ddr, qc%ddr)
    return
  end subroutine

  function dd_ddtoz (da, db)
    implicit none
    type (dd_complex):: dd_ddtoz
    real (ddknd), intent (in):: da, db
    complex (ddknd) xa
    xa = cmplx (da, db, ddknd)
    call ddxzc (xa, dd_ddtoz%ddc)
    return
  end function

  function dd_digamma_be (nb, rb, rc)
    implicit none
    integer, intent (in):: nb
    type (dd_real):: dd_digamma_be
    type (dd_real), intent (in):: rb(nb), rc
    call dddigammabe (nb, rb(1)%ddr, rc%ddr, dd_digamma_be%ddr)
    return
  end function

  function dd_dtoq (da)
    implicit none
    type (dd_real):: dd_dtoq
    real (ddknd), intent (in):: da
    dd_dtoq%ddr(1) = da
    dd_dtoq%ddr(2) = 0.d0
    return
  end function

  function dd_dtoz (da)
    implicit none
    type (dd_complex):: dd_dtoz
    real (ddknd), intent (in):: da
    complex (ddknd) xa
    xa = da
    call ddxzc (xa, dd_dtoz%ddc)
    return
  end function

  subroutine dd_eform (qa, n1, n2, b)
    implicit none
    type (dd_real), intent (in):: qa
    integer, intent (in):: n1, n2
    character(1), intent (out):: b(n1)
    call ddeformat (qa%ddr, n1, n2, b)
    return
  end subroutine

  function dd_egamma ()
    implicit none
    type (dd_real):: dd_egamma
    call ddegamc (dd_egamma%ddr)
    return
  end function

  function dd_erf (qa)
    implicit none
    type (dd_real):: dd_erf
    type (dd_real), intent (in):: qa
    call dderfr (qa%ddr, dd_erf%ddr)
    return
  end function

  function dd_erfc (qa)
    implicit none
    type (dd_real):: dd_erfc
    type (dd_real), intent (in):: qa
    call dderfcr (qa%ddr, dd_erfc%ddr)
    return
  end function

  function dd_exp (qa)
    implicit none
    type (dd_real):: dd_exp
    type (dd_real), intent (in):: qa
    call ddexp (qa%ddr, dd_exp%ddr)
    return
  end function

  function dd_expint (qa)
    implicit none
    type (dd_real):: dd_expint
    type (dd_real), intent (in):: qa
    call ddexpint (qa%ddr, dd_expint%ddr)
    return
  end function

  function dd_expz (za)
    implicit none
    type (dd_complex):: dd_expz
    type (dd_complex), intent (in):: za
    real (ddknd) qt1(4), qt2(4), qt3(4)
    call ddexp (za%ddc, qt1)
    call ddcssnr (za%ddc(3:4), qt2, qt3)
    call ddmul (qt1, qt2, dd_expz%ddc)
    call ddmul (qt1, qt3, dd_expz%ddc(3:4))
    return
  end function

  subroutine dd_fform (qa, n1, n2, b)
    implicit none
    type (dd_real), intent (in):: qa
    integer, intent (in):: n1, n2
    character(1), intent (out):: b(n1)
    call ddfformat (qa%ddr, n1, n2, b)
    return
  end subroutine

  function dd_gamma (qa)
    implicit none
    type (dd_real):: dd_gamma
    type (dd_real), intent (in):: qa
    call ddgammar (qa%ddr, dd_gamma%ddr)
    return
  end function

  function dd_hurwitz_zetan (ia, rb)
    implicit none
    type (dd_real):: dd_hurwitz_zetan
    integer, intent (in):: ia
    type (dd_real), intent (in):: rb
    call ddhurwitzzetan (ia, rb%ddr, dd_hurwitz_zetan%ddr)
    return
  end function

  function dd_hurwitz_zetan_be (nb, rb, is, aa)
    implicit none
    type (dd_real):: dd_hurwitz_zetan_be
    integer, intent (in):: nb, is
    type (dd_real), intent (in):: rb(nb), aa
    call ddhurwitzzetanbe (nb, rb(1)%ddr, is, aa%ddr, &
      dd_hurwitz_zetan_be%ddr)
    return
  end function

  function dd_hypergeom_pfq (np, nq, aa, bb, xx)
    implicit none
    type (dd_real):: dd_hypergeom_pfq
    integer, intent (in):: np, nq
    type (dd_real), intent (in):: aa(np), bb(nq), xx
    call ddhypergeompfq (np, nq, aa(1)%ddr, bb(1)%ddr, &
      xx%ddr, dd_hypergeom_pfq%ddr)
    return
  end function

  function dd_hypot (ra, rb)
    implicit none
    type (dd_real):: dd_hypot
    type (dd_real), intent (in):: ra, rb
    type (dd_real) r1, r2, r3
    call ddmul (ra%ddr, ra%ddr, r1%ddr)
    call ddmul (rb%ddr, rb%ddr, r2%ddr)
    call ddadd (r1%ddr, r2%ddr, r3%ddr)
    call ddsqrt (r3%ddr, dd_hypot%ddr)
    return
  end function

  function dd_incgamma (ra, rb)
    implicit none
    type (dd_real):: dd_incgamma
    type (dd_real), intent (in):: ra, rb
    call ddincgammar (ra%ddr, rb%ddr, dd_incgamma%ddr)
    return
  end function

  subroutine dd_inpq (iu, q1, q2, q3, q4, q5)
    implicit none
    integer, intent (in):: iu
    type (dd_real), intent (out):: q1, q2, q3, q4, q5
    optional:: q2, q3, q4, q5
    call ddinp (iu, q1%ddr)
    if (present (q2)) call ddinp (iu, q2%ddr)
    if (present (q3)) call ddinp (iu, q3%ddr)
    if (present (q4)) call ddinp (iu, q4%ddr)
    if (present (q5)) call ddinp (iu, q5%ddr)
    return
  end subroutine

  subroutine dd_inpz (iu, z1, z2, z3, z4, z5)
    implicit none
    integer, intent (in):: iu
    type (dd_complex), intent (out):: z1, z2, z3, z4, z5
    optional:: z2, z3, z4, z5
    call ddinp (iu, z1%ddc)
    call ddinp (iu, z1%ddc(3:4))
    if (present (z2)) call ddinp (iu, z2%ddc)
    if (present (z2)) call ddinp (iu, z2%ddc(3:4))
    if (present (z3)) call ddinp (iu, z3%ddc)
    if (present (z3)) call ddinp (iu, z3%ddc(3:4))
    if (present (z4)) call ddinp (iu, z4%ddc)
    if (present (z4)) call ddinp (iu, z4%ddc(3:4))
    if (present (z5)) call ddinp (iu, z5%ddc)
    if (present (z5)) call ddinp (iu, z5%ddc(3:4))
    return
  end subroutine

  function dd_itoq (ia)
    implicit none
    type (dd_real):: dd_itoq
    integer, intent (in):: ia
    dd_itoq%ddr(1) = ia
    dd_itoq%ddr(2) = 0.d0
    return
  end function

  function dd_log (qa)
    implicit none
    type (dd_real):: dd_log
    type (dd_real), intent (in):: qa
    call ddlog (qa%ddr, dd_log%ddr)
    return
  end function

  function dd_logz (za)
    implicit none
    type (dd_complex):: dd_logz
    type (dd_complex), intent (in):: za
    real (ddknd) qt1(4), qt2(4), qt3(4), qt4(4)
    call ddmul (za%ddc, za%ddc, qt1)
    call ddmul (za%ddc(3:4), za%ddc(3:4), qt2)
    call ddadd (qt1, qt2, qt3)
    call ddlog (qt3, qt4)
    call ddmuld (qt4, 0.5d0, dd_logz%ddc)
    call ddang (za%ddc, za%ddc(3:4), dd_logz%ddc(3:4))
    return
  end function

  function dd_log10 (qa)
    implicit none
    type (dd_real):: dd_log10
    type (dd_real), intent (in):: qa
    real (ddknd) qt1(4), qt2(4), qt3(4)
    call ddlog (qa%ddr, qt1)
    qt2(1) = 10.d0
    qt2(2) = 0.d0
    call ddlog (qt2, qt3)
    call dddiv (qt1, qt3, dd_log10%ddr)
    return
  end function

  function dd_log2 ()
    implicit none
    type (dd_real):: dd_log2
    call ddlog2c (dd_log2%ddr)
    return
  end function    

  function dd_maxq (qa, qb)
    implicit none
    type (dd_real):: dd_maxq
    type (dd_real), intent (in):: qa, qb
    integer ic
    call ddcpr (qa%ddr, qb%ddr, ic)
    if (ic >= 0) then
      call ddeq (qa%ddr, dd_maxq%ddr)
    else
      call ddeq (qb%ddr, dd_maxq%ddr)
    endif
    return
  end function

  function dd_maxq3 (qa, qb, qc)
    implicit none
    type (dd_real):: dd_maxq3
    type (dd_real), intent (in):: qa, qb, qc
    integer ic
    real (ddknd) qt0(4)
    call ddcpr (qa%ddr, qb%ddr, ic)
    if (ic >= 0) then
      call ddeq (qa%ddr, qt0)
    else
      call ddeq (qb%ddr, qt0)
    endif
    call ddcpr (qt0, qc%ddr, ic)
    if (ic >= 0) then
      call ddeq (qt0, dd_maxq3%ddr)
    else
      call ddeq (qc%ddr, dd_maxq3%ddr)
    endif
    return
  end function

  function dd_minq (qa, qb)
    implicit none
    type (dd_real):: dd_minq
    type (dd_real), intent (in):: qa, qb
    integer ic
    call ddcpr (qa%ddr, qb%ddr, ic)
    if (ic < 0) then
      call ddeq (qa%ddr, dd_minq%ddr)
    else
      call ddeq (qb%ddr, dd_minq%ddr)
    endif
    return
  end function

  function dd_minq3 (qa, qb, qc)
    implicit none
    type (dd_real):: dd_minq3
    type (dd_real), intent (in):: qa, qb, qc
    integer ic
    real (ddknd) qt0(4)
    call ddcpr (qa%ddr, qb%ddr, ic)
    if (ic < 0) then
      call ddeq (qa%ddr, qt0)
    else
      call ddeq (qb%ddr, qt0)
    endif
    call ddcpr (qt0, qc%ddr, ic)
    if (ic < 0) then
      call ddeq (qt0, dd_minq3%ddr)
    else
      call ddeq (qc%ddr, dd_minq3%ddr)
    endif
    return
  end function

  function dd_modq (qa, qb)
    implicit none
    type (dd_real):: dd_modq
    type (dd_real), intent (in):: qa, qb
    real (ddknd) qt1(4), qt2(4), qt3(4)
    call dddiv (qa%ddr, qb%ddr, qt1)
    call ddinfr (qt1, qt2, qt3)
    call ddmul (qb%ddr, qt2, qt1)
    call ddsub (qa%ddr, qt1, dd_modq%ddr)
    return
  end function

  function dd_nrt (qa, ib)
    implicit none
    type (dd_real):: dd_nrt
    type (dd_real), intent (in):: qa
    integer, intent (in):: ib
    call ddnrtf (qa%ddr, ib, dd_nrt%ddr)
    return
  end function

  subroutine dd_outq (iu, n1, n2, q1, q2, q3, q4, q5)
    implicit none
    integer, intent (in):: iu, n1, n2
    type (dd_real), intent (in):: q1, q2, q3, q4, q5
    optional:: q2, q3, q4, q5
    call ddout (iu, n1, n2, q1%ddr)
    if (present (q2)) call ddout (iu, n1, n2, q2%ddr)
    if (present (q3)) call ddout (iu, n1, n2, q3%ddr)
    if (present (q4)) call ddout (iu, n1, n2, q4%ddr)
    if (present (q5)) call ddout (iu, n1, n2, q5%ddr)
     return
  end subroutine

  subroutine dd_outz (iu, n1, n2, z1, z2, z3, z4, z5)
    implicit none
    integer, intent (in):: iu, n1, n2
    type (dd_complex), intent (in):: z1, z2, z3, z4, z5
    optional:: z2, z3, z4, z5
    call ddout (iu, n1, n2, z1%ddc)
    call ddout (iu, n1, n2, z1%ddc(3:4))
    if (present (z2)) call ddout (iu, n1, n2, z2%ddc)
    if (present (z2)) call ddout (iu, n1, n2, z2%ddc(3:4))
    if (present (z3)) call ddout (iu, n1, n2, z3%ddc)
    if (present (z3)) call ddout (iu, n1, n2, z3%ddc(3:4))
    if (present (z4)) call ddout (iu, n1, n2, z4%ddc)
    if (present (z4)) call ddout (iu, n1, n2, z4%ddc(3:4))
    if (present (z5)) call ddout (iu, n1, n2, z5%ddc)
    if (present (z5)) call ddout (iu, n1, n2, z5%ddc(3:4))
    return
  end subroutine

  function dd_pi ()
    implicit none
    type (dd_real):: dd_pi
    call ddpic (dd_pi%ddr)
    return
  end function    

  subroutine dd_poly (ia, qa, qb, qc)
    implicit none
    integer, intent (in):: ia
    type (dd_real), intent (in) :: qa(0:ia), qb
    type (dd_real), intent (out):: qc
    call ddpolyr (ia, qa(0)%ddr, qb%ddr, qc%ddr)
    return
  end subroutine

  function dd_polygamma (nn, ra)
    implicit none
    integer, intent (in):: nn
    type (dd_real), intent (in):: ra
    type (dd_real) dd_polygamma
    call ddpolygamma (nn, ra%ddr, dd_polygamma%ddr)
    return
  end function

  function dd_polygamma_be (nb, rb, nn, ra)
    implicit none
    integer, intent (in):: nb, nn
    type (dd_real), intent (in):: ra, rb(nb)
    type (dd_real) dd_polygamma_be
    call ddpolygammabe (nb, rb(1)%ddr, nn, ra%ddr, dd_polygamma_be%ddr)
    return
  end function

  subroutine dd_polylog_ini (nn, arr)
    implicit none
    integer, intent (in):: nn
    type (dd_real), intent (out):: arr(abs(nn))
    call ddpolylogini (nn, arr(1)%ddr)
    return
  end subroutine

  function dd_polylog_neg (nn, arr, ra)
    implicit none
    integer, intent (in):: nn
    type (dd_real), intent (in):: arr(abs(nn))
    type (dd_real), intent (in):: ra
    type (dd_real) dd_polylog_neg
    call ddpolylogneg (nn, arr(1)%ddr, ra%ddr, dd_polylog_neg%ddr)
    return
  end function

  function dd_polylog_pos (nn, ra)
    implicit none
    integer, intent (in):: nn
    type (dd_real), intent (in):: ra
    type (dd_real) dd_polylog_pos
    call ddpolylogpos (nn, ra%ddr, dd_polylog_pos%ddr)
    return
  end function

  function dd_qqtoz (qa, qb)
    implicit none
    type (dd_complex):: dd_qqtoz
    type (dd_real), intent (in):: qa, qb
    call ddqqc (qa%ddr, qb%ddr, dd_qqtoz%ddc)
    return
  end function

  function dd_qtod (qa)
    implicit none
    real (ddknd):: dd_qtod
    type (dd_real), intent (in):: qa
    dd_qtod = qa%ddr(1)
    return
  end function

  function dd_qtoq (qa)
    implicit none
    type (dd_real):: dd_qtoq
    type (dd_real), intent (in):: qa
    call ddeq (qa%ddr, dd_qtoq%ddr)
    return
  end function

  function dd_qtox (qa, qb)
    implicit none
    complex (ddknd):: dd_qtox
    type (dd_real), intent (in):: qa, qb
    real (ddknd) da, db
    da = qa%ddr(1)
    db = qb%ddr(1)
    dd_qtox = cmplx (da, db, ddknd)
    return
  end function

  function dd_qtoz (qa)
    implicit none
    type (dd_complex):: dd_qtoz
    type (dd_real), intent (in):: qa
    call ddmzc (qa%ddr, dd_qtoz%ddc)
    return
  end function

  function dd_signq (qa, qb)
    implicit none
    type (dd_real):: dd_signq
    type (dd_real), intent (in):: qa, qb
    call ddeq (qa%ddr, dd_signq%ddr)
    dd_signq%ddr(1) = sign (dd_signq%ddr(1), qb%ddr(1))
    if (qa%ddr(1) /= dd_signq%ddr(1)) dd_signq%ddr(2) = - dd_signq%ddr(2)
    return
  end function

  function dd_sin (qa)
    implicit none
    type (dd_real):: dd_sin
    type (dd_real), intent (in):: qa
    real (ddknd) qt1(4)
    call ddcssnr (qa%ddr, qt1, dd_sin%ddr)
    return
  end function

  function dd_sinz (za)
    implicit none
    type (dd_complex):: dd_sinz
    type (dd_complex), intent (in):: za
    real (ddknd) qt1(4), qt2(4), qt3(4), qt4(4), qt5(4), qt6(4)
    call ddeq (za%ddc(3:4), qt2)
    qt2(1) = - qt2(1)
    qt2(2) = - qt2(2)
    call ddexp (qt2, qt1)
    qt3(1) = 1.d0
    qt3(2) = 0.d0
    call dddiv (qt3, qt1, qt2)
    call ddcssnr (za%ddc, qt3, qt4)
    call ddadd (qt1, qt2, qt5)
    call ddmuld (qt5, 0.5d0, qt6)
    call ddmul (qt6, qt4, dd_sinz%ddc)
    call ddsub (qt1, qt2, qt5)
    call ddmuld (qt5, -0.5d0, qt6)
    call ddmul (qt6, qt3, dd_sinz%ddc(3:4))
    return
  end function

  function dd_sinh (qa)
    implicit none
    type (dd_real):: dd_sinh
    type (dd_real), intent (in):: qa
    real (ddknd) qt1(4)
    call ddcsshr (qa%ddr, qt1, dd_sinh%ddr)
    return
  end function

  function dd_sqrtq (qa)
    implicit none
    type (dd_real):: dd_sqrtq
    type (dd_real), intent (in):: qa
    call ddsqrt (qa%ddr, dd_sqrtq%ddr)
    return
  end function

  function dd_sqrtz (za)
    implicit none
    type (dd_complex):: dd_sqrtz
    type (dd_complex), intent (in):: za
    call ddcsqrt (za%ddc, dd_sqrtz%ddc)
    return
  end function

  function dd_struve_hn (nu, ra)
    implicit none
    integer, intent (in):: nu
    type (dd_real):: dd_struve_hn
    type (dd_real), intent (in):: ra
    call ddstruvehn (nu, ra%ddr, dd_struve_hn%ddr)
    return
  end function

  function dd_tan (qa)
    implicit none
    type (dd_real):: dd_tan
    type (dd_real), intent (in):: qa
    real (ddknd) qt1(4), qt2(4)
    call ddcssnr (qa%ddr, qt1, qt2)
    call dddiv (qt2, qt1, dd_tan%ddr)
    return
  end function

  function dd_tanh (qa)
    implicit none
    type (dd_real):: dd_tanh
    type (dd_real), intent (in):: qa
    real (ddknd) qt1(4), qt2(4)
    call ddcsshr (qa%ddr, qt1, qt2)
    call dddiv (qt2, qt1, dd_tanh%ddr)
    return
  end function

  function dd_xtoq (xa)
    implicit none
    type (dd_real):: dd_xtoq
    complex (ddknd), intent (in):: xa
    dd_xtoq%ddr(1) = xa
    dd_xtoq%ddr(2) = 0.d0
    return
  end function

  function dd_xtoz (xa)
    implicit none
    type (dd_complex):: dd_xtoz
    complex (ddknd), intent (in):: xa
    call ddxzc (xa, dd_xtoz%ddc)
    return
  end function

  function dd_zeta (ra)
    implicit none
    type (dd_real):: dd_zeta
    type (dd_real), intent (in):: ra
    call ddzetar (ra%ddr, dd_zeta%ddr)
    return
  end function

  function dd_zeta_be (nb, rb, rc)
    implicit none
    integer, intent (in):: nb
    type (dd_real):: dd_zeta_be
    type (dd_real), intent (in):: rb(nb), rc
    call ddzetabe (nb, rb(1)%ddr, rc%ddr, dd_zeta_be%ddr)
    return
  end function

  function dd_zeta_int (ia)
    implicit none
    type (dd_real):: dd_zeta_int
    integer, intent (in):: ia
    call ddzetaintr (ia, dd_zeta_int%ddr)
    return
  end function

  function dd_ztod (za)
    implicit none
    real (ddknd):: dd_ztod
    type (dd_complex), intent (in):: za
    dd_ztod = za%ddc(1)
    return
  end function

  function dd_ztoq (za)
    implicit none
    type (dd_real):: dd_ztoq
    type (dd_complex), intent (in):: za
    call ddeq (za%ddc, dd_ztoq%ddr)
    return
  end function

  function dd_ztox (za)
    implicit none
    complex (ddknd):: dd_ztox
    type (dd_complex), intent (in):: za
    real (ddknd) da, db
    da = za%ddc(1)
    db = za%ddc(3)
    dd_ztox = cmplx (da, db, ddknd)
    return
  end function

  function dd_ztoz (za)
    implicit none
    type (dd_complex):: dd_ztoz
    type (dd_complex), intent (in):: za
    call ddceq (za%ddc, dd_ztoz%ddc)
    return
  end function

end module ddmodule
