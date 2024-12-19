
DDFUN: A double-double package with special functions

Revision date: 16 Mar 2023

AUTHOR:
  David H. Bailey
  Lawrence Berkeley National Lab (retired) and University of California, Davis
  Email: dhbailey@lbl.gov
 
COPYRIGHT AND DISCLAIMER

All software in this package (c) 2023 David H. Bailey. By downloading or using this software you agree to the copyright, disclaimer and license agreement in the accompanying file DISCLAIMER.txt.

INDEX OF THIS README FILE:

I. PURPOSE OF PACKAGE
II. INSTALLING CODING ENVIRONMENT (FOR MAC OS X SYSTEMS, IF NEEDED)
III. INSTALLING FORTRAN COMPILER (IF NEEDED)
IV. DOWNLOADING AND COMPILING DDFUN 
V. BRIEF SUMMARY OF CODING INSTRUCTIONS AND USAGE
VI. SAMPLE APPLICATION PROGRAMS AND TESTS
VII. RECENT UPDATES
VIII. APPENDIX: TRANSCENDENTAL, SPECIAL AND MISCELLANEOUS FUNCTIONS

+++++
I. PURPOSE OF PACKAGE

This package permits one to perform floating-point computations (real and complex) to double-double precision (approximately 28 digits), by making only relatively minor changes to existing Fortran programs.  It is intended as a substitute for IEEE 128-bit quad precision on systems where quad precision is not available. All basic arithmetic operations and transcendental functions are supported, together with numerous special functions. The package should run correctly on any Unix-based system supporting a Fortran-2008 compiler and IEEE 64-bit floating-point arithmetic. Note however that results are NOT guaranteed to the last bit.

In addition to fast execution times, one key feature of this package is a 100% THREAD-SAFE design, which means that user-level applications can be easily converted for parallel execution, say using a threaded parallel environment such as OpenMP.

Two related software packages by the same author are DQFUN (double-quad real and complex) and MPFUN-2020 (arbitrary precision real and complex). They are available from the same site as this package:

http://www.davidhbailey.com/dhbsoftware/


II. INSTALLING CODING ENVIRONMENT (FOR MAC OS X SYSTEMS, IF NEEDED)

For Apple Mac OS X systems (highly recommended for DDFUN), first install the latest supported version of Xcode, which is available for free from the Apple Developer website:

  https://developer.apple.com/

Here click on Account, then enter your Apple ID and password, then go to

  https://developer.apple.com/download/more/?=xcode

From this list, select Xcode. As of the above date, the latest version of Xcode is 13.3. Download this package on your system, double click to decompress, and place the resulting Xcode app somewhere on your system (typically in the Applications folder). Double-click on the app to run Xcode, allowing it install additional components, then quit Xcode. Open a terminal window, using the Terminal application in the Utilities folder, and type

  xcode-select --install

which installs various command-line tools. The entire process of downloading Xcode and installing command-line tools takes roughly 30 minutes. When this is completed, you should be ready to continue with the installation.


III. INSTALLING FORTRAN COMPILER (IF NEEDED)

Running DDFUN is relatively straightforward, provided that one has a Unix-based system, such as Linux or Apple OS X, and a Fortran-2008 compliant compiler that supports IEEE 64-bit arithmetic. These requirements are met by the GNU gfortran compiler, the Intel ifort compiler, and others.

The gfortran compiler (highly recommended for DDFUN) is available free for a variety of systems at this website:

  https://gcc.gnu.org/wiki/GFortranBinaries

For Apple Mac OS X systems, download the installer file here:

  https://github.com/fxcoudert/gfortran-for-macOS/releases

The gfortran compiler is normally placed in /usr/local/lib and /usr/local/bin. Thus before one uses gfortran, one must insert a line in one's shell initialization file (if the Z shell is used, as on most Apple OS X systems, the shell initialization file is ~/.zshrc). The line to be included is:

  PATH=/usr/local/lib:/usr/local/bin:$PATH

The following line is also recommended for gfortran compiler users:

  GFORTRAN_UNBUFFERED_ALL=yes; export GFORTRAN_UNBUFFERED_ALL

The following line is recommended for inclusion in the shell initialization file, no matter what compiler is used (it prevents stack overflow system errors):

  ulimit -s unlimited

On most Unix systems (including Apple Mac OS X systems), the shell initialization file must be manually executed upon initiating a terminal shell, typically by typing "source .zshrc".


IV. DOWNLOADING AND COMPILING DDFUN

From the website https://www.davidhbailey.com/dhbsoftware, download the file "ddfun-vnn.tar.gz" (replace "vnn" by whatever is the current version on the website, such as "v03"). If the file is not decompressed by your browser, use gunzip at the shell level to do this. Some browsers (such as the Apple Safari browser) do not drop the ".gz" suffix after decompression; if so, remove this suffix manually at the shell level. Then type

  tar xfv ddfun-vnn.tar 

(where again "vnn" is replaced by the downloaded version). This should create the directory and unpack all files. Compiling the library (in directory "fortran") is straightforward -- just compile the four library code files, ddfuna.f90, ddfune.f90, ddmodule.f90 and second.f90. Compile/link scripts are available for the GNU gfortran compiler. To compile the library with this script, type

  ./gnu-complib.scr

[Note: The GNU gfortran script to compile ddfun.f90 includes the flag "-fno-expensive-optimizations". This flag prevents the compiler from optimizing away key operations in several multiply and divide routines. A similar flag may be required with other compilers.]

Then to compile and link the application program tpslqm1dd.f90, using the GNU gfortran compiler, producing the executable file tpslqm1dd, type

  ./gnu-complink.scr tpslqm1dd

To execute the program, with output to tpslq1dd.txt, type

  ./tpslqm1dd > tpslqm1dd.txt

These scripts assume that the user program is in the same directory as the library files; this can easily be changed by editing the script files.

Several sample test programs, together with reference output files, are included in fortran directory -- see Section VI below.


V. BRIEF SUMMARY OF CODING INSTRUCTIONS AND USAGE

First place the following line in every subprogram of the user's application code that contains a double-double variable or array, at the beginning of the declaration section, before any implicit or type statements:

  use ddmodule

To designate a variable or array as double-double real (DDR) in your application code, use the Fortran-90 type statement with the type "dd_real", as in this example:

  type (dd_real) a, b(m), c(m,n)

Similarly, to designate a variable or array as double-double complex (DDC), use a type statement with "dd_complex". Thereafter when one of these variables or arrays appears in code, e.g.,

  d = a + b(i) * sqrt(3.d0 - c(i,j))

the proper underlying routines from are DDFUN package are automatically called by the Fortran compiler.

Most common mixed-mode combinations (arithmetic operations, comparisons and assignments) involving DDR, DDC, double precision (DP) and double complex (DC) variables are supported. Some conversion functions are provided between DP and DDR; and between DC and DDC -- see table below. Along this line, constants should be written with a d exponent, such as 3.25d0, to specify double precision; alternatively, one may write 3.25_ddknd.

However, there are some hazards. For example, the code r1 = 3.14159d0, where r1 is DDR, does NOT produce the true double-double equivalent of 3.14159, unless the numerical value is a modest-sized whole number or exact binary fraction. To obtain the full double-double converted value, write this as r1 = '3.14159' instead. Similarly, the code r2 = r1 + 3.d0 * sqrt (2.d0), where r1 and r2 are DDR, does NOT produce the true double-double value, since the expression 3.d0 * sqrt (2.d0) will be performed in double precision (according to standard Fortran-90 precedence rules). To obtain the fully accurate DDR result, write this as r2 = r1 + 3.d0 * sqrt (ddreal (2.d0)).

Input and output of DDR and DDC data are performed using the subroutines ddread and ddwrite. For example, to output the variable r1 in E format to Fortran unit 6 (standard output), to 30-digit accuracy, in a field of width 40 characters, use the line of code

  call ddwrite (6, 40, 30, r1)

The second argument (40 in the above example) must be at least 10 larger than the third argument (30 in the above example). To read the variable r1 from Fortran unit 5 (standard input), use the line of code

  call ddread (5, r1)

Most Fortran-2008 intrinsic functions are supported with DDR and DDC arguments, as appropriate, and numerous special functions are also supported. A complete list of supported functions and subroutines is summarized in the Appendix below (section VIII).


VI. SAMPLE APPLICATION PROGRAMS AND TESTS

Four application programs are included, together with corresponding reference output files (e.g., tpphixdd.ref.txt) for comparison with user results:

testddfun.f90  Tests most arithmetic and transcendental functions.

tpslqm1dd.f90  Performs the 1-level multipair PSLQ integer relation algorithm.

tpphixdd.f90  Performs a Poisson polynomial application, using 1-level multipair PSLQ.

tquaddd.f90   Evaluates definite integrals, using tanh-sinh, exp-sinh and sinh-sinh algorithms.

In addition, the package includes a test script that compiles the library and run each of the above sample programs above. To run this script, type:

  ./gnu-ddfun-tests.scr

For each test program, the script outputs either TEST PASSED or TEST FAILED. If all tests pass, then one can be fairly confident that the DDFUN software is working properly. 


VII. RECENT UPDATES

6 Jan 2023: This release (version v02) fixed problem with ddeformat and ddfformat; corrected several other minor bugs; added the acosh, asinh, atanh, ddeform and ddfform routines; changed the calling sequence for ddwrite to include width and field parameters; declared the subroutine names in ddfun.f90 as private in ddmodule.f90; inserted intent statements for all routines; changed all parameter statements to the object oriented style; changed all data statements to parameter statements; added several conversion functions; inserted test messages in the test programs and included a test script.

27 Feb 2023: Included a large set of special functions, by means of a special translator program to convert the file mpfune.f90 from MPFUN20-Fort.


VIII. APPENDIX: TRANSCENDENTAL, SPECIAL AND MISCELLANEOUS FUNCTIONS:

In these tables, "F" denotes function, "S" denotes subroutine, "DDR" denotes double-double real, "DDC" denotes double-double complex, "DP" denotes double precision, "DC" denotes double complex, and "Int" denotes integer. The variable names r1,r2,r3 are DDR; z1 is DDC; d1 is DP, dc1 is DC, i1,i2,i3,n,nb,np,nq are integers; s1 is character(1); sn is character(n); rb is an DDR vector of length nb; ss is an DDR vector of length n; aa is an DDR vector of length np; bb is an DDR vector of length nq.

1. Standard Fortran-2008 transcendental functions:

Type   Name              Description
DDR    abs(r1)           Absolute value
DDR    abs(z1)           Absolute value of complex arg
DDR    acos(r1)          Inverse cosine
DDR    acosh(r1)         Inverse hyperbolic cosine
DDR    aimag(z1)         Imaginary part of complex arg
DDR    aint(r1)          Truncates to integer
DDR    anint(r1)         Rounds to closest integer
DDR    asin(r1)          Inverse sine
DDR    asinh(r1)         Inverse hyperbolic sine
DDR    atan(r1)          Inverse tangent
DDR    atan2(r1,r2)      Arctangent with two args
DDR    atanh(r1)         Inverse hyperbolic tangent
DDR    bessel_j0(r1)     Bessel function of the first kind, order 0
DDR    bessel_j1(r1)     Bessel function of the first kind, order 1
DDR    bessel_jn(n,r1)   Besel function of the first kind, order n
DDR    bessel_y0(r1)     Bessel function of the second kind, order 0
DDR    bessel_y1(r1)     Bessel function of the second kind, order 1
DDR    bessel_yn(n,r1)   Besel function of the second kind, order n
DDC    conjg(z1)         Complex conjugate
DDR    cos(r1)           Cosine of real arg
DDC    cos(z1)           Cosine of complex arg
DDR    cosh(r1)          Hyperbolic cosine
DP     dble(r1)          Converts DDR argument to DP
DC     dcmplx(z1)        Converts DDC argument to DC
DDR    erf(r1)           Error function
DDR    erfc(r1)          Complementary error function
DDR    exp(r1)           Exponential function of real arg
DDC    exp(z1)           Exponential function of complex arg
DDR    gamma(r1)         Gamma function
DDR    hypot(r1,r2)      Hypotenuse of two args
DDR    log(r1)           Natural logarithm of real arg
DDC    log(z1)           Natural logarithm of complex arg
DDR    log10(r1)         Base-10 logarithm of real arg
DDR    max(r1,r2)        Maximum of two (or three) args
DDR    min(r1,r2)        Minimum of two (or three) args
DDR    mod(r1,r2)        Mod function = r1 - r2*aint(r1/r2)
DDR    sign(r1,r2)       Transfers sign from r2 to r1
DDR    sin(r1)           Sine function of real arg
DDC    sin(z1)           Sine function of complex arg
DDR    sinh(r1)          Hyperbolic sine
DDR    sqrt(r1)          Square root of real arg
DDC    sqrt(z1)          Square root of complex arg
DDR    tan(r1)           Tangent function
DDR    tanh(r1)          Hyperbolic tangent function

2. Special functions:

Type   Name                 Description 
hline
F:DDR  agm(r1,r2)           Arithmetic-geometric mean of r1 and f2
S      ddberne(nb,rb)       Initialize array rb of length nb with first nb even 
                              Bernoulli numbers [3] 
F:DDR  bessel_i(r1,r2)      BesselI function, order r1, of r2  
F:DDR  bessel_in(n,r1)      BesselI function, order n, of r1 
F:DDR  bessel_j(r1,r2)      BesselJ function, order r1, of r2  
F:DDR  bessel_jn(n,r1)      BesselJ function, order n, of r1 
F:DDR  bessel_i(r1,r2)      BesselK function, order r1, of r2  
F:DDR  bessel_kn(n,r1)      BesselK function, order n, of r1 
F:DDR  bessel_y(r1,r2)      BesselY function, order r1, of r2  
F:DDR  bessel_yn(n,r1)      BesselY function, order n, of r1 
F:DDR  digamma_be(nb,rb,r1) Digamma function of r1, using nb even Bernoulli 
                              numbers in rb [3] 
F:DDR  erf(r1)              Error function 
F:DDR  erfc(r1)             Complementary error function 
F:DDR  expint(r1)           Exponential integral function 
F:DDR  gamma(r1)            Gamma function 
F:DDR  hurwitz_zetan(k,r1)  Hurwitz zeta function, order n >= 2, of r1 [4] 
F:DDR  hurwitz_zetan_be (nb,rb,n,r1)  Hurwitz zeta function, order n >= 2, of r1, 
                              using nb even Bernoulli numbers in rb [3] 
F:DDR  hypergeom_pfq(np,nq,aa,bb,r1)  Hypergeometric pFq function of aa, bb and r1; 
                              dimensions are aa(np) and bb(nq) [5] 
F:DDR  incgamma(r1,r2)      Incomplete gamma function [6] 
F:DDR  polygamma(k,r1)      Polygamma function, order k >= 1, of r1 [4] 
F:DDR  polygamma(nb,rb,k,r1)  Polygamma function, order k >= 1, of r1, using 
                              nb even Bernoulli numbers in rb [3] 
S      polylog_ini(n,ss)    Initialize array ss, of size |n|, for computing
                              polylogarithms of order n when n < 0 [6]
F:DDR  polylog_neg(n,ss,r1) Polylogarithm function of r1, for n < 0, using 
                              precomputed data in ss [7]
F:DDR  polylog_pos(n,r1)    Polylogarithm function, order n >= 0, of r1 [8]
F:DDR  struve_hn(n,r1)      StruveH function, order n >= 0, of r1 [9]
F:DDR  zeta(r1)             Zeta function of r1
F:DDR  zeta_be(nb,rb,r1)    Zeta function of r1, using nb even Bernoulli
                              numbers in rb [3]
F:DDR  zeta_int(n)          Zeta function of integer argument n

Notes:
[3]: For most applications, set nb >= 140; see ddberne above.
[4]: For hurwitz_zetan and polygamma, the argument r1 is limited to the range (0, 1).
[5]: For hypergeom_pfq, the integers np and nq must not exceed 10.
[6]: For incgamma, r1 must not be zero, and must not be negative unless r2 = 0.
[7]: For polylog_ini and polylog_neg, the integer n is limited to the range [-1000, -1]. 
[8]: For polylog_pos, the argument r1 is limited to the range (-1,1). 
[9]: For struve_hn, the argument r1 is limited to the range [-1000, 1000].

3. Miscellaneous conversion, I/O and transcendental functions:

Type   Name                 Description
F:DP   dble(r1)             Converts DDR argument to DP
F:DC   dcmplx(z1)           Converts DDC argument to DC
F:DDC  ddcmplx(r1,r2)       Converts (r1,r2) to DDC
F:DDC  ddcmplx(dc1)         Converts DC arg to DDC
S      ddcssh(r1,r2,r3)     Returns both cosh and sinh of r1, in the same 
                              time as calling just cosh or just sinh
S      ddcssn(r1,r2,r3)     Returns both cos and sin of r1, in the same
                              time as calling just cos or just sin
S      dddecmd(r1,d1,i1)    Converts r1 to the form d1*10^i1
S      ddeform(r1,i1,i2,s1) Converts r1 to char(1) string in Ei1.i2
                              suitable for output
S      ddfform(r1,i1,i2,s1) Converts r1 to char(1) string in Fi1.i2
                              suitable for output
F:DDR  ddegamma()           Returns Euler's gamma constant
S      ddinit()             Initializes for extra-high precision
F:DDR  ddlog2()             Returns log(2)
F:DDR  ddnrt(r1,i1)         Returns the i1-th root of r1
F:DDR  ddpi()               Returns pi
F:DDR  ddrand(r1)           Returns pseudorandom number, based on r1
                              Start with an irrational, say r1 = ddlog2()
                              Typical iterated usage: r1 = ddrand(r1)
S      ddread(i1,r1)        Inputs r1 from Fortran unit i1; up to five
                              DDR args may be listed
S      ddread(i1,z1)        Inputs z1 from Fortran unit i1; up to five
                              DDC args may be listed
F:DDR  ddreal(i1)           Converts Int arg to DDR
F:DDR  ddreal(d1)           Converts DP arg to DDR
F:DDR  ddreal(z1)           Converts DDC arg to DDR
F:DDR  ddreal(s1,i1)        Converts character(1) string of length i1 to DDR
F:DDR  ddreal(sn)           Converts character(n) string to DDR
F:Int  ddwprec(r1)          Returns precision in words assigned to r1
F:Int  ddwprec(z1)          Returns precision in words assigned to z1
S      ddwrite(i1,i2,i3,r1) Outputs r1 in Ei2.i3 format to unit i1; up to
                              five DDR args may be listed
S      ddwrite(i1,i2,i3,z1) Outputs z1 in Ei2.i3 format to unit i1; up to
                              five DDC args may be listed


