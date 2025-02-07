module test_ddcdiv
  use ddfuna  ! Import the module containing the ddcdiv subroutine
  implicit none

contains

  subroutine unittest_ddcdiv(filename)
   implicit none
   character(len=*), intent(in) :: filename
   real(8) :: a_real_hi, a_real_lo, a_imag_hi, a_imag_lo
   real(8) :: b_real_hi, b_real_lo, b_imag_hi, b_imag_lo
   real(8) :: expected_real_hi, expected_real_lo, expected_imag_hi, expected_imag_lo
   real(8) :: a(4), b(4), c(4), expected_c(4)
   integer :: tolerance
   logical :: test_passed
   integer :: iunit, ios, total_tests, passed_tests
   integer :: scale_diff_real, scale_diff_imag

   tolerance = 29  ! Double-double precision tolerance
   total_tests = 0
   passed_tests = 0

   ! Open the binary file
   open(unit=10, file=filename, form="unformatted", access="stream", status="old", action="read")
   write(*,*) "Running tests from:", filename

   do
      ! Read eight double-precision values (one test case)
      read(10, iostat=ios) a_real_hi, a_real_lo, a_imag_hi, a_imag_lo, &
                           b_real_hi, b_real_lo, b_imag_hi, b_imag_lo, &
                           expected_real_hi, expected_real_lo, expected_imag_hi, expected_imag_lo
      if (ios /= 0) exit  ! Exit loop at end of file

      total_tests = total_tests + 1

      ! Pack inputs
      a(1) = a_real_hi
      a(2) = a_real_lo
      a(3) = a_imag_hi
      a(4) = a_imag_lo

      b(1) = b_real_hi
      b(2) = b_real_lo
      b(3) = b_imag_hi
      b(4) = b_imag_lo

      expected_c(1) = expected_real_hi
      expected_c(2) = expected_real_lo
      expected_c(3) = expected_imag_hi
      expected_c(4) = expected_imag_lo

      ! Call the ddcdiv subroutine
      call ddcdiv(a, b, c)

      ! Calculate scale difference
      call dd_calc_scale_diff(c, expected_c, scale_diff_real)
      call dd_calc_scale_diff(c(3), expected_c(3), scale_diff_imag)

      ! Compare results with expected values
      test_passed = (scale_diff_real >= tolerance .or. scale_diff_real == 0) .and. & 
                    (scale_diff_imag >= tolerance .or. scale_diff_imag == 0)

      ! Print results
      if (test_passed) then
         passed_tests = passed_tests + 1
      else
         write(*,*) "Test Failed:", &
                    "inputs: A [", a(1), ", ", a(2), "] / B [", b(1), ", ", b(2), "]", &
                    "result: [", c(1), ", ", c(2), "]", &
                    "expected: [", expected_c(1), ", ", expected_c(2), "]", &
                    "error: [", abs(c(1) - expected_c(1)), ", ", abs(c(2) - expected_c(2)), "]", &
                    "scale_diff: [", scale_diff_real, ", ", scale_diff_imag, "]"
      end if
   end do

   close(10)

   ! Summary
   write(*,*) "Tests completed:", passed_tests, "passed out of", total_tests

  end subroutine unittest_ddcdiv

end module test_ddcdiv