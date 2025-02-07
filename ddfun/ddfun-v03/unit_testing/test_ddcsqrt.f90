module test_ddcsqrt
  use ddfuna  ! Import the module containing the ddcsqrt subroutine
  implicit none

contains

  subroutine unittest_ddcsqrt(filename)
   implicit none
   character(len=*), intent(in) :: filename
   real(8) :: a(4), expected_b(4), b(4)
   integer :: tolerance
   logical :: test_passed
   integer :: iunit, ios, total_tests, passed_tests
   integer :: scale_diff(2)

   tolerance = 30  ! Double-double precision tolerance
   total_tests = 0
   passed_tests = 0

   ! Open the binary file
   open(unit=10, file=filename, form="unformatted", access="stream", status="old", action="read")
   write(*,*) "Running tests from:", filename

   do
      ! Read eight double-precision values (one test case)
      read(10, iostat=ios) a(1), a(2), a(3), a(4), expected_b(1), expected_b(2), expected_b(3), expected_b(4)
      if (ios /= 0) exit  ! Exit loop at end of file

      total_tests = total_tests + 1

      ! Call the ddcsqrt subroutine
      call ddcsqrt(a, b)

      ! Calculate scale difference for real and imaginary parts
      call dd_calc_scale_diff(b(1:2), expected_b(1:2), scale_diff(1))
      call dd_calc_scale_diff(b(3:4), expected_b(3:4), scale_diff(2))

      ! Compare results with expected values
      test_passed = (scale_diff(1) >= tolerance .or. scale_diff(1) == 0) .and. &
                    (scale_diff(2) >= tolerance .or. scale_diff(2) == 0)

      ! Print results
      if (test_passed) then
         passed_tests = passed_tests + 1
      else
         write(*,*) "Test Failed:", &
                    "input: [", a(1), ", ", a(2), ", ", a(3), ", ", a(4), "]", &
                    "result: [", b(1), ", ", b(2), ", ", b(3), ", ", b(4), "]", &
                    "expected: [", expected_b(1), ", ", expected_b(2), ", ", expected_b(3), ", ", expected_b(4), "]", &
                    "error real: [", abs(b(1) - expected_b(1)), ", ", abs(b(2) - expected_b(2)), "]", &
                    "error imag: [", abs(b(3) - expected_b(3)), ", ", abs(b(4) - expected_b(4)), "]", &
                    "scale_diff real: [", scale_diff(1), "]", &
                    "scale_diff imag: [", scale_diff(2), "]"
      end if
   end do

   close(10)

   ! Summary
   write(*,*) "Tests completed:", passed_tests, "passed out of", total_tests

  end subroutine unittest_ddcsqrt

end module test_ddcsqrt