module test_ddcsshr
  use ddfuna  ! Import the module containing the ddcsshr subroutine
  implicit none

contains

  subroutine unittest_ddcsshr(filename)
   implicit none
   character(len=*), intent(in) :: filename
   real(8) :: hi_a, lo_a, expected_hi_x, expected_lo_x, expected_hi_y, expected_lo_y
   real(8) :: a(2), x(2), y(2), expected_x(2), expected_y(2)
   integer :: tolerance
   logical :: test_passed
   integer :: ios, total_tests, passed_tests
   integer :: scale_diff_x, scale_diff_y

   tolerance = 29  ! Double-double precision tolerance
   total_tests = 0
   passed_tests = 0

   ! Open the binary file
   open(unit=10, file=filename, form="unformatted", access="stream", status="old", action="read")
   write(*,*) "Running tests from:", filename

   do
      ! Read six double-precision values (one test case)
      read(10, iostat=ios) hi_a, lo_a, expected_hi_x, expected_lo_x, expected_hi_y, expected_lo_y
      if (ios /= 0) exit  ! Exit loop at end of file

      total_tests = total_tests + 1

      ! Pack input
      a(1) = hi_a
      a(2) = lo_a

      ! Expected results
      expected_x(1) = expected_hi_x
      expected_x(2) = expected_lo_x
      expected_y(1) = expected_hi_y
      expected_y(2) = expected_lo_y

      ! Call the ddcsshr subroutine
      call ddcsshr(a, x, y)

      ! Calculate scale differences
      call dd_calc_scale_diff(x, expected_x, scale_diff_x)
      call dd_calc_scale_diff(y, expected_y, scale_diff_y)

      ! Compare results with expected values
      test_passed = (scale_diff_x >= tolerance .or. scale_diff_x == 0) .and. &
                    (scale_diff_y >= tolerance .or. scale_diff_y == 0)

      ! Print results
      if (test_passed) then
         passed_tests = passed_tests + 1
      else
         write(*,*) "Test Failed:", &
                    "input: [", a(1), ", ", a(2), "]", &
                    "result X: [", x(1), ", ", x(2), "]", &
                    "expected X: [", expected_x(1), ", ", expected_x(2), "]", &
                    "error X: [", abs(x(1) - expected_x(1)), ", ", abs(x(2) - expected_x(2)), "]", &
                    "scale_diff X: [", scale_diff_x, "]", &
                    "result Y: [", y(1), ", ", y(2), "]", &
                    "expected Y: [", expected_y(1), ", ", expected_y(2), "]", &
                    "error Y: [", abs(y(1) - expected_y(1)), ", ", abs(y(2) - expected_y(2)), "]", &
                    "scale_diff Y: [", scale_diff_y, "]"
      end if
   end do

   close(10)

   ! Summary
   write(*,*) "Tests completed:", passed_tests, "passed out of", total_tests

  end subroutine unittest_ddcsshr

end module test_ddcsshr