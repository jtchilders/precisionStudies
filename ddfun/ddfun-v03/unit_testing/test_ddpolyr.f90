module test_ddpolyr
  use ddfuna  ! Import the module containing the ddpolyr subroutine and other ddsubroutines
  implicit none

contains

  subroutine unittest_ddpolyr(filename)
   implicit none
   character(len=*), intent(in) :: filename
   integer :: n, i, ios, total_tests, passed_tests, iunit, padding
   real(8) :: x0_hi, x0_lo, expected_x_hi, expected_x_lo
   real(8) :: x(2), x0(2), expected_x(2), a(2,0:10)  ! Assume maximum degree 10
   real(8) :: tolerance
   logical :: test_passed

   tolerance = 1.0e-30  ! Double-double precision tolerance
   total_tests = 0
   passed_tests = 0

   ! Open the binary file
   open(unit=10, file=filename, form="unformatted", access="stream", status="old", action="read")
   write(*,*) "Running tests from:", filename

   do
      ! Read n and the coefficients for the polynomial
      read(10, iostat=ios) n, padding
      if (ios /= 0) exit  ! Exit loop at end of file

      ! Read coefficients
      do i = 0, n
         read(10) a(1,i), a(2,i)
      end do

      ! Read initial guess x0 and expected result
      read(10) x0_hi, x0_lo, expected_x_hi, expected_x_lo
      read(10, iostat=ios)
      if (ios /= 0) exit  ! Exit loop at end of file

      total_tests = total_tests + 1

      ! Pack initial guess
      x0(1) = x0_hi
      x0(2) = x0_lo

      ! Expected result
      expected_x(1) = expected_x_hi
      expected_x(2) = expected_x_lo

      ! Call the ddpolyr subroutine
      call ddpolyr(n, a, x0, x)

      ! Compare results with expected values
      test_passed = abs(x(1) - expected_x(1)) < tolerance .and. &
                    abs(x(2) - expected_x(2)) < tolerance

      ! Print results
      if (test_passed) then
         passed_tests = passed_tests + 1
      else
         write(*,*) "Test Failed:", &
                    "n: ", n, &
                    "initial guess: [", x0(1), ", ", x0(2), "]", &
                    "result: [", x(1), ", ", x(2), "]", &
                    "expected: [", expected_x(1), ", ", expected_x(2), "]", &
                    "error: [", abs(x(1) - expected_x(1)), ", ", abs(x(2) - expected_x(2)), "]"
      end if
   end do

   close(10)

   ! Summary
   write(*,*) "Tests completed:", passed_tests, "passed out of", total_tests

  end subroutine unittest_ddpolyr

end module test_ddpolyr