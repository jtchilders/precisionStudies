module test_ddpower
  use ddfuna  ! Import the module containing the ddpower subroutine
  implicit none

contains

  subroutine unittest_ddpower(filename)
    implicit none
    character(len=*), intent(in) :: filename
    real(8) :: hi_a, lo_a, hi_b, lo_b, expected_hi, expected_lo
    real(8) :: a(2), b(2), c(2), expected_c(2)
    real(8) :: tolerance
    logical :: test_passed
    integer :: iunit, ios, total_tests, passed_tests

    tolerance = 2.0e-29  ! Double-double precision tolerance
    total_tests = 0
    passed_tests = 0

    ! Open the binary file
    open(unit=10, file=filename, form="unformatted", access="stream", status="old", action="read")
    write(*,*) "Running tests from:", filename

    do
      ! Read six double-precision values (one test case)
      read(10, iostat=ios) hi_a, lo_a, hi_b, lo_b, expected_hi, expected_lo
      if (ios /= 0) exit  ! Exit loop at end of file

      total_tests = total_tests + 1

      ! Pack inputs
      a(1) = hi_a
      a(2) = lo_a
      b(1) = hi_b
      b(2) = lo_b

      ! Expected result
      expected_c(1) = expected_hi
      expected_c(2) = expected_lo

      ! Call the ddpower subroutine
      call ddpower(a, b, c)

      ! Compare results with expected values
      test_passed = abs(c(1) - expected_c(1)) < tolerance .and. &
                    abs(c(2) - expected_c(2)) < tolerance

      ! Print results
      if (test_passed) then
        passed_tests = passed_tests + 1
      else
        write(*,*) "Test Failed:", &
                   "input: [", a(1), ", ", a(2), "; ", b(1), ", ", b(2), "]", &
                   "result: [", c(1), ", ", c(2), "]", &
                   "expected: [", expected_c(1), ", ", expected_c(2), "]", &
                   "error: [", abs(c(1) - expected_c(1)), ", ", abs(c(2) - expected_c(2)), "]"
      end if
    end do

    close(10)

    ! Summary
    write(*,*) "Tests completed:", passed_tests, "passed out of", total_tests

  end subroutine unittest_ddpower

end module test_ddpower