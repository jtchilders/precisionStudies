module test_ddcmul
  use ddfuna  ! Import the module containing the ddcmul subroutine
  implicit none

contains

  subroutine unittest_ddcmul(filename)
    implicit none
    character(len=*), intent(in) :: filename
    real(8) :: a(4), b(4), expected_c(4)
    real(8) :: result_c(4)
    integer :: tolerance
    logical :: test_passed
    integer :: iunit, ios, total_tests, passed_tests
    integer :: scale_diff

    tolerance = 29  ! Double-double precision tolerance
    total_tests = 0
    passed_tests = 0

    ! Open the binary file
    open(unit=10, file=filename, form="unformatted", access="stream", status="old", action="read")
    write(*,*) "Running tests from:", filename

    do
      ! Read eight double-precision values (one test case)
      read(10, iostat=ios) a(1), a(2), a(3), a(4), b(1), b(2), b(3), b(4), &
                           expected_c(1), expected_c(2), expected_c(3), expected_c(4)
      if (ios /= 0) exit  ! Exit loop at end of file

      total_tests = total_tests + 1

      ! Call the ddcmul subroutine
      call ddcmul(a, b, result_c)

      ! Calculate scale difference
      call dd_calc_scale_diff(result_c, expected_c, scale_diff)

      ! Compare results with expected values
      test_passed = scale_diff >= tolerance .or. scale_diff == 0

      ! Print results
      if (test_passed) then
        passed_tests = passed_tests + 1
      else
        write(*,*) "Test Failed:", &
                   "inputs: A[", a(1), ", ", a(2), ", ", a(3), ", ", a(4), "] * B[", b(1), ", ", b(2), ", ", &
                              b(3), ", ", b(4), "]", &
                   "result: [", result_c(1), ", ", result_c(2), ", ", result_c(3), ", ", result_c(4), "]", &
                   "expected: [", expected_c(1), ", ", expected_c(2), ", ", expected_c(3), ", ", expected_c(4), "]", &
                   "error: [", abs(result_c(1) - expected_c(1)), ", ", abs(result_c(2) - expected_c(2)), ", ", &
                            abs(result_c(3) - expected_c(3)), ", ", abs(result_c(4) - expected_c(4)), "]", &
                   "scale_diff: [", scale_diff, "]"
      end if
    end do

    close(10)

    ! Summary
    write(*,*) "Tests completed:", passed_tests, "passed out of", total_tests

  end subroutine unittest_ddcmul

end module test_ddcmul