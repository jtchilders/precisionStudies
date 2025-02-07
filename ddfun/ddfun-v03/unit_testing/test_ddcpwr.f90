module test_ddcpwr
  use ddfuna  ! Import the module containing the ddcpwr subroutine
  implicit none

contains

  subroutine unittest_ddcpwr(filename)
    implicit none
    character(len=*), intent(in) :: filename
    real(ddknd) :: a(4), expected_b(4), b(4)
    integer :: n, tolerance, placeholder
    logical :: test_passed
    integer :: ios, total_tests, passed_tests
    integer :: scale_diff_real, scale_diff_imag

    tolerance = 30  ! Double-double precision tolerance
    total_tests = 0
    passed_tests = 0

    ! Open the binary file
    open(unit=10, file=filename, form="unformatted", access="stream", status="old", action="read")
    write(*,*) "Running tests from:", filename

    do
      ! Read eight double-precision values and one integer (one test case)
      read(10, iostat=ios) a(1), a(2), a(3), a(4), n, placeholder, expected_b(1), expected_b(2), expected_b(3), expected_b(4)
      if (ios /= 0) exit  ! Exit loop at end of file

      total_tests = total_tests + 1


      write(*,*) "Test Failed:", &
      "inputs: A = [", a(1), ",", a(2), ",", a(3), ",", a(4), "]", &
      "n = ", n, &
      "result: [", b(1), ",", b(2), ",", b(3), ",", b(4), "]", &
      "expected: [", expected_b(1), ",", expected_b(2), ",", expected_b(3), ",", expected_b(4), "]", &
      "scale_diff: [", scale_diff_real, ",", scale_diff_imag, "]"

      ! Call the ddcpwr subroutine
      call ddcpwr(a, n, b)

      ! Calculate scale difference
      call dd_calc_scale_diff(b, expected_b, scale_diff_real)
      call dd_calc_scale_diff(b(3), expected_b(3), scale_diff_imag)

      ! Compare results with expected values
      test_passed = (scale_diff_real >= tolerance .or. scale_diff_real == 0) .and. &
                    (scale_diff_imag >= tolerance .or. scale_diff_imag == 0)

      ! Print results
      if (test_passed) then
        passed_tests = passed_tests + 1
      else
        write(*,*) "Test Failed:", &
                   "inputs: A = [", a(1), ",", a(2), ",", a(3), ",", a(4), "]", &
                   "n = ", n, &
                   "result: [", b(1), ",", b(2), ",", b(3), ",", b(4), "]", &
                   "expected: [", expected_b(1), ",", expected_b(2), ",", expected_b(3), ",", expected_b(4), "]", &
                   "scale_diff: [", scale_diff_real, ",", scale_diff_imag, "]"
      end if
    end do

    close(10)

    ! Summary
    write(*,*) "Tests completed:", passed_tests, "passed out of", total_tests

  end subroutine unittest_ddcpwr

end module test_ddcpwr