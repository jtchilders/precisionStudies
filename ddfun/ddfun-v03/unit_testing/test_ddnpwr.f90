module test_ddnpwr
  use ddfuna  ! Import the module containing the ddnpwr subroutine
  implicit none

contains

  subroutine unittest_ddnpwr(filename)
    implicit none
    character(len=*), intent(in) :: filename
    real(8) :: hi_a, lo_a, expected_hi, expected_lo
    integer :: n, padding
    real(8) :: a(2), b(2), expected_b(2)
    integer :: tolerance
    logical :: test_passed
    integer :: iunit, ios, total_tests, passed_tests
    integer :: scale_diff

    tolerance = 30  ! Double-double precision tolerance
    total_tests = 0
    passed_tests = 0

    ! Open the binary file
    open(unit=10, file=filename, form="unformatted", access="stream", status="old", action="read")
    write(*,*) "Running tests from:", filename

    do
      ! Read six double-precision values and one integer (one test case)
      read(10, iostat=ios) hi_a, lo_a, n, padding, expected_hi, expected_lo
      if (ios /= 0) exit  ! Exit loop at end of file

      total_tests = total_tests + 1

      ! Pack inputs
      a(1) = hi_a
      a(2) = lo_a

      ! Expected result
      expected_b(1) = expected_hi
      expected_b(2) = expected_lo

      ! Call the ddnpwr subroutine
      call ddnpwr(a, n, b)

      ! Calculate scale difference
      call dd_calc_scale_diff(b, expected_b, scale_diff)

      ! Compare results with expected values
      test_passed = scale_diff >= tolerance .or. scale_diff == 0

      ! Print results
      if (test_passed) then
        passed_tests = passed_tests + 1
      else
        write(*,*) "Test Failed:", &
                   "input: [", a(1), ", ", a(2), "] ^ ", n, &
                   "result: [", b(1), ", ", b(2), "]", &
                   "expected: [", expected_b(1), ", ", expected_b(2), "]", &
                   "error: [", abs(b(1) - expected_b(1)), ", ", abs(b(2) - expected_b(2)), "]", &
                   "scale_diff: [", scale_diff, "]"
      end if
    end do

    close(10)

    ! Summary
    write(*,*) "Tests completed:", passed_tests, "passed out of", total_tests

  end subroutine unittest_ddnpwr

end module test_ddnpwr