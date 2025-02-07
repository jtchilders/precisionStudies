module test_ddatanh
  use ddfuna  ! Import the module containing the ddatanh subroutine
  implicit none

contains

  subroutine unittest_ddatanh(filename)
    implicit none
    character(len=*), intent(in) :: filename
    real(8) :: hi_a, lo_a, expected_hi, expected_lo
    real(8) :: dda(2), ddc(2), expected_ddc(2)
    integer :: tolerance
    logical :: test_passed
    integer :: ios, total_tests, passed_tests
    integer :: scale_diff

    tolerance = 29  ! Double-double precision tolerance
    total_tests = 0
    passed_tests = 0

    ! Open the binary file
    open(unit=10, file=filename, form="unformatted", access="stream", status="old", action="read")
    write(*,*) "Running tests from:", filename

    do
      ! Read four double-precision values (one test case)
      read(10, iostat=ios) hi_a, lo_a, expected_hi, expected_lo
      if (ios /= 0) exit  ! Exit loop at end of file

      total_tests = total_tests + 1

      ! Pack inputs
      dda(1) = hi_a
      dda(2) = lo_a

      ! Expected result
      expected_ddc(1) = expected_hi
      expected_ddc(2) = expected_lo

      ! Call the ddatanh subroutine
      call ddatanh(dda, ddc)

      ! Calculate scale difference
      call dd_calc_scale_diff(ddc, expected_ddc, scale_diff)

      ! Compare results with expected values
      test_passed = scale_diff >= tolerance .or. scale_diff == 0

      ! Print results
      if (test_passed) then
        passed_tests = passed_tests + 1
      else
        write(*,*) "Test Failed:", &
                   "input: [", dda(1), ", ", dda(2), "]", &
                   "result: [", ddc(1), ", ", ddc(2), "]", &
                   "expected: [", expected_ddc(1), ", ", expected_ddc(2), "]", &
                   "error: [", abs(ddc(1) - expected_ddc(1)), ", ", abs(ddc(2) - expected_ddc(2)), "]", &
                   "scale_diff: [", scale_diff, "]"
      end if
    end do

    close(10)

    ! Summary
    write(*,*) "Tests completed:", passed_tests, "passed out of", total_tests

  end subroutine unittest_ddatanh

end module test_ddatanh