module test_ddcssnr
  use ddfuna  ! Import the module containing the ddcssnr subroutine
  implicit none

contains

  subroutine unittest_ddcssnr(filename)
    implicit none
    character(len=*), intent(in) :: filename
    real(8) :: hi_a, lo_a, expected_cos_hi, expected_cos_lo, expected_sin_hi, expected_sin_lo
    real(8) :: ddr_a(2), ddr_cos(2), ddr_sin(2), expected_ddr_cos(2), expected_ddr_sin(2)
    integer :: tolerance
    logical :: test_passed
    integer :: iunit, ios, total_tests, passed_tests
    integer :: scale_diff_cos, scale_diff_sin

    tolerance = 28  ! Double-double precision tolerance
    total_tests = 0
    passed_tests = 0

    ! Open the binary file
    open(unit=10, file=filename, form="unformatted", access="stream", status="old", action="read")
    write(*,*) "Running tests from:", filename

    do
      ! Read six double-precision values (one test case)
      read(10, iostat=ios) hi_a, lo_a, expected_cos_hi, expected_cos_lo, expected_sin_hi, expected_sin_lo
      if (ios /= 0) exit  ! Exit loop at end of file

      total_tests = total_tests + 1

      ! Pack input
      ddr_a(1) = hi_a
      ddr_a(2) = lo_a

      ! Expected results
      expected_ddr_cos(1) = expected_cos_hi
      expected_ddr_cos(2) = expected_cos_lo
      expected_ddr_sin(1) = expected_sin_hi
      expected_ddr_sin(2) = expected_sin_lo

      ! Call the ddcssnr subroutine
      call ddcssnr(ddr_a, ddr_cos, ddr_sin)

      ! Calculate scale differences
      call dd_calc_scale_diff(ddr_cos, expected_ddr_cos, scale_diff_cos)
      call dd_calc_scale_diff(ddr_sin, expected_ddr_sin, scale_diff_sin)

      ! Compare results with expected values
      test_passed = (scale_diff_cos >= tolerance .or. scale_diff_cos == 0) .and. &
                    (scale_diff_sin >= tolerance .or. scale_diff_sin == 0)

      ! Print results
      if (test_passed) then
        passed_tests = passed_tests + 1
      else
        write(*,*) "Test Failed:", &
                   "input: [", ddr_a(1), ", ", ddr_a(2), "]", &
                   "result_cos: [", ddr_cos(1), ", ", ddr_cos(2), "]", &
                   "expected_cos: [", expected_ddr_cos(1), ", ", expected_ddr_cos(2), "]", &
                  !  "cos_error: [", abs(ddr_cos(1) - expected_ddr_cos(1)), ", ", abs(ddr_cos(2) - expected_ddr_cos(2)), "]", &
                   "scale_diff_cos: [", scale_diff_cos, "]", &
                   "result_sin: [", ddr_sin(1), ", ", ddr_sin(2), "]", &
                   "expected_sin: [", expected_ddr_sin(1), ", ", expected_ddr_sin(2), "]", &
                  !  "sin_error: [", abs(ddr_sin(1) - expected_ddr_sin(1)), ", ", abs(ddr_sin(2) - expected_ddr_sin(2)), "]", &
                   "scale_diff_sin: [", scale_diff_sin, "]"
      end if
    end do

    close(10)

    ! Summary
    write(*,*) "Tests completed:", passed_tests, "passed out of", total_tests

  end subroutine unittest_ddcssnr

end module test_ddcssnr