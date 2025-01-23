module test_ddsub
  use ddfuna  ! Import the module containing the ddsub subroutine
  implicit none

contains

  subroutine unittest_ddsub(filename)
    implicit none
    character(len=*), intent(in) :: filename
    real(8) :: hi_a, lo_a, hi_b, lo_b, expected_hi, expected_lo
    real(8) :: dda(2), ddb(2), ddc(2), expected_ddc(2)
    real(8) :: tolerance
    logical :: test_passed
    integer :: ios, total_tests, passed_tests

    tolerance = 1.0e-30  ! Double-double precision tolerance
    total_tests = 0
    passed_tests = 0

    ! Open the binary file
    open(unit=10, file=filename, form="unformatted", access="stream", status="old", action="read")
    write(*,*) "Running tests for ddsub from:", filename

    do
      ! Read six double-precision values (one test case)
      read(10, iostat=ios) hi_a, lo_a, hi_b, lo_b, expected_hi, expected_lo
      if (ios /= 0) exit  ! Exit loop at end of file

      total_tests = total_tests + 1

      ! Pack inputs
      dda(1) = hi_a
      dda(2) = lo_a
      ddb(1) = hi_b
      ddb(2) = lo_b

      ! Expected result
      expected_ddc(1) = expected_hi
      expected_ddc(2) = expected_lo

      ! Call the ddsub subroutine
      call ddsub(dda, ddb, ddc)

      ! Compare results with expected values
      test_passed = abs(ddc(1) - expected_ddc(1)) < tolerance .and. &
                    abs(ddc(2) - expected_ddc(2)) < tolerance

      ! Print results
      if (test_passed) then
        passed_tests = passed_tests + 1
      else
        write(*,*) "Test Failed:", &
                   "inputs: [", dda(1), ", ", dda(2), "] - [", ddb(1), ", ", ddb(2), "]", &
                   "result: [", ddc(1), ", ", ddc(2), "]", &
                   "expected: [", expected_ddc(1), ", ", expected_ddc(2), "]", &
                   "error: [", abs(ddc(1) - expected_ddc(1)), ", ", abs(ddc(2) - expected_ddc(2)), "]"
      end if
    end do

    close(10)

    ! Summary
    write(*,*) "Tests completed:", passed_tests, "passed out of", total_tests

  end subroutine unittest_ddsub

end module test_ddsub
