module test_ddcadd
  use ddfuna  ! Import the module containing the ddcadd and ddadd subroutines
  implicit none

contains

  subroutine unittest_ddcadd(filename)
   implicit none
   character(len=*), intent(in) :: filename
   real(ddknd) :: hi_a, lo_a, extra1_a, extra2_a, hi_b, lo_b, extra1_b, extra2_b
   real(ddknd) :: expected_hi, expected_lo, expected_extra1, expected_extra2
   real(ddknd) :: dda(4), ddb(4), ddc(4), expected_ddc(4)
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
      ! Read ten double-precision values (one test case)
      read(10, iostat=ios) hi_a, lo_a, extra1_a, extra2_a, hi_b, lo_b, extra1_b, extra2_b, expected_hi, expected_lo
      if (ios /= 0) exit  ! Exit loop at end of file
      read(10, iostat=ios) expected_extra1, expected_extra2

      total_tests = total_tests + 1

      ! Pack inputs
      dda(1) = hi_a
      dda(2) = lo_a
      dda(3) = extra1_a
      dda(4) = extra2_a

      ddb(1) = hi_b
      ddb(2) = lo_b
      ddb(3) = extra1_b
      ddb(4) = extra2_b

      ! Expected result
      expected_ddc(1) = expected_hi
      expected_ddc(2) = expected_lo
      expected_ddc(3) = expected_extra1
      expected_ddc(4) = expected_extra2

      ! Call the ddcadd subroutine
      call ddcadd(dda, ddb, ddc)

      ! Calculate scale difference
      call dd_calc_scale_diff(ddc(1:2), expected_ddc(1:2), scale_diff)
   
      ! Compare results with expected values
      test_passed = scale_diff >= tolerance .or. scale_diff == 0

      ! Check additional components
      call dd_calc_scale_diff(ddc(3:4), expected_ddc(3:4), scale_diff)

      ! Another comparison with expected values
      test_passed = test_passed .and. (scale_diff >= tolerance .or. scale_diff == 0)

      ! Print results
      if (test_passed) then
         passed_tests = passed_tests + 1
      else
         write(*,*) "Test Failed:", &
                    "inputs: [", dda(1), ", ", dda(2), ", ", dda(3), ", ", dda(4), "] + [", &
                                ddb(1), ", ", ddb(2), ", ", ddb(3), ", ", ddb(4), "]", &
                    "result: [", ddc(1), ", ", ddc(2), ", ", ddc(3), ", ", ddc(4), "]", &
                    "expected: [", expected_ddc(1), ", ", expected_ddc(2), ", ", expected_ddc(3), &
                                ", ", expected_ddc(4), "]", &
                    "error: [", abs(ddc(1) - expected_ddc(1)), ", ", abs(ddc(2) - expected_ddc(2)), &
                             ", ", abs(ddc(3) - expected_ddc(3)), ", ", abs(ddc(4) - expected_ddc(4)), "]", &
                    "scale_diff: [", scale_diff, "]"
      end if
   end do

   close(10)

   ! Summary
   write(*,*) "Tests completed:", passed_tests, "passed out of", total_tests

  end subroutine unittest_ddcadd

end module test_ddcadd