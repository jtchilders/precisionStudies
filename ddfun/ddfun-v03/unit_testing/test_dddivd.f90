module test_dddivd
   use ddfuna  ! Import the module containing the dddivd subroutine
   implicit none
 
 contains
 
   subroutine unittest_dddivd(filename)
     implicit none
     character(len=*), intent(in) :: filename
     real(8) :: dd_hi, dd_lo, divisor, expected_hi, expected_lo
     real(8) :: dd_value(2), result(2), expected_result(2)
     real(8) :: tolerance
     logical :: test_passed
     integer :: ios, total_tests, passed_tests
 
     tolerance = 1.0e-30  ! Double-double precision tolerance
     total_tests = 0
     passed_tests = 0
 
     ! Open the binary file
     open(unit=10, file=filename, form="unformatted", access="stream", status="old", action="read")
     write(*,*) "Running tests from:", filename
 
     do
       ! Read five double-precision values (one test case)
       read(10, iostat=ios) dd_hi, dd_lo, divisor, expected_hi, expected_lo
       if (ios /= 0) exit  ! Exit loop at end of file
 
       total_tests = total_tests + 1
 
       ! Pack inputs
       dd_value(1) = dd_hi
       dd_value(2) = dd_lo
 
       ! Expected result
       expected_result(1) = expected_hi
       expected_result(2) = expected_lo
 
       ! Call the dddivd subroutine
       call dddivd(dd_value, divisor, result)
 
       ! Compare results with expected values
       test_passed = abs(result(1) - expected_result(1)) < tolerance .and. &
                     abs(result(2) - expected_result(2)) < tolerance
 
       ! Print results
       if (test_passed) then
         passed_tests = passed_tests + 1
       else
         write(*,*) "Test Failed:", &
                    "inputs: dd=[", dd_value(1), ", ", dd_value(2), "] / divisor=", divisor, &
                    "result: [", result(1), ", ", result(2), "]", &
                    "expected: [", expected_result(1), ", ", expected_result(2), "]"
       end if
     end do
 
     close(10)
 
     ! Summary
     write(*,*) "Tests completed:", passed_tests, "passed out of", total_tests
 
   end subroutine unittest_dddivd
 
 end module test_dddivd
 