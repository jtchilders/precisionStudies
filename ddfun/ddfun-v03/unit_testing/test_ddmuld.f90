module test_ddmuld
   use ddfuna  ! Import the module containing the ddmuld subroutine
   implicit none
 
 contains
 
   subroutine unittest_ddmuld(filename)
     implicit none
     character(len=*), intent(in) :: filename
     real(8) :: double_value, dd_hi, dd_lo, expected_hi, expected_lo
     real(8) :: dda(2), result(2), expected_result(2)
     real(8) :: tolerance
     logical :: test_passed
     integer :: iunit, ios, total_tests, passed_tests
 
     tolerance = 1.0e-30  ! Double-double precision tolerance
     total_tests = 0
     passed_tests = 0
 
     ! Open the binary file
     open(unit=10, file=filename, form="unformatted", access="stream", status="old", action="read")
     write(*,*) "Running tests from:", filename
 
     do
       ! Read five double-precision values (one test case)
       read(10, iostat=ios) dd_hi, dd_lo, double_value, expected_hi, expected_lo
       if (ios /= 0) exit  ! Exit loop at end of file
 
       total_tests = total_tests + 1
 
       ! Pack inputs
       dda(1) = dd_hi
       dda(2) = dd_lo
 
       ! Expected result
       expected_result(1) = expected_hi
       expected_result(2) = expected_lo
 
       ! Call the ddmuld subroutine
       call ddmuld(dda, double_value, result)
 
       ! Compare results with expected values
       test_passed = abs(result(1) - expected_result(1)) < tolerance .and. &
                     abs(result(2) - expected_result(2)) < tolerance
 
       ! Print results
       if (test_passed) then
         passed_tests = passed_tests + 1
       else
         write(*,*) "Test Failed:", &
                    "inputs: double=", double_value, " dd=[", dda(1), ", ", dda(2), "]", &
                    "result: [", result(1), ", ", result(2), "]", &
                    "expected: [", expected_result(1), ", ", expected_result(2), "]", &
                    "error: [", abs(result(1) - expected_result(1)), ", ", abs(result(2) - expected_result(2)), "]"
       end if
     end do
 
     close(10)
 
     ! Summary
     write(*,*) "Tests completed:", passed_tests, "passed out of", total_tests
 
   end subroutine unittest_ddmuld
 
 end module test_ddmuld
 