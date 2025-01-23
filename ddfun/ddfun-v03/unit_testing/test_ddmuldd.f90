module test_ddmuldd
   use ddfuna  ! Import the module containing the ddmuldd subroutine
   implicit none
 
 contains
 
   subroutine unittest_ddmuldd(filename)
     implicit none
     character(len=*), intent(in) :: filename
     real(8) :: a, b, expected_hi, expected_lo
     real(8) :: result(2), expected_result(2)
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
       ! Read four double-precision values (one test case)
       read(10, iostat=ios) a, b, expected_hi, expected_lo
       if (ios /= 0) exit  ! Exit loop at end of file
 
       total_tests = total_tests + 1
 
       ! Expected result
       expected_result(1) = expected_hi
       expected_result(2) = expected_lo
 
       ! Call the ddmuldd subroutine
       call ddmuldd(a, b, result)
 
       ! Compare results with expected values
       test_passed = abs(result(1) - expected_result(1)) < tolerance .and. &
                     abs(result(2) - expected_result(2)) < tolerance
 
       ! Print results
       if (test_passed) then
         passed_tests = passed_tests + 1
       else
         write(*,*) "Test Failed:", &
                    "inputs: a=", a, " b=", b, &
                    "result: [", result(1), ", ", result(2), "]", &
                    "expected: [", expected_result(1), ", ", expected_result(2), "]", &
                    "error: [", abs(result(1) - expected_result(1)), ", ", abs(result(2) - expected_result(2)), "]"
       end if
     end do
 
     close(10)
 
     ! Summary
     write(*,*) "Tests completed:", passed_tests, "passed out of", total_tests
 
   end subroutine unittest_ddmuldd
 
 end module test_ddmuldd
 