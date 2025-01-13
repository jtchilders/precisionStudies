program main
    use test_ddadd
    use test_ddsub
    use test_ddmul
    use test_dddiv
    implicit none
  
    write(*,*) "Starting unit tests for DDFUN library..."
  
    ! Call unit tests with respective input files
    call unittest_ddadd("ddadd_test_cases.bin")
    call unittest_ddsub("ddsub_test_cases.bin")
    call unittest_ddmul("ddmul_test_cases.bin")
    call unittest_dddiv("dddiv_test_cases.bin")
  
    write(*,*) "All tests completed."
  
  end program main
  