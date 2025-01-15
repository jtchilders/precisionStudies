program main
    use test_ddadd
    use test_dddiv
    use test_dddivd
    use test_ddexp
    use test_ddlog
    use test_ddmul
    use test_ddmuld
    use test_ddmuldd
    use test_ddnpwr
    use test_ddpower
    use test_ddsqrt
    use test_ddsub
    implicit none
  
    write(*,*) "Starting unit tests for DDFUN library..."
  
    ! Call unit tests with respective input files
    call unittest_ddadd("ddadd_test_cases.bin")
    call unittest_dddiv("dddiv_test_cases.bin")
    call unittest_dddivd("dddivd_test_cases.bin")
    call unittest_ddexp("ddexp_test_cases.bin")
    call unittest_ddlog("ddlog_test_cases.bin")
    call unittest_ddmul("ddmul_test_cases.bin")
    call unittest_ddmuld("ddmuld_test_cases.bin")
    call unittest_ddmuldd("ddmuldd_test_cases.bin")
    call unittest_ddnpwr("ddnpwr_test_cases.bin")
    call unittest_ddpower("ddpower_test_cases.bin")
    call unittest_ddsqrt("ddsqrt_test_cases.bin")
    call unittest_ddsub("ddsub_test_cases.bin")
  
    write(*,*) "All tests completed."
  
  end program main
  