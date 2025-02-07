
# make list of test files
set(TEST_SOURCES
   test_ddabs.f90
   test_ddacosh.f90
   test_ddadd.f90
   test_ddasinh.f90
   test_ddatanh.f90
   test_ddcadd.f90
   test_ddcdiv.f90
   test_ddcmul.f90
   test_ddcpwr.f90
   test_ddcsqrt.f90
   test_ddcssnr.f90
   test_ddcsub.f90
   test_dddiv.f90
   test_dddivd.f90
   test_ddexp.f90
   test_ddlog.f90
   test_ddmul.f90
   test_ddmuld.f90
   test_ddmuldd.f90
   test_ddneg.f90
   test_ddnpwr.f90
   test_ddpolyr.f90
   test_ddpower.f90
   test_ddsqrt.f90
   test_ddsub.f90
)


# link into an executable
add_executable(unit_tests main.f90 ${TEST_SOURCES})
# target_link_libraries(unit_tests unit_test_lib)
target_link_libraries(unit_tests ddfun)
target_include_directories(unit_tests PRIVATE ${CMAKE_BINARY_DIR}/fortran)
