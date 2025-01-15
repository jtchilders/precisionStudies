

# make list of test files
set(TEST_SOURCES
    tests/test_ddadd.f90
    tests/test_dddiv.f90
    tests/test_dddivd.f90
    tests/test_ddexp.f90
    tests/test_ddlog.f90
    tests/test_ddmul.f90
    tests/test_ddmuld.f90
    tests/test_ddmuldd.f90
    tests/test_ddnpwr.f90
    tests/test_ddpower.f90
    tests/test_ddsub.f90
    tests/test_ddsqrt.f90
    )


# link into a library
add_library(unit_tests ${TEST_SOURCES})
# build commands for the $TEST_SOURCES to object files need to include the ddfun modules
target_include_directories(unit_tests PRIVATE ${CMAKE_BINARY_DIR}/fortran)


# link into an executable
add_executable(ddfun_tests_main main.f90)
target_link_libraries(ddfun_tests_main unit_tests)
target_link_libraries(ddfun_tests_main ddfun)

# Copy input files to the runtime directory
file(COPY 
     ${CMAKE_CURRENT_SOURCE_DIR}/python_generators/ddadd_test_cases.bin
     ${CMAKE_CURRENT_SOURCE_DIR}/python_generators/dddiv_test_cases.bin
     ${CMAKE_CURRENT_SOURCE_DIR}/python_generators/dddivd_test_cases.bin
     ${CMAKE_CURRENT_SOURCE_DIR}/python_generators/ddexp_test_cases.bin
     ${CMAKE_CURRENT_SOURCE_DIR}/python_generators/ddlog_test_cases.bin
     ${CMAKE_CURRENT_SOURCE_DIR}/python_generators/ddmul_test_cases.bin
     ${CMAKE_CURRENT_SOURCE_DIR}/python_generators/ddmuld_test_cases.bin
     ${CMAKE_CURRENT_SOURCE_DIR}/python_generators/ddmuldd_test_cases.bin
     ${CMAKE_CURRENT_SOURCE_DIR}/python_generators/ddnpwr_test_cases.bin
     ${CMAKE_CURRENT_SOURCE_DIR}/python_generators/ddpower_test_cases.bin
     ${CMAKE_CURRENT_SOURCE_DIR}/python_generators/ddsqrt_test_cases.bin
     ${CMAKE_CURRENT_SOURCE_DIR}/python_generators/ddsub_test_cases.bin
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
