set(CPP_UNIT_TEST_DEPS
   test_ddabs.h
   test_ddadd.h
   test_dddiv.h
   test_dddivd.h
   test_ddexp.h
   test_ddlog.h
   test_ddmul.h
   test_ddmuld.h
   test_ddmuldd.h
   test_ddneg.h
   test_ddnpwr.h
   test_ddpower.h
   test_ddsqrt.h
   test_ddsub.h
)


# compile main into a binary, depends on the unit tests
add_executable(unit_tests main.cpp ${CPP_UNIT_TEST_DEPS})
target_link_libraries(unit_tests ddmath)
