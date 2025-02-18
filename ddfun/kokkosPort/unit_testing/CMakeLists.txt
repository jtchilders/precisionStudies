set(CPP_UNIT_TEST_DEPS
   test_ddabs.h
   test_ddacosh.h
   test_ddadd.h
   test_ddasinh.h
   test_ddatanh.h
   test_ddcadd.h
   test_ddcdiv.h
   test_ddcmul.h
   test_ddcpwr.h
   test_ddcsqrt.h
   test_ddcsshr.h
   test_ddcssnr.h
   test_ddcsub.h
   test_dddiv.h
   test_dddivd.h
   test_ddexp.h
   test_ddlog.h
   test_ddmul.h
   test_ddmuld.h
   test_ddmuldd.h
   test_ddneg.h
   test_ddnpwr.h
   test_ddpolyr.h
   test_ddpower.h
   test_ddsqrt.h
   test_ddsub.h
)



# add_library(unit_tests SHARED ${CPP_UNIT_TEST_DEPS} ${CPP_SRC})
# target_include_directories(unit_tests PUBLIC ${CMAKE_SOURCE_DIR}/src)
# target_include_directories(unit_tests PUBLIC ${CMAKE_SOURCE_DIR}/unit_testing)
# # add kokkos dependencies
# target_link_libraries(unit_tests PUBLIC Kokkos::kokkos)