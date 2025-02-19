

set(CPP_SRC
   ddouble.hpp
   ddcomplex.hpp)

set(FTN_SRC
   ddfuna_c.f90)

add_library(ddfuna SHARED ${FTN_SRC})
target_include_directories(ddfuna PUBLIC ${CMAKE_SOURCE_DIR}/src)
# # add kokkos dependencies
# target_link_libraries(ddmath PUBLIC Kokkos::kokkos)