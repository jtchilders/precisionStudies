# define project using C++
cmake_minimum_required(VERSION 3.10)
project(KOKKOS_DDFUN LANGUAGES CXX)

# set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# depends on Kokkos
find_package(Kokkos REQUIRED)

# add directories
add_subdirectory(src)
add_subdirectory(unit_testing)

# compile main into a binary, depends on the unit tests
add_executable(unit_tests unit_tests.cpp)
target_link_libraries(unit_tests PRIVATE Kokkos::kokkos)
target_include_directories(unit_tests PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_include_directories(unit_tests PRIVATE ${CMAKE_SOURCE_DIR}/unit_testing)