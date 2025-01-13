
cmake_minimum_required(VERSION 3.10)
project(DDFUN_Tests LANGUAGES Fortran)


# Path to the directory containing `ddfun.mod`
set(DDFUN_MOD_DIR fortran) 
set(UNIT_TEST_DIR unit_testing)

# build libraries in subdirectories
add_subdirectory(${DDFUN_MOD_DIR})
add_subdirectory(${UNIT_TEST_DIR})



