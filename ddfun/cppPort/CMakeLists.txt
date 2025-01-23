# define project using C++
cmake_minimum_required(VERSION 3.10)
project(CPP_DDFUN LANGUAGES CXX)

# set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# add directories
add_subdirectory(src)
add_subdirectory(unit_testing)