

set(CPP_SRC 
   ddmath.cpp
   ddmath.h
   ddouble.h)

add_library(ddmath STATIC ${CPP_SRC})

target_include_directories(ddmath PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})