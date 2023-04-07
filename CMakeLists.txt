cmake_minimum_required(VERSION 3.16)
project(thesis)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_VERBOSE TRUE)

add_executable(main main.cpp utility_func.hpp utility_func.cpp)
add_executable(test unit_tests/test_runner.cpp unit_tests/tests.cpp utility_func.hpp utility_func.cpp)

target_include_directories(test PUBLIC "/opt/homebrew/Cellar/googletest/1.13.0/include")
target_link_libraries(test PUBLIC "/opt/homebrew/Cellar/googletest/1.13.0/lib/libgtest.a")