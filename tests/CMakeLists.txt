cmake_minimum_required (VERSION 3.9)

set (CMAKE_CXX_COMPILER mpiicpc)

# Find required C++ libraries

find_package(GTest REQUIRED)
include(GoogleTest)

add_executable(testing
	tests.cpp
)


# Set flags
target_compile_options(testing PRIVATE -ax=AVX2 -Wall -O2 -mkl -std=c++14)

# Add linkers
target_include_directories(testing PRIVATE ${CMAKE_BINARY_DIR}/../include/)
target_link_libraries(testing GTest::GTest GTest::Main ${CMAKE_BINARY_DIR}/../lib/cspice.a)
gtest_discover_tests(testing)