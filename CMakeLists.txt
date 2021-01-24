cmake_minimum_required (VERSION 3.9)

project (ACROBAT LANGUAGES CXX)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set (CMAKE_CXX_COMPILER mpiicpc)

# Find required C++ libraries

enable_testing()
find_package(GTest REQUIRED)
include(GoogleTest)

# Add source files
add_executable(ACROBAT
    src/main.cpp
)

add_executable(testing
	tests/tests.cpp
)


# Set flags

target_compile_options(ACROBAT PRIVATE -qopenmp -ax=AVX2 -Wall -O2 -mkl -std=c++14)
target_compile_options(testing PRIVATE -ax=AVX2 -Wall -O2 -mkl -std=c++14)

# Add linkers
target_include_directories(ACROBAT PRIVATE ${CMAKE_BINARY_DIR}/../include/)
target_include_directories(testing PRIVATE ${CMAKE_BINARY_DIR}/../include/)

target_link_libraries(ACROBAT PRIVATE ${CMAKE_BINARY_DIR}/../lib/cspice.a)
target_link_libraries(ACROBAT PRIVATE -qopenmp)
target_link_libraries(testing GTest::GTest GTest::Main)
gtest_discover_tests(testing)