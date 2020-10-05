cmake_minimum_required (VERSION 3.9)

project (SCROTAL LANGUAGES CXX)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set (CMAKE_CXX_COMPILER g++-8)

# Find required C++ libraries

#find_package(OpenMP REQUIRED)

# Add source files

add_executable(SCROTAL
    src/main.cpp
)

# Set flags

#target_compile_options(SCROTAL PRIVATE -qopenmp -ax=AVX2 -Wall -O2 -mkl -std=c++17)
target_compile_options(SCROTAL PRIVATE -fopenmp -Wall -O2 -std=c++17)

# Add linkers
target_include_directories(SCROTAL PRIVATE ${CMAKE_BINARY_DIR}/../include/)
target_include_directories(SCROTAL PRIVATE /user/local/include/)

target_link_libraries(SCROTAL PRIVATE ${CMAKE_BINARY_DIR}/../lib/cspice.a)
#target_link_libraries(SCROTAL PRIVATE -qopenmp)