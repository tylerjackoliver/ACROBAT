set(CMAKE_CXX_COMPILER mpiicpc)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set (CMAKE_CXX_COMPILER mpiicpc)

add_executable(ACROBAT
	main.cpp
)

# Set flags
target_compile_options(ACROBAT PRIVATE -qopenmp -O0 -Wall -g -mkl -std=c++14)

# Add linkers
target_include_directories(ACROBAT PRIVATE ${CMAKE_BINARY_DIR}/../include/)

target_link_libraries(ACROBAT PRIVATE ${CMAKE_BINARY_DIR}/../lib/cspice.a)
target_link_libraries(ACROBAT PRIVATE -qopenmp)