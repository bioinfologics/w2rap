# Minimum cmake version
cmake_minimum_required (VERSION 2.8.2)

# Name of project and that it is C++ only.
project (w2rap_lmp_proc)

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)

SET(CMAKE_CXX_FLAGS "-Wall -std=c++11 -lz")
add_definitions(${CMAKE_CXX_FLAGS})

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/seqan/share/cmake/Modules/")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin")

add_executable (dedup_fastq src/dedup_fastq.cpp)
