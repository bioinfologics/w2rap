# Minimum cmake version
cmake_minimum_required (VERSION 2.8.2)

# Name of project and that it is C++ only.
project (w2rap_lmp_proc CXX)

find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()

set(CMAKE_CXX_STANDARD 11)

add_subdirectory(deps/flash)
add_subdirectory(deps/nextclip)
add_subdirectory(deps/soap_scaffolder)


set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)

SET(CMAKE_CXX_FLAGS "-Wall")
add_definitions(${CMAKE_CXX_FLAGS})

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/seqan/share/cmake/Modules/")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin")

add_executable (dedup_fastq src/dedup_fastq.cpp)

if(ZLIB_FOUND)
    target_link_libraries(dedup_fastq ${ZLIB_LIBRARIES})
endif()

install(TARGETS dedup_fastq s_scaffold s_map s_prepare nextclip flash
        RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/bin
        )