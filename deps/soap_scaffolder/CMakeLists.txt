cmake_minimum_required (VERSION 2.8)

IF (MY_COMPILER_OPTION STREQUAL "Intel")
    FIND_PROGRAM(ICPC_PROGRAM icpc ...)
    IF(ICPC_PROGRAM)
        MESSAGE(STATUS "** Compiling with Intel settings **")
        IF(ENV{CC})
            MESSAGE(WARNING "Overwriting CC envvar")
        ENDIF()
        SET(ENV{CC} "${ICC_PROGRAM}")
    ELSE()
        MESSAGE(FATAL_ERROR "Intel compiler not found")
    ENDIF()
ENDIF ()


project (soap_scaffolder)

## Compiler flags
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fomit-frame-pointer -pthread -DMER127")
SET(CMAKE_C_FLAGS_DEBUG "-O0 -g")
SET(CMAKE_C_FLAGS_RELEASE "-O3 -w")

#SET(CMAKE_CXX_LINK_EXECUTABLE "/usr/bin/ld")
add_definitions(${CMAKE_C_FLAGS})

## To move the bins to the final location
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin")

#set(MALLOC_LIBRARY "" CACHE PATH "Malloc library path")

find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()

## Include directories
include_directories("inc/")

## Mac compilation issue
#include_directories("/usr/include/malloc")

## Set sources variables
add_library(prepare_lib OBJECT
        check.c)

add_library(map_lib OBJECT
        newhash.c
        kmer.c
        prlHashCtg.c
        seq.c
        check.c
        hashFunction.c
        prlRead2Ctg.c
        lib.c
        attachPEinfo.c
        connect.c
        mem_manager.c
        stack.c
        readseq1by1.c
        )

add_library(scaff_lib OBJECT
        orderContig.c
        fibHeap.c
        fib.c
        attachPEinfo.c
        check.c
        lib.c
        connect.c
        mem_manager.c
        stack.c
        output_contig.c
        loadPreGraph.c
        kmer.c
        arc.c
        darray.c
        loadGraph.c
        seq.c
        prlReadFillGap.c
        localAsm.c
        newhash.c
        searchPath.c
)
## add executables and respective sources
add_executable(s_prepare prepare.c $<TARGET_OBJECTS:prepare_lib>)

## add executables and respective sources
add_executable(s_map map.c $<TARGET_OBJECTS:map_lib>)

## add executables and respective sources
add_executable(s_scaff scaffold.c $<TARGET_OBJECTS:scaff_lib>)


##Zlib link
if (ZLIB_FOUND)
    set(ZLIB libz.so)
    target_link_libraries(s_prepare ${ZLIB_LIBRARIES})
    target_link_libraries(s_map ${ZLIB_LIBRARIES})
    target_link_libraries(s_scaff ${ZLIB_LIBRARIES})
endif()


