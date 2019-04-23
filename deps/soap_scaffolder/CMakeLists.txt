cmake_minimum_required(VERSION 3.13)
project(soap_scaffolder C)

set(CMAKE_C_STANDARD 11)
find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()

add_definitions(-DMER127)

include_directories(inc)

add_executable(s_prepare
               prepare.c
               check.c)

add_executable(s_map
               map.c
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
               readseq1by1.c)

add_executable(s_scaffold
               scaffold.c
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

if(ZLIB_FOUND)
    target_link_libraries(s_prepare ${ZLIB_LIBRARIES})
    target_link_libraries(s_map ${ZLIB_LIBRARIES})
    target_link_libraries(s_scaffold ${ZLIB_LIBRARIES})
endif()

find_package (Threads REQUIRED)
target_link_libraries(s_prepare ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(s_map ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(s_scaffold ${CMAKE_THREAD_LIBS_INIT})
