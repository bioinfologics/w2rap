cmake_minimum_required(VERSION 3.13)
project(flash C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()


add_executable(flash
               combine_reads.c
               flash.c
               iostream.c
               read_io.c
               read_queue.c
               read_util.c
               util.c)

if (ZLIB_FOUND)
    target_link_libraries(flash ${ZLIB_LIBRARIES})
endif()

find_package (Threads REQUIRED)
target_link_libraries(flash ${CMAKE_THREAD_LIBS_INIT})
