cmake_minimum_required(VERSION 3.13)
project(nextclip C)

if("${CMAKE_C_COMPILER_ID}" MATCHES "[clang][Clang]")
    message(FATAL_ERROR "This code contains implementation details not supported by clang, please use GNU GCC or another alternative compatible with it")
endif()

set(CMAKE_C_STANDARD 11)

include_directories(include)

SET(MAXK 63 CACHE STRING "Maximum size of K to be used (default 64)")

SET(NUMBER_OF_BITFIELDS_IN_BINARY_KMER 1)
if(MAXK GREATER_EQUAL 63)
    SET(NUMBER_OF_BITFIELDS_IN_BINARY_KMER 2)
elseif(MAXK GREATER_EQUAL 95)
    SET(NUMBER_OF_BITFIELDS_IN_BINARY_KMER 3)
else(MAXK GREATER_EQUAL 127)
    SET(NUMBER_OF_BITFIELDS_IN_BINARY_KMER 4)
endif()

add_definitions(-DNBITS)

add_executable(nextclip
           src/binary_kmer.c
           src/element.c
           src/hash_table.c
           src/hash_value.c
           src/logger.c
           src/nextclip.c)

target_link_libraries(nextclip m)
