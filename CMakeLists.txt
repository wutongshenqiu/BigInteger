cmake_minimum_required(VERSION 3.17)
project(BigInteger)

set(CMAKE_CXX_STANDARD 20)

include_directories(.)

add_executable(BigInteger
        src/big_integer.cpp
        src/big_integer.hpp
        )

add_subdirectory(tests)