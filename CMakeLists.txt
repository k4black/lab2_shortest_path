cmake_minimum_required(VERSION 3.15)
project(shortest)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED)

add_executable(shortest main.cpp)

################################
# Testing
################################

add_executable(tests_graph tests/tests_graph.cpp)
add_executable(tests_auxiliary tests/tests_auxiliary.cpp)