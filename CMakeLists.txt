cmake_minimum_required(VERSION 3.17)
project(gen_risk2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

add_executable(gen_risk2 main.cpp GWAS.hxx)