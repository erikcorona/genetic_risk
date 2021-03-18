cmake_minimum_required(VERSION 3.17)
project(gen_risk2)


include_directories(gen_risk_lib)
add_subdirectory(gen_risk_lib)

add_subdirectory(Google_tests)

set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Wall -std=c++20")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -std=c++20 -O3 -march=native")

set_target_properties(gen_risk_lib PROPERTIES LINKER_LANGUAGE CXX)
add_executable(gen_risk main.cpp)
target_link_libraries(gen_risk gen_risk_lib)
