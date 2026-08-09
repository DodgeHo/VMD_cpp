cmake_minimum_required(VERSION 3.20)
project(VMD_cpp_package_consumer LANGUAGES CXX)

find_package(VMD_cpp 2 CONFIG REQUIRED)

add_executable(vmd_package_consumer main.cpp)
target_link_libraries(vmd_package_consumer PRIVATE VMD_cpp::vmd)
