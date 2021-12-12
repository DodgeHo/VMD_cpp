cmake_minimum_required (VERSION 3.8)

project ("VMD_CPP")

find_package(Eigen3)

include_directories($SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/eigen-3.4.0 )
include_directories($SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/eigen-3.4.0/unsupported )
add_executable (VMD_CPP "VMD.cpp" "VMD.h" "VMD_Utils.cpp")

add_custom_target(build_all_configs
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ALL_BUILD --config Debug
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ALL_BUILD --config RelWithDebInfo
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ALL_BUILD --config Release
)