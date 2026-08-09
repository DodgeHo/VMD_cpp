cmake_minimum_required(VERSION 3.20)

project(VMD_cpp VERSION 2.0.0 LANGUAGES CXX)

include(CMakePackageConfigHelpers)
include(CTest)
include(FetchContent)
include(GNUInstallDirs)

option(VMD_BUILD_CLI "Build the vmd_cli executable" ON)
option(VMD_BUILD_EXAMPLES "Build examples" ON)
option(VMD_BUILD_TESTS "Build tests" ${BUILD_TESTING})
option(VMD_BUILD_BENCHMARKS "Build benchmarks" ON)
option(VMD_BUILD_DOCS "Add the optional Doxygen docs target" ON)
option(VMD_ENABLE_WARNINGS "Enable high compiler warnings" ON)
option(VMD_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(VMD_USE_VENDORED_EIGEN "Use the repository's Eigen checkout before downloading Eigen" OFF)

set(VMD_FFT_BACKEND "Eigen" CACHE STRING "FFT backend: Eigen, FFTW, or MKL")
set_property(CACHE VMD_FFT_BACKEND PROPERTY STRINGS Eigen FFTW MKL)

find_package(Eigen3 3.3 CONFIG QUIET)
if(NOT TARGET Eigen3::Eigen AND VMD_USE_VENDORED_EIGEN AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Eigen/Eigen/Core")
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/Eigen"
        INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/Eigen")
endif()

if(NOT TARGET Eigen3::Eigen)
    FetchContent_Declare(
        eigen
        URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
        URL_HASH SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    if(POLICY CMP0169)
        cmake_policy(SET CMP0169 OLD)
    endif()
    FetchContent_GetProperties(eigen)
    if(NOT eigen_POPULATED)
        FetchContent_Populate(eigen)
    endif()
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${eigen_SOURCE_DIR}"
        INTERFACE_INCLUDE_DIRECTORIES "${eigen_SOURCE_DIR}")
endif()

add_library(vmd
    src/fft_backend.cpp
    src/vmd.cpp
    VMD_Utils.cpp
)
add_library(VMD_cpp::vmd ALIAS vmd)

target_compile_features(vmd PUBLIC cxx_std_17)
target_include_directories(vmd
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE src
)
target_link_libraries(vmd PUBLIC
    $<BUILD_INTERFACE:Eigen3::Eigen>
    $<INSTALL_INTERFACE:Eigen3::Eigen>
)

string(TOUPPER "${VMD_FFT_BACKEND}" VMD_FFT_BACKEND_UPPER)
if(VMD_FFT_BACKEND_UPPER STREQUAL "EIGEN")
    target_compile_definitions(vmd PRIVATE VMD_FFT_BACKEND_EIGEN=1)
elseif(VMD_FFT_BACKEND_UPPER STREQUAL "FFTW")
    find_path(FFTW3_INCLUDE_DIR fftw3.h REQUIRED)
    find_library(FFTW3_LIBRARY NAMES fftw3 libfftw3-3 REQUIRED)
    if(NOT TARGET FFTW3::fftw3)
        add_library(FFTW3::fftw3 UNKNOWN IMPORTED)
        set_target_properties(FFTW3::fftw3 PROPERTIES
            IMPORTED_LOCATION "${FFTW3_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}"
        )
    endif()
    target_link_libraries(vmd PRIVATE FFTW3::fftw3)
    target_compile_definitions(vmd PRIVATE VMD_FFT_BACKEND_FFTW=1)
elseif(VMD_FFT_BACKEND_UPPER STREQUAL "MKL")
    find_package(MKL CONFIG REQUIRED)
    target_link_libraries(vmd PRIVATE MKL::MKL)
    target_compile_definitions(vmd PRIVATE VMD_FFT_BACKEND_MKL=1)
else()
    message(FATAL_ERROR "VMD_FFT_BACKEND must be Eigen, FFTW, or MKL")
endif()

function(vmd_enable_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /utf-8)
    endif()
    if(NOT VMD_ENABLE_WARNINGS)
        return()
    endif()
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive-)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
    endif()
endfunction()
vmd_enable_warnings(vmd)

if(VMD_ENABLE_SANITIZERS AND NOT MSVC)
    target_compile_options(vmd PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(vmd PUBLIC -fsanitize=address,undefined)
endif()

if(VMD_BUILD_CLI)
    add_library(vmd_cli_support STATIC app/cli_io.cpp app/report.cpp)
    target_compile_features(vmd_cli_support PUBLIC cxx_std_17)
    target_include_directories(vmd_cli_support PUBLIC app)
    target_link_libraries(vmd_cli_support PUBLIC vmd)
    vmd_enable_warnings(vmd_cli_support)

    add_executable(vmd_cli app/vmd_cli.cpp)
    target_link_libraries(vmd_cli PRIVATE vmd_cli_support)
    vmd_enable_warnings(vmd_cli)
endif()

if(VMD_BUILD_EXAMPLES)
    add_executable(vmd_example examples/vmd_example.cpp)
    target_link_libraries(vmd_example PRIVATE vmd)
    vmd_enable_warnings(vmd_example)
endif()

if(VMD_BUILD_TESTS AND BUILD_TESTING)
    add_executable(vmd_tests tests/vmd_tests.cpp)
    target_link_libraries(vmd_tests PRIVATE vmd)
    vmd_enable_warnings(vmd_tests)
    add_test(NAME vmd.unit COMMAND vmd_tests)

    if(VMD_BUILD_CLI)
        add_test(
            NAME vmd.cli.e2e
            COMMAND ${CMAKE_COMMAND}
                -DCLI=$<TARGET_FILE:vmd_cli>
                -DOUTPUT_DIR=${CMAKE_CURRENT_BINARY_DIR}/cli-e2e-output
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cli_e2e.cmake
        )
    endif()
endif()

if(VMD_BUILD_BENCHMARKS)
    add_executable(vmd_benchmark benchmarks/vmd_benchmark.cpp)
    target_link_libraries(vmd_benchmark PRIVATE vmd)
    vmd_enable_warnings(vmd_benchmark)
endif()

if(VMD_BUILD_DOCS)
    find_package(Doxygen QUIET)
    if(Doxygen_FOUND)
        set(DOXYGEN_USE_MDFILE_AS_MAINPAGE README.md)
        doxygen_add_docs(docs README.md include/vmd docs COMMENT "Generate VMD_cpp API documentation")
    else()
        message(STATUS "Doxygen not found; the optional docs target will not be generated")
    endif()
endif()

install(TARGETS vmd EXPORT VMD_cppTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES VMD.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_package_config_file(
    cmake/VMD_cppConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/VMD_cppConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VMD_cpp
)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/VMD_cppConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)
install(EXPORT VMD_cppTargets
    FILE VMD_cppTargets.cmake
    NAMESPACE VMD_cpp::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VMD_cpp
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/VMD_cppConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/VMD_cppConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/VMD_cpp
)
