cmake_minimum_required (VERSION 3.31)

project (matrix_library LANGUAGES C CXX)

# Compiler requirements
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Find libraries
find_package(nlohmann_json REQUIRED)
# HDF5 / commented out bc was causing problems
# find_package(HDF5 REQUIRED)
# find_package(HighFive REQUIRED)
find_package(Armadillo REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(chemfiles REQUIRED)

add_library(common INTERFACE)

target_include_directories(common INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(common INTERFACE
    Eigen3::Eigen
    armadillo
    nlohmann_json::nlohmann_json
    chemfiles
    BLAS::BLAS
    LAPACK::LAPACK
)

add_compile_options(-Wa,--gsframe=no) # suppresses a useless warning on x86 machines

# executables

add_executable(main src/main.cpp)
target_link_libraries(main PRIVATE common)

# google test stuff

include(FetchContent)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)

FetchContent_MakeAvailable(googletest)

add_executable(my_tests test/test_matrix.cpp)
target_link_libraries(my_tests PRIVATE common gtest_main)

enable_testing()
add_test(NAME my_tests COMMAND my_tests)
