cmake_minimum_required (VERSION 3.31)

project (MatrixLibrary VERSION 1.0.0
  DESCRIPTION "Matrix Library with Linear Algebra Features"
  LANGUAGES CXX)

# Compiler requirements
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(GNUInstallDirs)

# Find libraries
find_package(nlohmann_json REQUIRED)
find_package(Armadillo REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Threads REQUIRED) # for google benchmarking
#find_package(chemfiles REQUIRED)

# Library
add_library(MatrixLibrary)

add_library(MatrixLibrary::MatrixLibrary ALIAS MatrixLibrary)

target_sources(
  MatrixLibrary
  PRIVATE src/matrix_utilfuncs.cpp
          src/matrix_basic_func.cpp
          src/matrix_operators.cpp
          src/matrix_eigendecomp.cpp
          src/helper_func.cpp
)

target_include_directories(
  MatrixLibrary
  PRIVATE src
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
			$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_link_libraries(MatrixLibrary PUBLIC
    Eigen3::Eigen
    armadillo
    nlohmann_json::nlohmann_json
    #chemfiles
    BLAS::BLAS
    LAPACK::LAPACK
)

add_compile_options(-Wa,--gsframe=no) # suppresses a useless warning on x86 machines

# executables
add_executable(main src/main.cpp) # need to remove for actual library when we're done testing
target_link_libraries(main PRIVATE MatrixLibrary)

# Google Test
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
FetchContent_MakeAvailable(googletest)

add_executable(my_tests test/test_matrix.cpp)
target_link_libraries(my_tests PRIVATE MatrixLibrary gtest_main)

enable_testing()
add_test(NAME my_tests COMMAND my_tests)

# Google Benchmark
FetchContent_Declare(
  googlebenchmark
  URL https://github.com/google/benchmark/archive/refs/tags/v1.9.0.zip
)
FetchContent_MakeAvailable(googlebenchmark)

add_executable(matrix_benchmarks
  benchmarking/benchmark_main.cpp
  benchmarking/benchmark_addition.cpp
  benchmarking/benchmark_multiplication.cpp
)

target_link_libraries(matrix_benchmarks 
  PRIVATE MatrixLibrary 
  benchmark::benchmark
  Threads::Threads)

################
# BEGIN INSTALL SECTION
#################

include(GNUInstallDirs)

install(TARGETS MatrixLibrary
				EXPORT MatrixLibraryTargets
				LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
				ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
				RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
				INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
	${CMAKE_CURRENT_SOURCE_DIR}/cmake/MatrixLibraryConfig.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/MatrixLibraryConfig.cmake
	INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/MatrixLibrary
)
	
write_basic_package_version_file(
	${CMAKE_CURRENT_BINARY_DIR}/MatrixLibraryConfigVersion.cmake
	VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(
	FILES ${CMAKE_CURRENT_BINARY_DIR}/MatrixLibraryConfig.cmake
				${CMAKE_CURRENT_BINARY_DIR}/MatrixLibraryConfigVersion.cmake
	DESTINATION ${CMAKE_INSTALLL_DATAROOTDIR}/cmake/MatrixLibrary
)