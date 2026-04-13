cmake_minimum_required(VERSION 3.12)
project(cct)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

set(CMAKE_CXX_STANDARD 20)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if (COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif ()

add_library(cct_lib INTERFACE)

target_include_directories(cct_lib INTERFACE include sdsl-lite/include)
target_include_directories(cct_lib INTERFACE include bit_vector/include)
target_include_directories(cct_lib INTERFACE include bit_vector/_deps/pasta_utils-src/include)
target_include_directories(cct_lib INTERFACE include bit_vector/_deps/tlx-src)

target_include_directories(cct_lib INTERFACE include libsais)

set(LIBSAIS_DIR ${CMAKE_SOURCE_DIR}/libsais)
set(LIBSAIS_LIB ${LIBSAIS_DIR}/libsais.a)

add_library(libsais STATIC IMPORTED)
set_target_properties(libsais PROPERTIES IMPORTED_LOCATION ${LIBSAIS_LIB})

add_subdirectory(tests)

add_subdirectory(experiments)

add_subdirectory(googletest)
