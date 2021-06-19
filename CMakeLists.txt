# Mostly adapted from https://github.com/pabloariasal/modern-cmake-sample, but simplified.
cmake_minimum_required(VERSION 3.17)
project(DES-evolutionary-algorithm LANGUAGES C)

add_library(lib src/lib.c)
target_include_directories(lib
    PUBLIC 
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

add_executable(bin bin/main.c)
target_link_libraries(bin lib m)

target_compile_options(lib PRIVATE
    $<$<OR:$<C_COMPILER_ID:Clang>,$<C_COMPILER_ID:AppleClang>,$<C_COMPILER_ID:GNU>>:
    -Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:
    /W4>)
target_compile_options(bin PRIVATE
    $<$<OR:$<C_COMPILER_ID:Clang>,$<C_COMPILER_ID:AppleClang>,$<C_COMPILER_ID:GNU>>:
    -Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:
    /W4>)
