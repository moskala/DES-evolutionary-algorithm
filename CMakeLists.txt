# Mostly adapted from https://github.com/pabloariasal/modern-cmake-sample, but simplified.
cmake_minimum_required(VERSION 3.17)
project(DES-evolutionary-algorithm LANGUAGES C)

add_library(desopt src/lib.c) 
target_include_directories(desopt
    PUBLIC 
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

add_executable(desopttest bin/main.c bin/tests.c bin/functions.c)
target_link_libraries(desopttest desopt)

# From https://hsf-training.github.io/hsf-training-cmake-webpage/aio/index.html
find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
    target_link_libraries(desopttest ${MATH_LIBRARY})
endif()

target_compile_options(desopt PRIVATE
    $<$<OR:$<C_COMPILER_ID:Clang>,$<C_COMPILER_ID:AppleClang>,$<C_COMPILER_ID:GNU>>:
    -Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:
    /W4>)
target_compile_options(desopttest PRIVATE
    $<$<OR:$<C_COMPILER_ID:Clang>,$<C_COMPILER_ID:AppleClang>,$<C_COMPILER_ID:GNU>>:
    -Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:
    /W4>)
