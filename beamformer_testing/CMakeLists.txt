cmake_minimum_required(VERSION 3.10)
project(fold_vdif C)

if(APPLE AND EXISTS "/opt/homebrew")
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew")
endif()

find_library(FFTW3_LIB NAMES fftw3 REQUIRED)
find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h REQUIRED)

set(SOURCES
        fold_vdif.c
        plot.c
        stats.c
)

add_executable(fold_vdif ${SOURCES})

target_link_libraries(fold_vdif PRIVATE ${FFTW3_LIB} m)
target_include_directories(fold_vdif PRIVATE ${FFTW3_INCLUDE_DIR})
target_compile_options(fold_vdif PRIVATE -Wall -Wextra -O3)
