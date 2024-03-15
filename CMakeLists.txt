cmake_minimum_required(VERSION 3.10)
project(fft)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(lib/flgl)

include_directories(${PROJECT_SOURCE_DIR}/lib/flgl/inc)
include_directories(${PROJECT_SOURCE_DIR}/lib/fftw/inc)
include_directories(${PROJECT_SOURCE_DIR}/lib/sw)

file(GLOB_RECURSE PROJECT_SOURCES "lib/sw/*.cpp" "src/*.cpp")

add_executable(${PROJECT_NAME} ${PROJECT_SOURCES})
set_target_properties(${PROJECT_NAME} PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

target_link_libraries(${PROJECT_NAME} FLGL)

find_library(FFTW
             NAMES fftw3
             HINTS lib/fftw/bin/
             NO_DEFAULT_PATH)

if(FFTW)
    target_link_libraries(${PROJECT_NAME} ${FFTW})
else()
    message(FATAL_ERROR "FFTW library not found")
endif()
