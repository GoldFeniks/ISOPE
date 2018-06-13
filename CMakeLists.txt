cmake_minimum_required(VERSION 3.10)
project(ISOPE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")

add_executable(ISOPE main.cpp fft.hpp fft.cpp model.hpp types.hpp utils.hpp)
target_link_libraries(ISOPE fftw3)