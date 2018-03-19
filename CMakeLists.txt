cmake_minimum_required(VERSION 3.10)
project(ISOPE)

set(CMAKE_CXX_STANDARD 17)

add_executable(ISOPE main.cpp fft.hpp fft.cpp)
target_link_libraries(ISOPE fftw3)