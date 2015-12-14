cmake_minimum_required(VERSION 2.8.5 FATAL_ERROR)
project("Situs")

set(CMAKE_CXX_FLAGS "-fpermissive")

add_subdirectory(fftw)
add_subdirectory(src)
