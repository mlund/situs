cmake_minimum_required(VERSION 2.8.5 FATAL_ERROR)
project("Situs")

set(CMAKE_CXX_FLAGS "-fpermissive")

option(ENABLE_OPENMP "Try to use OpenMP parallization" off)
if (ENABLE_OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

add_subdirectory(fftw)
add_subdirectory(src)
