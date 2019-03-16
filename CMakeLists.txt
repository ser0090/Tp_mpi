cmake_minimum_required(VERSION 3.12)
project(Tp_mpi C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_COMPILER /usr/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/bin/mpic++)
include_directories(/usr/include/mpi/)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS}  -fopenmp")
# -O3 -mno-avx -march=native
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}  -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -fopenmp")
add_executable(Tp_mpi main.c)
add_executable(Tp_mpi_as main_asinc.c)