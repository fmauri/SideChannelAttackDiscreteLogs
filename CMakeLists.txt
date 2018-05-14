cmake_minimum_required(VERSION 3.10)
project(SideChannelAttackKangaroo)

set(CMAKE_CXX_STANDARD 11)

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

add_executable(SideChannelAttackKangaroo main.cpp DiscreteLog.h WeinerKangaroo/WeinerKangaroo.h WeinerKangaroo/WeinerKangaroo.cpp Chronometer.h)

target_link_libraries(SideChannelAttackKangaroo ntl m gmp)