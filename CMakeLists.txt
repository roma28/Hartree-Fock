cmake_minimum_required(VERSION 3.22)
project(SCF C)

set(CMAKE_C_STANDARD 99)
find_package(GSL REQUIRED)

include_directories(${GSL_INCLUDE_DIRS})

add_executable(SCF src/main.c src/atom_basis.c src/atom.c src/primitive_integrals.c include/primitive_integrals.h)
target_link_libraries(SCF ${GSL_LIBRARIES})
