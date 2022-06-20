cmake_minimum_required(VERSION 3.20)
project(SCF C)

set(CMAKE_C_STANDARD 99)
find_package(GSL REQUIRED)

include_directories(${GSL_INCLUDE_DIRS})

add_executable(SCF src/main.c src/atom_basis.c src/atom.c src/primitive_integrals.c src/integrals.c src/integral_tools.c src/xyz_reader.c src/c-logger/logger.c src/c-logger/loggerconf.c)
target_link_libraries(SCF ${GSL_LIBRARIES})
