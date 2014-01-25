cmake_minimum_required(VERSION 2.6)

project(Elegent)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${Elegent_SOURCE_DIR}/../cmake/modules/")

find_package(ROOT REQUIRED)
find_package(HepMC REQUIRED)

add_subdirectory(libElegent)
add_subdirectory(utilities)
