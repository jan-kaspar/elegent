cmake_minimum_required(VERSION 2.6)

# project name and version
project(Elegent)
set(Elegent_VERSION "1.3.0")

# search for installations of required packages
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${Elegent_SOURCE_DIR}/cmake/modules/")
find_package(ROOT REQUIRED)
find_package(GSL 2.0 REQUIRED)
find_package(HepMC)

message(STATUS "Elegent version: " ${Elegent_VERSION})

# configuration file
configure_file(
	"${PROJECT_SOURCE_DIR}/interface/Config.h.in"
	"${PROJECT_BINARY_DIR}/interface/Config.h"
)

include_directories("${PROJECT_BINARY_DIR}")

# compiler flags
add_definitions("-Wall -Wextra -Werror")
set(CMAKE_CXX_FLAGS "--std=c++1y -O3")

# list of subdirectories
add_subdirectory(interface)
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(tools)
add_subdirectory(test/functional)
