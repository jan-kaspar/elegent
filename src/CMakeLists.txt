cmake_minimum_required(VERSION 2.6)

# project name and version
project(Elegent)
set(Elegent_VERSION "2.0pre")

# search for installations of required packages
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${Elegent_SOURCE_DIR}/../cmake/modules/")
find_package(ROOT REQUIRED)
find_package(HepMC REQUIRED)
find_package(GSL REQUIRED)

# include SVN revision in project version, if available
include(FindSubversion)
if (SUBVERSION_FOUND)
	Subversion_WC_INFO(${Elegent_SOURCE_DIR} SVN)
	set(Elegent_VERSION "${Elegent_VERSION} (revision ${SVN_WC_REVISION})")
endif()

message(STATUS "Elegent version: " ${Elegent_VERSION})

# configuration file
configure_file(
	"${PROJECT_SOURCE_DIR}/libElegent/interface/Config.h.in"
	"${PROJECT_BINARY_DIR}/libElegent/interface/Config.h"
)

include_directories("${PROJECT_BINARY_DIR}/libElegent")

# compiler flags
add_definitions("-Wall -Wextra -Werror")
set(CMAKE_CXX_FLAGS "-O3")

# list of subdirectories
add_subdirectory(libElegent)
add_subdirectory(utilities)
