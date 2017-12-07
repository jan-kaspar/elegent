cmake_minimum_required(VERSION 2.6)

include_directories(${Elegent_SOURCE_DIR})
include_directories(${ROOT_INCLUDE_DIR})

if (HEPMC_FOUND)
	set(LIBS Elegent ${HEPMC_LIB} ${ROOT_LIBRARIES} ${GSL_LIBRARIES})
else()
	set(LIBS Elegent ${ROOT_LIBRARIES} ${GSL_LIBRARIES})
endif()

add_executable(example1 example1.cc)
target_link_libraries(example1 ${LIBS})
