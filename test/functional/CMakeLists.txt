cmake_minimum_required(VERSION 2.6)

include_directories(${Elegent_SOURCE_DIR})
include_directories(${ROOT_INCLUDE_DIR})
include_directories(${HEPMC_INCLUDE})

set(LIBS Elegent ${HEPMC_LIB} ${ROOT_LIBRARIES} ${GSL_LIBRARIES})

add_executable(debugger debugger.cc)
target_link_libraries(debugger ${LIBS})
