cmake_minimum_required(VERSION 2.6)

include_directories(${Elegent_SOURCE_DIR})
include_directories(${ROOT_INCLUDE_DIR})

if (HEPMC_FOUND)
	set(LIBS Elegent ${HEPMC_LIB} ${ROOT_LIBRARIES} ${GSL_LIBRARIES})
else()
	set(LIBS Elegent ${ROOT_LIBRARIES} ${GSL_LIBRARIES})
endif()

add_executable(debugger debugger.cc)
target_link_libraries(debugger ${LIBS})

add_executable(coulomb_I_function_interpolation coulomb_I_function_interpolation.cc)
target_link_libraries(coulomb_I_function_interpolation ${LIBS})

add_executable(coulomb_I_function_interpolation_t coulomb_I_function_interpolation_t.cc)
target_link_libraries(coulomb_I_function_interpolation_t ${LIBS})
