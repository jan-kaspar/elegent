find_program(ROOT_CONFIG_EXECUTABLE root-config
	PATHS ${ROOT_PREFIX}/bin $ENV{ROOTSYS}/bin)

if (NOT ROOT_CONFIG_EXECUTABLE)
	message(FATAL_ERROR "Could NOT find ROOT")
	set(ROOT_FOUND FALSE)
else ()    
	set(ROOT_FOUND TRUE)

	execute_process(
	  COMMAND ${ROOT_CONFIG_EXECUTABLE} --prefix 
	  OUTPUT_VARIABLE ROOTSYS 
	  OUTPUT_STRIP_TRAILING_WHITESPACE)

	execute_process(
	  COMMAND ${ROOT_CONFIG_EXECUTABLE} --version 
	  OUTPUT_VARIABLE ROOT_VERSION
	  OUTPUT_STRIP_TRAILING_WHITESPACE)

	execute_process(
	  COMMAND ${ROOT_CONFIG_EXECUTABLE} --incdir
	  OUTPUT_VARIABLE ROOT_INCLUDE_DIR
	  OUTPUT_STRIP_TRAILING_WHITESPACE)

	execute_process(
	  COMMAND ${ROOT_CONFIG_EXECUTABLE} --libs
	  OUTPUT_VARIABLE ROOT_LIBRARIES
	  OUTPUT_STRIP_TRAILING_WHITESPACE)

	message(STATUS "Found ROOT ${ROOT_VERSION} in ${ROOTSYS}")
	message(STATUS "\t ROOT_INCLUDE_DIR=${ROOT_INCLUDE_DIR}")
	message(STATUS "\t ROOT_LIBRARIES=${ROOT_LIBRARIES}")
endif()
