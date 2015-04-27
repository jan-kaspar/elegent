cmake_minimum_required(VERSION 2.6)

include_directories(${Elegent_SOURCE_DIR}/libElegent)
include_directories(${ROOT_INCLUDE_DIR})
include_directories(${HEPMC_INCLUDE})
include_directories(${GSL_INCLUDE_DIR})

add_library(Elegent SHARED Constants.cc Math.cc
	Model.cc BHModel.cc BSWModel.cc DLModel.cc ExpModel.cc FerreiraModel.cc GodizovModel.cc InterpolationModel.cc IslamModel.cc JenkovszkyModel.cc PPPModel.cc
	CoulombInterference.cc ModelFactory.cc Generator.cc)

install(TARGETS Elegent DESTINATION lib)
