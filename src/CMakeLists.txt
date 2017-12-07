cmake_minimum_required(VERSION 2.6)

include_directories(${Elegent_SOURCE_DIR})
include_directories(${ROOT_INCLUDE_DIR})
include_directories(${GSL_INCLUDE_DIR})

if (HEPMC_FOUND)
	include_directories(${HEPMC_INCLUDE})

	add_library(Elegent SHARED Constants.cc Math.cc
		Model.cc BHModel.cc BSWModel.cc DLModel.cc ExpModel.cc FerreiraModel.cc GNLModel.cc
		GodizovModel.cc InterpolationModel.cc IFunctionInterpolator.cc IslamModel.cc IslamModel2015.cc
		JenkovszkyModel.cc PPPModel.cc CoulombInterference.cc ModelFactory.cc
		Generator.cc)
else()
	add_library(Elegent SHARED Constants.cc Math.cc
		Model.cc BHModel.cc BSWModel.cc DLModel.cc ExpModel.cc FerreiraModel.cc GNLModel.cc
		GodizovModel.cc InterpolationModel.cc IFunctionInterpolator.cc IslamModel.cc IslamModel2015.cc
		JenkovszkyModel.cc PPPModel.cc CoulombInterference.cc ModelFactory.cc)
endif()

install(TARGETS Elegent DESTINATION lib)
