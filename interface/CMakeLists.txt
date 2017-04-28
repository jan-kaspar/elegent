cmake_minimum_required(VERSION 2.6)

install(FILES
	BHModel.h
	BSWModel.h
	Constants.h
	CoulombInterference.h
	DLModel.h
	ExpModel.h
	FerreiraModel.h
	Generator.h
	GodizovModel.h
	InterpolationModel.h
	IFunctionInterpolator.h
	IslamModel.h
	IslamModel2015.h
	JenkovszkyModel.h
	Math.h
	ModelFactory.h
	Model.h
	PPPModel.h
	"${PROJECT_BINARY_DIR}/interface/Config.h"
		DESTINATION "include/Elegent")
