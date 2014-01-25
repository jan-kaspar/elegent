# try to find HepMC in user defined path
find_library(HEPMC_LIB
	NAMES
		HepMC
	PATHS
		${HEPMC_PREFIX}/lib
)

# if not try to find HepMC in standard instalation paths
if (${HEPMC_LIB} MATCHES "HEPMC_LIB-NOTFOUND")
	find_library(HEPMC_LIB
		NAMES
			HepMC
		PATHS
			/usr/lib
			/usr/local/lib
	)
endif()

if (NOT ${HEPMC_LIB} MATCHES "HEPMC_LIB-NOTFOUND")
	find_path(HEPMC_INCLUDE
		HepMC/GenEvent.h
			/usr/include
			/usr/local/include
			${HEPMC_PREFIX}/include
	)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HepMC DEFAULT_MSG
		HEPMC_LIB HEPMC_INCLUDE)
