cmake_minimum_required(VERSION 2.6)

include_directories(${Elegent_SOURCE_DIR})
include_directories(${ROOT_INCLUDE_DIR})
include_directories(${HEPMC_INCLUDE})

set(LIBS Elegent ${HEPMC_LIB} ${ROOT_LIBRARIES} ${GSL_LIBRARIES})

add_executable(ElegentBDistributionSampler ElegentBDistributionSampler.cc)
target_link_libraries(ElegentBDistributionSampler ${LIBS})

add_executable(ElegentSDistributionSampler ElegentSDistributionSampler.cc)
target_link_libraries(ElegentSDistributionSampler ${LIBS})

add_executable(ElegentTDistributionSampler ElegentTDistributionSampler.cc)
target_link_libraries(ElegentTDistributionSampler ${LIBS})

add_executable(ElegentGeneratorTest ElegentGeneratorTest.cc)
target_link_libraries(ElegentGeneratorTest ${LIBS})

add_executable(debugger debugger.cc)
target_link_libraries(debugger ${LIBS})

install(TARGETS ElegentBDistributionSampler ElegentSDistributionSampler ElegentTDistributionSampler ElegentGeneratorTest DESTINATION bin)
