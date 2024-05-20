################################################################################################
# Project WAC/Pythia
################################################################################################

set(SOURCE_FILES
  PythiaAnalysisConfiguration.cpp
  PythiaConfiguration.cpp
  PythiaEventGenerator.cpp
  )

set(HEADER_FILES
  PythiaAnalysisConfiguration.hpp
  PythiaConfiguration.hpp
  PythiaEventGenerator.hpp
  )

ROOT_GENERATE_DICTIONARY(G__Pythia ${HEADER_FILES} LINKDEF PythiaLinkDef.h)

################################################################################################
# Create a shared library with geneated dictionary
################################################################################################

add_library(WacPythia SHARED ${SOURCE_FILES} G__Pythia.cxx)
target_link_libraries(WacPythia Base  PYTHIA8LIB EGPythia8  ${ROOT_LIBRARIES} ${EXTRA_LIBS} )
target_include_directories(WacPythia  PUBLIC Base Pythia ${EXTRA_INCLUDES})

add_executable(TestPythiaAnalysisConfiguration TestPythiaAnalysisConfiguration.cpp)
target_include_directories(TestPythiaAnalysisConfiguration PUBLIC WacPythia ${EXTRA_INCLUDES})
target_link_libraries(TestPythiaAnalysisConfiguration Base WacPythia EG ${ROOT_LIBRARIES} ${EXTRA_LIBS} ${CMAKE_DL_LIBS})

add_executable(RunPythiaSimulationSingleParticles RunPythiaSimulationSingleParticles.cpp)
target_include_directories(RunPythiaSimulationSingleParticles PUBLIC Base WacPythia ${EXTRA_INCLUDES})
target_link_libraries(RunPythiaSimulationSingleParticles Base WacPythia Base EG PYTHIA8LIB EGPythia8 ${ROOT_LIBRARIES} ${EXTRA_LIBS} ${CMAKE_DL_LIBS})

add_executable(RunPythiaSimulationTwoParticles RunPythiaSimulationTwoParticles.cpp)
target_include_directories(RunPythiaSimulationTwoParticles PUBLIC Base WacPythia ${EXTRA_INCLUDES})
target_link_libraries(RunPythiaSimulationTwoParticles Base WacPythia Base EG PYTHIA8LIB EGPythia8 ${ROOT_LIBRARIES} ${EXTRA_LIBS} ${CMAKE_DL_LIBS})

add_executable(RunPythiaSimulationTwoParticlesDiff RunPythiaSimulationTwoParticlesDiff.cpp)
target_include_directories(RunPythiaSimulationTwoParticlesDiff PUBLIC Base WacPythia ${EXTRA_INCLUDES})
target_link_libraries(RunPythiaSimulationTwoParticlesDiff Base WacPythia Base EG PYTHIA8LIB EGPythia8 ${ROOT_LIBRARIES} ${EXTRA_LIBS} ${CMAKE_DL_LIBS})

add_executable(RunPythiaSimulationTwoParticlesDiffME RunPythiaSimulationTwoParticlesDiffME.cpp)
target_include_directories(RunPythiaSimulationTwoParticlesDiffME PUBLIC Base WacPythia ${EXTRA_INCLUDES})
target_link_libraries(RunPythiaSimulationTwoParticlesDiffME Base WacPythia Base EG PYTHIA8LIB EGPythia8 ${ROOT_LIBRARIES} ${EXTRA_LIBS} ${CMAKE_DL_LIBS})

add_executable(statUncertainPythia statUncertainPythia.cpp)
target_include_directories(statUncertainPythia PUBLIC Base WacPythia ${EXTRA_INCLUDES})
target_link_libraries(statUncertainPythia Base WacPythia EG ${WAC_ROOT_LIBRARIES} ${EXTRA_LIBS} ${CMAKE_DL_LIBS})

install(FILES  "${CMAKE_CURRENT_BINARY_DIR}/libPythia.rootmap" "${CMAKE_CURRENT_BINARY_DIR}/libPythia_rdict.pcm" DESTINATION "$ENV{WAC_LIB}")
install(TARGETS WacPythia  LIBRARY DESTINATION "$ENV{WAC_LIB}") 
install(TARGETS TestPythiaAnalysisConfiguration DESTINATION "$ENV{WAC_BIN}")
install(TARGETS RunPythiaSimulationSingleParticles DESTINATION "$ENV{WAC_BIN}")
install(TARGETS RunPythiaSimulationTwoParticles DESTINATION "$ENV{WAC_BIN}")
install(TARGETS RunPythiaSimulationTwoParticlesDiff DESTINATION "$ENV{WAC_BIN}")
install(TARGETS RunPythiaSimulationTwoParticlesDiffME DESTINATION "$ENV{WAC_BIN}")
install(TARGETS statUncertainPythia DESTINATION "$ENV{WAC_BIN}")

