################################################################################################
# Project WAC/Pythia
################################################################################################

ROOT_GENERATE_DICTIONARY(G__Pythia PythiaEventGenerator.hpp PythiaConfiguration.hpp LINKDEF PythiaLinkDef.h)

################################################################################################
# Create a shared library with geneated dictionary
################################################################################################
add_library(WacPythia SHARED PythiaEventGenerator.cpp PythiaConfiguration.cpp G__Pythia.cxx)
target_link_libraries(WacPythia Base  EGLib PythiaLib EGPythiaLib  ${ROOT_LIBRARIES} ${EXTRA_LIBS} )
target_include_directories(WacPythia  PUBLIC Base Pythia ${EXTRA_INCLUDES})

add_executable(RunPythiaSimulationNudyn RunPythiaSimulationNudyn.cpp)
target_include_directories(RunPythiaSimulationNudyn PUBLIC Base WacPythia ${EXTRA_INCLUDES})
target_link_libraries(RunPythiaSimulationNudyn Base WacPythia  EGLib PythiaLib EGPythiaLib   ${ROOT_LIBRARIES} ${EXTRA_LIBS} ${CMAKE_DL_LIBS})

add_executable(RunPythiaSimulationSingleParticles RunPythiaSimulationSingleParticles.cpp)
target_include_directories(RunPythiaSimulationSingleParticles PUBLIC Base WacPythia ${EXTRA_INCLUDES})
target_link_libraries(RunPythiaSimulationSingleParticles Base WacPythia Base EGLib PythiaLib EGPythiaLib   ${ROOT_LIBRARIES} ${EXTRA_LIBS} ${CMAKE_DL_LIBS})

add_executable(RunPythiaSimulationTwoParticles RunPythiaSimulationTwoParticles.cpp)
target_include_directories(RunPythiaSimulationTwoParticles PUBLIC Base WacPythia ${EXTRA_INCLUDES})
target_link_libraries(RunPythiaSimulationTwoParticles Base WacPythia Base EGLib PythiaLib EGPythiaLib   ${ROOT_LIBRARIES} ${EXTRA_LIBS} ${CMAKE_DL_LIBS})


install(FILES  "${CMAKE_CURRENT_BINARY_DIR}/libPythia.rootmap" "${CMAKE_CURRENT_BINARY_DIR}/libPythia_rdict.pcm" DESTINATION "$ENV{WAC_LIB}")
install(TARGETS WacPythia  LIBRARY DESTINATION "$ENV{WAC_LIB}") 
install(TARGETS RunPythiaSimulationSingleParticles DESTINATION "$ENV{WAC_BIN}") 
install(TARGETS RunPythiaSimulationTwoParticles DESTINATION "$ENV{WAC_BIN}") 
install(TARGETS RunPythiaSimulationNudyn DESTINATION "$ENV{WAC_BIN}") 

