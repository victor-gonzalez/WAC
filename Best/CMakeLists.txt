################################################################################################
# Project WAC/Best
################################################################################################

set(SOURCE_FILES
  BestAnalysisConfiguration.cpp
  BestConfiguration.cpp
  BestEventReader.cpp
  )

set(HEADER_FILES
  BestAnalysisConfiguration.hpp
  BestConfiguration.hpp
  BestEventReader.hpp
  )

ROOT_GENERATE_DICTIONARY(G__Best ${HEADER_FILES} LINKDEF BestLinkDef.h)

################################################################################################
# Create a shared library with geneated dictionary
################################################################################################

add_library(WacBest SHARED ${SOURCE_FILES} G__Best.cxx)
target_link_libraries(WacBest Base  ${ROOT_LIBRARIES} ${EXTRA_LIBS} )
target_include_directories(WacBest  PUBLIC Base ${EXTRA_INCLUDES} )

add_executable(TestBestAnalysisConfiguration TestBestAnalysisConfiguration.cpp)
target_include_directories(TestBestAnalysisConfiguration PUBLIC WacBest)
target_link_libraries(TestBestAnalysisConfiguration Base WacBest ${ROOT_LIBRARIES} ${CMAKE_DL_LIBS})

add_executable(RunBestSimulationSingleParticles RunBestSimulationSingleParticles.cpp)
target_include_directories(RunBestSimulationSingleParticles PUBLIC Base WacBest)
target_link_libraries(RunBestSimulationSingleParticles Base WacBest ${ROOT_LIBRARIES} ${CMAKE_DL_LIBS})

add_executable(RunBestSimulationTwoParticlesDiff RunBestSimulationTwoParticlesDiff.cpp)
target_include_directories(RunBestSimulationTwoParticlesDiff PUBLIC Base WacBest)
target_link_libraries(RunBestSimulationTwoParticlesDiff Base WacBest ${ROOT_LIBRARIES} ${CMAKE_DL_LIBS})

add_executable(statUncertainBest statUncertainBest.cpp)
target_include_directories(statUncertainBest PUBLIC Base WacBest ${EXTRA_INCLUDES})
target_link_libraries(statUncertainBest Base WacBest ${WAC_ROOT_LIBRARIES} ${CMAKE_DL_LIBS})

install(FILES  "${CMAKE_CURRENT_BINARY_DIR}/libBest.rootmap" "${CMAKE_CURRENT_BINARY_DIR}/libBest_rdict.pcm" DESTINATION "$ENV{WAC_LIB}")
install(TARGETS WacBest  LIBRARY DESTINATION "$ENV{WAC_LIB}") 
install(TARGETS TestBestAnalysisConfiguration DESTINATION "$ENV{WAC_BIN}")
install(TARGETS RunBestSimulationSingleParticles DESTINATION "$ENV{WAC_BIN}")
install(TARGETS RunBestSimulationTwoParticlesDiff DESTINATION "$ENV{WAC_BIN}")
install(TARGETS statUncertainBest DESTINATION "$ENV{WAC_BIN}")

