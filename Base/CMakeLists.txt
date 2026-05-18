################################################################################################
# Project WAC/Base
################################################################################################

#include(${ROOT_USE_FILE})
#include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
#add_definitions(${ROOT_CXX_FLAGS})

set(SOURCE_FILES
  TaskConfiguration.cpp
  Factory.cpp
  HistogramCollection.cpp
  Histograms.cpp
  Task.cpp
  TaskCollection.cpp
  MessageLogger.cpp
  AnalysisConfiguration.cpp
  Event.cpp
  EventPool.cpp
  EventFilter.cpp
  EventHistos.cpp
  EventLoop.cpp
  GeneralizedGaussianConfiguration.cpp
  GeneralizedGaussianProfile.cpp
  GeneratorConfiguration.cpp
  TwoPartCorrelationAnalyzer.cpp
  TwoPartDiffCorrelationAnalyzer.cpp
  TwoPartDiffCorrelationAnalyzerME.cpp
  Particle.cpp
  ParticleFilter.cpp
  ParticleHistos.cpp
  ParticleAnalyzer.cpp
  ParticlePairCombinedHistos.cpp
  ParticlePairBalanceFunctionDiffHistos.cpp
  ParticlePairCombinedDiffHistos.cpp
  ParticlePairDerivedHistos.cpp
  ParticlePairDerivedDiffHistos.cpp
  ParticlePairFilter.cpp
  ParticlePairHistos.cpp
  ParticlePairDiffHistos.cpp
  TrackAndPairConfiguration.cpp
  DetectorEffectsTask.cpp
  )

set(HEADER_FILES
  TaskConfiguration.hpp
  Factory.hpp
  HistogramCollection.hpp
  Histograms.hpp
  Task.hpp
  TaskCollection.hpp
  MessageLogger.hpp
  AnalysisConfiguration.hpp
  Event.hpp
  EventPool.hpp
  EventFilter.hpp
  EventHistos.hpp
  EventLoop.hpp
  GeneralizedGaussianConfiguration.hpp
  GeneralizedGaussianProfile.hpp
  GeneratorConfiguration.hpp
  TwoPartCorrelationAnalyzer.hpp
  TwoPartDiffCorrelationAnalyzer.hpp
  TwoPartDiffCorrelationAnalyzerME.hpp
  Particle.hpp
  ParticleFilter.hpp
  ParticleHistos.hpp
  ParticleAnalyzer.hpp
  ParticlePairCombinedHistos.hpp
  ParticlePairBalanceFunctionDiffHistos.hpp
  ParticlePairCombinedDiffHistos.hpp
  ParticlePairDerivedHistos.hpp
  ParticlePairDerivedDiffHistos.hpp
  ParticlePairFilter.hpp
  ParticlePairHistos.hpp
  ParticlePairDiffHistos.hpp
  TrackAndPairConfiguration.hpp
  DetectorEffectsTask.hpp
  )

root_generate_dictionary(G__Base ${HEADER_FILES} LINKDEF BaseLinkDef.h)

################################################################################################
# Create a shared library with geneated dictionary
################################################################################################
add_compile_options(-Wall -Wextra -Wpedantic)
add_library(Base SHARED ${HEADER_FILES} BaseLinkDef.h ${SOURCE_FILES} G__Base.cxx)

# optimization for big histogram access within the pair inner loop
# WARNING: the sumw2 structure is not used -> estimate uncertainties using sub-samples
add_definitions(-DOPTIMIZEADDBINCONTENT)

target_link_libraries(Base ${WAC_ROOT_LIBRARIES} ${EXTRA_LIBS} )
target_include_directories(Base  PUBLIC Base ${EXTRA_INCLUDES} ) 

install(FILES  "${CMAKE_CURRENT_BINARY_DIR}/libBase.rootmap" "${CMAKE_CURRENT_BINARY_DIR}/libBase_rdict.pcm" DESTINATION "$ENV{WAC_LIB}")
install(TARGETS Base  LIBRARY DESTINATION "$ENV{WAC_LIB}")


