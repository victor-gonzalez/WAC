# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT montecarlo/pythia8 package
# @author Pere Mato, CERN
############################################################################

set(HEADER_FILES
  TPythia8.h
  TPythia8Decayer.h
  )

set(SOURCE_FILES
  TPythia8.cxx
  TPythia8Decayer.cxx
  )

root_generate_dictionary(G__EGPythia8 ${HEADER_FILES} LINKDEF LinkDef.h)

add_compile_options(-Wno-delete-non-virtual-dtor -Wno-misleading-indentation -Wno-overloaded-virtual)
add_library(EGPythia8 SHARED ${HEADER_FILES} LinkDef.h ${SOURCE_FILES} G__EGPythia8.cxx)


target_link_libraries(EGPythia8 PUBLIC PYTHIA8LIB)
target_include_directories(EGPythia8 PUBLIC Home ${EXTRA_INCLUDES})
install(FILES  "${CMAKE_CURRENT_BINARY_DIR}/libEGPythia8.rootmap" "${CMAKE_CURRENT_BINARY_DIR}/libEGPythia8_rdict.pcm" DESTINATION "$ENV{WAC_LIB}")
install(TARGETS EGPythia8  LIBRARY DESTINATION "$ENV{WAC_LIB}")

