################################################################################################
#
# CMakeLists.txt for WAC package.
#
# Creates  libraries with dictionary and various executables
#
################################################################################################
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

################################################################################################
# Project WAC
################################################################################################
project(WAC VERSION 1.0 DESCRIPTION "Wayne State Code" LANGUAGES CXX)

set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(WAC_ROOT "${WAC_ROOT}")

#option(USE_ECCENTRICITY "Use eccentricity lib" OFF)

################################################################################################
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT CONFIG REQUIRED)
find_library(EG EG)
find_library(PYTHIA8LIB pythia8 PATHS $ENV{PYTHIA8}/lib)

add_compile_options(-pthread -std=c++20 -m64)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} $ENV{MYHOME} $ENV{ATGSI} "Base" "Best" "Pythia" "Urqmd" $ENV{PYTHIA8}/include)
set(WAC_ROOT_LIBRARIES ROOT::Core ROOT::MathCore ROOT::Hist ROOT::Tree ROOT::Graf ROOT::Physics ROOT::Gpad)

add_library(EG SHARED IMPORTED)
add_library(PYTHIA8LIB SHARED IMPORTED)
# address of Pythia8 xmldoc
add_definitions(-DPYTHIA8XMLPATH="$ENV{PYTHIA8}/share/Pythia8/xmldoc")


if (NOT DEFINED ENV{MYHOME} AND NOT DEFINED ENV{ATGSI})
find_library(EGPythia8 EGPythia8)
add_library(EGPythia8 SHARED IMPORTED)
set_target_properties(EGPythia8 PROPERTIES IMPORTED_LOCATION ${EGPythia8})
endif()

set_target_properties(EG PROPERTIES IMPORTED_LOCATION ${EG})
set_target_properties(PYTHIA8LIB PROPERTIES IMPORTED_LOCATION ${PYTHIA8LIB})


if (DEFINED ENV{MYHOME})
add_subdirectory(Home)
endif()
if (DEFINED ENV{ATGSI})
add_subdirectory(GSI)
endif()
add_subdirectory(Base)
add_subdirectory(Best)
if (NOT DEFINED ENV{ATWSU})
add_subdirectory(Pythia)
endif()
add_subdirectory(Urqmd)
add_subdirectory(Utils)

