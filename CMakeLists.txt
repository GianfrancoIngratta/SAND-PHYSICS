cmake_minimum_required(VERSION 3.0.0)
project(SANDPHYSICS VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_PLATFORM_INDEPENDENT_CODE ON)

# set output directory
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_INSTALL_LIBDIR "${CMAKE_INSTALL_PREFIX}/lib")

# find_package(ROOT REQUIRED COMPONENTS RIO Net)
find_package(ROOT REQUIRED COMPONENTS Geom Physics Matrix MathCore RIO Net Tree)
find_package(EDepSim)

# find_package(SandReco)
include_directories(${PROJECT_SOURCE_DIR}/include /storage/gpfs_data/neutrino/users/gi/sand-reco/include)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include" "${CMAKE_CURRENT_SOURCE_DIR}")

# Creates a libGenieUtils shared library
add_library(GenieUtils SHARED src/GenieUtils.cpp)
target_include_directories(GenieUtils PUBLIC /opt/exp_software/neutrino/GENIEv2/Generator/src)
target_link_directories(GenieUtils PUBLIC /opt/exp_software/neutrino/GENIEv2/Generator/lib /opt/exp_software/neutrino/LHAPDF5/lib) 
target_link_libraries(GenieUtils PUBLIC ROOT::RIO ROOT::Net ROOT::Tree ROOT::TreePlayer ROOT::ROOTDataFrame
                      ROOT::Physics ROOT::Geom ROOT::EG ROOT::EGPythia6
                      GMessenger GRegistry GAlgorithm GInteraction GHEP GBase GNumerical GUtils
                      GPDG GBaryonResonance GEVGCore GEVGDrivers GNtuple GGeo GFluxDrivers GPDF
                      GElFF GDecay GFragmentation GNuclear GLlewellynSmith GCrossSections GCharm
                      GElas GGiBUU GReinSehgal GQPM GBodekYang GEVGModules GQEL GRES GDIS GCoh
                      GDfrc GMEC GNuE GNuGamma GAlvarezRuso GHadronTransp GBergerSehgal GSingleKaon GMuELoss  
                      xml2
                      log4cpp
                      EDepSim::edepsim_io)

add_library(GeoUtils SHARED src/GeoUtils.cpp)
target_link_libraries(GeoUtils PUBLIC ROOT::RIO ROOT::Net ROOT::ROOTDataFrame EDepSim::edepsim_io)

# add_library(SelectSignalUtils SHARED src/SelectSignal.cpp)
# target_link_libraries(SelectSignalUtils PUBLIC ROOT::RIO ROOT::Net ROOT::ROOTDataFrame EDepSim::edepsim_io GenieUtils)

# Creates a libUtils shared library
add_library(RDataFrameUtils SHARED src/RDataFrameUtils.cpp src/GenieUtils.cpp src/GeoUtils.cpp)
# add_library(RDataFrameUtils SHARED src/RDataFrameUtils.cpp)
target_include_directories(RDataFrameUtils PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
  "$<INSTALL_INTERFACE:include>"
  /opt/exp_software/neutrino/GENIEv2/Generator/src
  /opt/exp_software/neutrino/LHAPDF5/lib)
target_link_libraries(RDataFrameUtils PUBLIC ROOT::RIO ROOT::Net ROOT::Tree ROOT::TreePlayer ROOT::ROOTDataFrame
ROOT::Physics ROOT::Geom ROOT::EG ROOT::EGPythia6 GenieUtils GeoUtils)
ROOT_GENERATE_DICTIONARY(RDataFrameUtilsDict RDataFrameUtils.h GHEP/GHepParticle.h 
MODULE RDataFrameUtils LINKDEF include/RDataFrameUtilsLinkDef.h)

add_executable(AnalyseGenie src/AnalyseGenie.cpp)
target_link_libraries(AnalyseGenie PUBLIC ROOT::RIO ROOT::Net ROOT::Tree ROOT::TreePlayer ROOT::ROOTDataFrame ROOT::Physics RDataFrameUtils)

add_executable(AnalyseEDepSim src/AnalyseEDepSim.cpp)
target_link_libraries(AnalyseEDepSim PUBLIC ROOT::RIO ROOT::Net ROOT::Tree ROOT::TreePlayer ROOT::ROOTDataFrame ROOT::Physics GenieUtils RDataFrameUtils GeoUtils)

add_executable(AnalyseDigit src/AnalyseDigitize.cpp)
target_link_libraries(AnalyseDigit PUBLIC ROOT::RIO ROOT::Net ROOT::Tree ROOT::TreePlayer ROOT::ROOTDataFrame ROOT::Physics GenieUtils RDataFrameUtils GeoUtils)

add_executable(AnalyseDriftReco src/AnalyseDriftReco.cpp)
target_link_libraries(AnalyseDriftReco PUBLIC ROOT::RIO ROOT::Net ROOT::Tree ROOT::TreePlayer ROOT::ROOTDataFrame ROOT::Physics GenieUtils RDataFrameUtils GeoUtils)

add_executable(SelectSignal src/SelectSignal.cpp)
target_link_libraries(SelectSignal PUBLIC ROOT::RIO ROOT::Net ROOT::Tree ROOT::TreePlayer ROOT::ROOTDataFrame ROOT::Physics GenieUtils RDataFrameUtils GeoUtils)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "default install path" FORCE )
endif()

install(TARGETS  GenieUtils RDataFrameUtils GeoUtils AnalyseGenie
        EXPORT SandPhysicsTargets
        RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")

install(
  DIRECTORY "${CMAKE_BINARY_DIR}/."
  DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
  FILES_MATCHING
  PATTERN "*.pcm"
  PATTERN "*.rootmap"
  PATTERN "lib" EXCLUDE
  PATTERN "bin" EXCLUDE
  PATTERN "CMakeFiles" EXCLUDE)

install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include
  DESTINATION "${CMAKE_INSTALL_PREFIX}"
  PATTERN "RDataFrameUtilsLinkDef.h" EXCLUDE)

  # ############################################################
#
# Prepare the package so that it can be used with the find_package interface.
#
# ############################################################

# Include module with function 'write_basic_package_version_file'
include(CMakePackageConfigHelpers)

# Build the targets description so that the package can be configured
# using find_package.
install(EXPORT SandPhysicsTargets
  NAMESPACE SandPhysics::
  DESTINATION lib/cmake/SandPhysics)

# Write the 'SandPhysicsConfigVersion.cmake' file which can be used to
# check if a version meets the requested properties.
write_basic_package_version_file(
  SandPhysicsConfigVersion.cmake
  COMPATIBILITY SameMajorVersion)

# Write the 'SandPhysicsConfig.cmake' file so that a user package can
# access this with find_package.
# configure_package_config_file(
#   PackageConfig.cmake.in
#   SandPhysicsConfig.cmake
#   PATH_VARS CMAKE_INSTALL_PREFIX
#   INSTALL_DESTINATION lib/cmake/SandPhysics)

# Install the config files.
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/SandPhysicsConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/SandPhysicsConfigVersion.cmake
  DESTINATION lib/cmake/SandPhysics)
