cmake_minimum_required(VERSION 3.0.0)
project(SANDPHYSICS VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_PLATFORM_INDEPENDENT_CODE ON)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# find_package(ROOT REQUIRED COMPONENTS RIO Net)
find_package(ROOT REQUIRED COMPONENTS Geom Physics Matrix MathCore RIO Net Tree)
find_package(EDepSim)

find_package(SandReco)
include_directories(${PROJECT_SOURCE_DIR}/include /storage/gpfs_data/neutrino/users/gi/sand-reco/include)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include" "${CMAKE_CURRENT_SOURCE_DIR}" "/storage/gpfs_data/neutrino/users/gi/STT-DESIGN-STUDY/FastReco")

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
                      log4cpp)

# Creates a libSandRecoUtils shared library
add_library(SandRecoUtils SHARED src/SandRecoUtils.cpp)
target_include_directories(SandRecoUtils PUBLIC /storage/gpfs_data/neutrino/users/gi/sand-reco/include)
target_link_directories(SandRecoUtils PUBLIC /storage/gpfs_data/neutrino/users/gi/sand-reco/lib)
target_link_libraries(SandRecoUtils PUBLIC ROOT::RIO ROOT::Net ROOT::EG)

# Creates a libFastRecoUtils shared library
add_library(evtinfo SHARED IMPORTED) # or STATIC instead of SHARED
set_target_properties(evtinfo PROPERTIES
  IMPORTED_LOCATION "/storage/gpfs_data/neutrino/users/gi/STT-DESIGN-STUDY/FastReco/build/libevtinfo.so"
  INTERFACE_INCLUDE_DIRECTORIES "/storage/gpfs_data/neutrino/users/gi/STT-DESIGN-STUDY/FastReco"
)

# Creates a libUtils shared library
add_library(RDataFrameUtils SHARED src/RDataFrameUtils.cpp)
target_link_libraries(RDataFrameUtils PUBLIC ROOT::RIO ROOT::Net ROOT::ROOTDataFrame GenieUtils SandRecoUtils evtinfo)

add_library(GeoUtils SHARED src/GeoUtils.cpp)
target_link_libraries(GeoUtils PUBLIC ROOT::RIO ROOT::Net ROOT::ROOTDataFrame EDepSim::edepsim_io)

add_executable(AnalyseGenie src/AnalyseGenie.cpp)
target_link_libraries(AnalyseGenie PUBLIC ROOT::RIO ROOT::Net ROOT::Tree ROOT::TreePlayer ROOT::ROOTDataFrame ROOT::Physics GenieUtils RDataFrameUtils GeoUtils)

add_executable(AnalyseSandReco src/AnalyseSandReco.cpp)
target_link_libraries(AnalyseSandReco PUBLIC ROOT::RIO ROOT::Net ROOT::Tree ROOT::TreePlayer ROOT::ROOTDataFrame ROOT::Physics 
GenieUtils RDataFrameUtils GeoUtils SandRecoUtils)

install(TARGETS  GenieUtils RDataFrameUtils AnalyseGenie AnalyseSandReco
        RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")