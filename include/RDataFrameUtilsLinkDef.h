#ifdef __CINT__

#include "RDataFrameUtils.h"
#include "GenieUtils.h"
#include "EDepSimUtils.h"

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ namespace genie;

#pragma link C++ class ROOT::VecOps::RVec < genie::GHepParticle> + ;
#pragma link C++ class ROOT::VecOps::RVec < std::string> + ;
#pragma link C++ class ROOT::VecOps::RVec < TLorentzVector> + ;
#pragma link C++ class ROOT::VecOps::RVec < event_topology> + ;
#pragma link C++ class ROOT::VecOps::RVec < track_hits> + ;
#pragma link C++ class event_topology + ;
#pragma link C++ class track_hits + ;

#endif