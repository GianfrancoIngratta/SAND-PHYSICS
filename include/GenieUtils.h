#ifndef GENIE_UTILS_H
#define GENIE_UTILS_H

#include <string>

#include "TString.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"

namespace GenieUtils{

std::string PDG2Name(int pdg);

genie::GHepParticle GenieParticle(int pdg, 
                                  genie::GHepStatus_t status, 
                                  const TLorentzVector& momentum);

// class GenieEvent
// {
//     private:
        
//         int gNofMuons; 
//         int gNofProton; 
//         int gNofNeutrons; 
//         int gNofChargedPions; 
//         int gNofPi0;
//         int gNofGamma;
//         int gNofElectronsPositrons;
//         int gNofHyperons;
//         int gNofHeavyMesons;
//         int gNofNuclei;
        
//         const char* gEventName;

//         genie::GHepParticle gMuon;
//         std::vector<genie::GHepParticle> gHadronSystem;

//     public:
// }                                  

}

#endif