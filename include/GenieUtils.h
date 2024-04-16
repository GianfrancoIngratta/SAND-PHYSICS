#ifndef GENIE_UTILS_H
#define GENIE_UTILS_H

#include <string>
#include <mutex>

#include "TString.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"

namespace GenieUtils{

auto database = genie::PDGLibrary::Instance();

std::string PDG2Name(int pdg);

genie::GHepParticle GenieParticle(int pdg, 
                                  genie::GHepStatus_t status, 
                                  const TLorentzVector& momentum);

struct event_topology
{
    int NofMuons = 0;

    int NofProtons = 0;

    int NofNeutrons = 0;
    // pi+- pi0
    int NofPions = 0;
    // e+- gamma
    int NofElPosGamma = 0;
    // mesons : K0, KL, K+- ..., hadrons : Sigma, lambda...
    int NofExhotic;
    // C12, 016, ...
    int NofRecoiledNuclei = 0;

    TString GetTopologyName(){
        return  TString::Format("%dmu_%dpr_%dne_%dpi_%dem_%dex_%dnu",
        NofMuons,NofProtons,NofNeutrons,NofPions,NofElPosGamma,NofExhotic,NofRecoiledNuclei);
    };
};                              

}

#endif