#include "GenieUtils.h"

std::string GenieUtils::PDG2Name(int pdg){
    return database->Find(pdg)->GetName();
}

genie::GHepParticle GenieUtils::GenieParticle(int pdg, 
                                              genie::GHepStatus_t status, 
                                              const TLorentzVector& momentum){
    genie::GHepParticle p;
    p.SetPdgCode(pdg);
    p.SetStatus(status);
    p.SetMomentum(momentum);
    return p;
}