#include "GenieUtils.h"

std::string GenieUtils::PDG2Name(int pdg){
    genie::PDGLibrary* instance = genie::PDGLibrary::Instance();
    std::string name = instance->Find(pdg)->GetName();
    return name;
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