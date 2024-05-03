#include "GenieUtils.h"

std::string GenieUtils::PDG2Name(int pdg){
    return database->Find(pdg)->GetName();
}

void GenieUtils::UpdateTopology(GenieUtils::event_topology& t, int pdg){
    if(pdg == genie::kPdgProton || pdg == genie::kPdgAntiProton ){
        t.NofProtons++;
    }else if(pdg == genie::kPdgNeutron || pdg == genie::kPdgAntiNeutron){
        t.NofNeutrons++;
    }else if(pdg == genie::kPdgPiP){
        t.NofPiP++;
    }else if(pdg == genie::kPdgPiM){
        t.NofPiM++;
    }else if(pdg == genie::kPdgPi0){
        t.NofPi0++;
    }else if(pdg == genie::kPdgGamma || pdg == genie::kPdgElectron || pdg == genie::kPdgPositron){
        t.NofEM++;
    }else if(genie::pdg::IsKaon(pdg)){
        t.NofKaons++; 
    }else if(genie::pdg::IsBaryonResonance(pdg)){
        t.NofResonance++;
    }else{
        t.NofOther++;
    }
}