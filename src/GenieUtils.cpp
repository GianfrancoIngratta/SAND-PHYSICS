#include "GenieUtils.h"

std::string GenieUtils::PDG2Name(int pdg){
    return database->Find(pdg)->GetName();
}

Double_t GenieUtils::GetMass(int pdg){
    return database->Find(pdg)->Mass();;
}

Double_t GenieUtils::GetCharge(int pdg){
    return database->Find(pdg)->Charge();;
}

std::string GenieUtils::ShortEventType(const TString& s){

    if(s.Contains("Weak[CC]")){
        if(s.Contains(",RES;")){
            return "CC_RES";
        }else if(s.Contains(",DIS;")){
            return "CC_DIS";
        }else if(s.Contains(",QES;")){
            return "CC_QES";
        }else if(s.Contains(",COH;")){
            return "CC_COH";
        }else{
            return "CC_other";
        }
    }else if(s.Contains("Weak[NC]")){
        if(s.Contains(",RES;")){
            return "NC_RES";
        }else if(s.Contains(",DIS;")){
            return "NC_DIS";
        }else if(s.Contains(",QES;")){
            return "NC_QES";
        }else if(s.Contains(",COH;")){
            return "NC_COH";
        }else{
            return "NC_other";
        }
    }else{
        return "Unknown";
    }
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