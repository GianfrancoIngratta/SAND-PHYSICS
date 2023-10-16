#include "RDataFrameUtils.h"
#include "GenieUtils.h"
#include "GeoUtils.h"
#include "TChain.h"

TGeoManager* geo = nullptr;

//RDFUtils______________________________________________________________________

ROOT::RDataFrame RDFUtils::InitDF(const char* production, const char* tree_name){
    
    TChain* chain = new TChain(tree_name, tree_name);

    chain->Add(production);

    ROOT::RDataFrame* df = new ROOT::RDataFrame(*chain);

    if(!df){
        std::cout<<"DF NOT initialized in file "<<__FILE__<<" line "<<__LINE__<<"\n";
        throw "";
        }

    return *df;
}

void RDFUtils::PrintColumns(ROOT::RDataFrame& df){

    auto colNames = df.GetColumnNames();
    
    for (auto &&colName : colNames){
        std::cout << "colName : " <<colName<<"\n";
        // std::cout << "colName : " <<colName<<", colType : " <<df.GetColumnType(colName)<<"\n";
        }
}

template<int coord>
ROOT::VecOps::RVec<double> RDFUtils::GetComponent(const ROOT::VecOps::RVec<TLorentzVector>& vTL){

    ROOT::VecOps::RVec<double> v;

    for(auto& TL : vTL) v.push_back(TL[coord]);

    return v;
}

ROOT::VecOps::RVec<double> RDFUtils::VectorDifference(const ROOT::VecOps::RVec<double>& v1,
                                                      const ROOT::VecOps::RVec<double>& v2){
    ROOT::VecOps::RVec<double> diff;
    for(auto i=0u; i<v1.size(); i++) diff.push_back(v1[i]-v2[i]);
    return diff;                                                        
}

ROOT::VecOps::RVec<double> RDFUtils::VectorSubtractConst(const ROOT::VecOps::RVec<double>& v1, double c){
    ROOT::VecOps::RVec<double> diff;
    for(auto i=0u; i<v1.size(); i++) diff.push_back(v1[i]-c);
    return diff;                                                        
}

TLorentzVector RDFUtils::VectorFilterByHighest(const ROOT::VecOps::RVec<double>& filter,
                                               const ROOT::VecOps::RVec<TLorentzVector>& v){
    // return value of v at index position corresponding to index of filter's highest
    if(v.size()==0){
        return {0.,0.,0.,0.};
    }else{                                      
        auto pointer2max = std::max_element(filter.begin(),filter.end());
        auto index2max = std::find(filter.begin(),filter.end(),*pointer2max);
        return v[index2max-filter.begin()];
    }
}

//RDFUtils::GENIE_________________________________________________________________

std::string RDFUtils::GENIE::InteractionTarget(const ROOT::VecOps::RVec<int>& pdg){
    return GenieUtils::PDG2Name(pdg[1]);
}

ROOT::VecOps::RVec<int> RDFUtils::GENIE::StableFinalStateParticles(const ROOT::VecOps::RVec<int>& pdg,
                                                                   const ROOT::VecOps::RVec<int>& status){
// get final state particles (to be tracked by detector-level)
    ROOT::VecOps::RVec<int> fsp;
    for (auto i = 0u; i < status.size(); i++)
    {
        if(status[i]==1) fsp.push_back(pdg[i]);
        // if(status[i]==1) fsp.push_back(GenieUtils::PDG2Name(pdg[i]));
    }

    return fsp;
}

int RDFUtils::GENIE::NofFinalStateParticles(const ROOT::VecOps::RVec<int>& pdg){
    return pdg.size();
}

std::string RDFUtils::GENIE::EventType(TObjString& s){
    
    if(s.GetString().Contains("DIS")){
        return "DIS";
    }else if(s.GetString().Contains("RES")){
        return "RES";
    }else if(s.GetString().Contains("QES")){
        return "QES";
    }else if(s.GetString().Contains("COH")){
        return "COH";
    }else{
        return "Other";
    }
}

ROOT::VecOps::RVec<genie::GHepParticle> RDFUtils::GENIE::AllGenieParticles(const ROOT::VecOps::RVec<int>& pdg,
                                                                           const ROOT::VecOps::RVec<int>& status,
                                                                           const ROOT::VecOps::RVec<double>& P4){
    ROOT::VecOps::RVec<genie::GHepParticle> particles;

    for(auto i=0u; i<pdg.size(); i++){
        
        TLorentzVector momentum = {P4[4*i],P4[4*i+1],P4[4*i+2],P4[4*i+3]};

        auto flag = static_cast<genie::GHepStatus_t>(status[i]);
        
        auto particle = GenieUtils::GenieParticle(pdg[i], flag, momentum);
        
        particles.push_back(particle);
    }

    return particles;                                                                            
}

ROOT::VecOps::RVec<genie::GHepParticle> RDFUtils::GENIE::StableFinalStateParticles(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    
    ROOT::VecOps::RVec<genie::GHepParticle> stables;
    
    for(auto& p : particles){
        if(p.Status()==genie::kIStStableFinalState) stables.push_back(p);
    }

    return stables;
}

template<int PDG>
ROOT::VecOps::RVec<TLorentzVector> RDFUtils::GENIE::GetMomentum(const ROOT::VecOps::RVec<int>& pdg,
                                                                const ROOT::VecOps::RVec<int>& status,
                                                                const ROOT::VecOps::RVec<double>& P4){
    ROOT::VecOps::RVec<TLorentzVector> v;

    for (auto i = 0u; i < pdg.size(); i++)
    {
        if(pdg[i]==PDG && status[i]==1){
            TLorentzVector p = {P4[4*i],P4[4*i+1],P4[4*i+2],P4[4*i+3]};
            v.push_back(p);
        }
    }

    return v;
}

TLorentzVector RDFUtils::GENIE::GetMomomentumHadronSystem(const ROOT::VecOps::RVec<int>& pdg,
                                                          const ROOT::VecOps::RVec<int>& status,
                                                          const ROOT::VecOps::RVec<double>& P4){
    TLorentzVector hadron_vector; 
    for (auto i = 0u; i < pdg.size(); i++)
    {
        // exclude muons, electrons, positrons, neutrinos, photons
        bool isNotLepton = (abs(pdg[i])!=13) && (abs(pdg[i])==14) && 
                           (abs(pdg[i])==11) && (abs(pdg[i])==22);
        bool isStableFinalState = (status[i]==1);
        if(isNotLepton && isStableFinalState){
            TLorentzVector running_hadron = {P4[4*i],P4[4*i+1],P4[4*i+2],P4[4*i+3]};
            hadron_vector += running_hadron;
        }
    }

    return hadron_vector;                                                               
}

double RDFUtils::GENIE::PImbalanceTrasverse2beam(TLorentzVector& p1, TLorentzVector& p2){
    /*
    take 2 particles and subtract their momentum transverse to the beam
    and return the magnitude of the difference
    */
   auto diff = p1 - p2;

   return sqrt(diff.X()*diff.X() + diff.Y()*diff.Y()); 
}

ROOT::RDF::RNode RDFUtils::GENIE::AddColumnsFromGENIE(ROOT::RDataFrame& df){
    return df.Define("Interaction_vtxX","EvtVtx[0]")
             .Define("Interaction_vtxY","EvtVtx[1]")
             .Define("Interaction_vtxZ","EvtVtx[2]")
             .Define("Interaction_vtxT","EvtVtx[3]")
             .Define("InteractionTarget", RDFUtils::GENIE::InteractionTarget, {"StdHepPdg"})
             .Define("InteractionTargetFromGEO", RDFUtils::GEO::GetMaterialFromCoordinates, {"Interaction_vtxX","Interaction_vtxY","Interaction_vtxZ"})
             //
             .Define("Particles", RDFUtils::GENIE::AllGenieParticles, {"StdHepPdg","StdHepStatus","StdHepP4"})
             .Define("StableFinalStateParticles", RDFUtils::GENIE::StableFinalStateParticles, {"Particles"})
             //
             .Define("EventType", RDFUtils::GENIE::EventType, {"EvtCode"})
             // muons
             .Define("MuonMomentum", RDFUtils::GENIE::GetMomentum<13>, {"StdHepPdg","StdHepStatus","StdHepP4"})
             .Define("MuonMomentumPX", RDFUtils::GetComponent<0>, {"MuonMomentum"})
             .Define("MuonMomentumPY", RDFUtils::GetComponent<1>, {"MuonMomentum"})
             .Define("MuonMomentumPZ", RDFUtils::GetComponent<2>, {"MuonMomentum"})
             .Define("MuonMomentumP",  RDFUtils::GetComponent<3>, {"MuonMomentum"})
             .Define("MuonHighestP", RDFUtils::VectorFilterByHighest, {"MuonMomentumP","MuonMomentum"}) // useless here, just 1 mu per event but ok
             // protons
             .Define("ProtonMomentum", RDFUtils::GENIE::GetMomentum<2122>, {"StdHepPdg","StdHepStatus","StdHepP4"})
             .Define("ProtonMomentumPX", RDFUtils::GetComponent<0>, {"ProtonMomentum"})
             .Define("ProtonMomentumPY", RDFUtils::GetComponent<1>, {"ProtonMomentum"})
             .Define("ProtonMomentumPZ", RDFUtils::GetComponent<2>, {"ProtonMomentum"})
             .Define("ProtonMomentumP", RDFUtils::GetComponent<3>, {"ProtonMomentum"})
             .Define("ProtonHighestP", RDFUtils::VectorFilterByHighest, {"ProtonMomentumP","ProtonMomentum"})
             // final hadron system
             .Define("HadronSystemMomentum", RDFUtils::GENIE::GetMomomentumHadronSystem, {"StdHepPdg","StdHepStatus","StdHepP4"})
             .Define("HadronSystemMomentumPX", "HadronSystemMomentum[0]")
             .Define("HadronSystemMomentumPY", "HadronSystemMomentum[1]")
             .Define("HadronSystemMomentumPZ", "HadronSystemMomentum[2]")
             .Define("HadronSystemMomentumP", "HadronSystemMomentum[3]")
             .Define("MuonHadronSystKinImbalance", RDFUtils::GENIE::PImbalanceTrasverse2beam, {"MuonHighestP","HadronSystemMomentum"})
             ;
}

//RDFUtils::GEO___________________________________________________________________

std::string RDFUtils::GEO::GetMaterialFromCoordinates(double x, double y, double z){
    
    if(!geo){
        std::cout<<"GEO not initialized in "<<__FILE__<<" "<<__LINE__<<"\n";
        throw "";
    }
   
    return geo->FindNode(x*100., y*100., z*100.)->GetVolume()->GetMaterial()->GetName();
}