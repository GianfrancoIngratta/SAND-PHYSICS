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

double RDFUtils::GetColumnSum(const ROOT::VecOps::RVec<double>& v){
    double sum = 0.;
    for(auto& value : v) sum+=value;
    return sum;
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

template<genie::GHepStatus_t STATUS>
ROOT::VecOps::RVec<genie::GHepParticle> RDFUtils::GENIE::GetParticlesWithStatus(const ROOT::VecOps::RVec<genie::GHepParticle>& particles){
    
    ROOT::VecOps::RVec<genie::GHepParticle> selected;

    for(auto& p : particles){
        if(p.Status() == STATUS) selected.push_back(p);
    }

    return selected;
}

ROOT::VecOps::RVec<int> RDFUtils::GENIE::GetPDG(const ROOT::VecOps::RVec<genie::GHepParticle>& particles)
{
    ROOT::VecOps::RVec<int> pdgs;
    for(auto& p : particles) pdgs.push_back(p.Pdg());
    return pdgs;
}

template<int component>
ROOT::VecOps::RVec<double> RDFUtils::GENIE::GetMomentum(const ROOT::VecOps::RVec<genie::GHepParticle>& particles)
{
    // input : particles
    // output : momentum component of each particle
    ROOT::VecOps::RVec<double> momenta;
    for(auto& p : particles){
        TLorentzVector momentum = *p.GetP4();
        momenta.push_back(momentum[component]);
        }
    return momenta;
}

ROOT::RDF::RNode RDFUtils::GENIE::AddColumnsFromGENIE(ROOT::RDataFrame& df){
    return df.Define("Interaction_vtxX","EvtVtx[0]")
             .Define("Interaction_vtxY","EvtVtx[1]")
             .Define("Interaction_vtxZ","EvtVtx[2]")
             .Define("Interaction_vtxT","EvtVtx[3]")
             .Define("EventType", RDFUtils::GENIE::EventType, {"EvtCode"})
             .Define("InteractionTarget", RDFUtils::GENIE::InteractionTarget, {"StdHepPdg"})
             .Define("InteractionTargetFromGEO", RDFUtils::GEO::GetMaterialFromCoordinates, {"Interaction_vtxX","Interaction_vtxY","Interaction_vtxZ"})
             // particles
             .Define("Particles", RDFUtils::GENIE::AllGenieParticles, {"StdHepPdg","StdHepStatus","StdHepP4"})
             //initial state
             .Define("InitialStateParticles", RDFUtils::GENIE::GetParticlesWithStatus<genie::kIStInitialState>, {"Particles"})
             .Define("InitialStateParticlesPDG", RDFUtils::GENIE::GetPDG, {"InitialStateParticles"})
             .Define("InitialStateParticlesE", RDFUtils::GENIE::GetMomentum<3>, {"InitialStateParticles"})
             .Define("InitialStateEnergy", RDFUtils::GetColumnSum, {"InitialStateParticlesE"})
             //final state state particles
             .Define("StableFinalStateParticles", RDFUtils::GENIE::GetParticlesWithStatus<genie::kIStStableFinalState>, {"Particles"})
             .Define("StableFinalStateParticlesPDG", RDFUtils::GENIE::GetPDG, {"StableFinalStateParticles"})
             .Define("StableFinalStateParticlesE", RDFUtils::GENIE::GetMomentum<3>, {"StableFinalStateParticles"})
             .Define("StableFinalStateEnergy", RDFUtils::GetColumnSum, {"StableFinalStateParticlesE"})
             // final state nuclear remnant
             .Define("FinalStateNuclearRemnant", RDFUtils::GENIE::GetParticlesWithStatus<genie::kIStFinalStateNuclearRemnant>, {"Particles"})
             .Define("FinalStateNuclearRemnantPDG", RDFUtils::GENIE::GetPDG, {"FinalStateNuclearRemnant"})
             .Define("FinalStateNuclearRemnantE", RDFUtils::GENIE::GetMomentum<3>, {"FinalStateNuclearRemnant"})
             .Define("FinalStateNuclearRemnantEnergy", RDFUtils::GetColumnSum, {"FinalStateNuclearRemnantE"})
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