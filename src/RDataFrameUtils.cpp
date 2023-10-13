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
    
    for (auto &&colName : colNames) 
        std::cout << "colName : " <<colName<<"\n";
        // std::cout << "colName : " <<colName<<", colType : " <<df.GetColumnType(colName)<<std::endl;
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
    }

    return fsp;
}

int RDFUtils::GENIE::NofFinalStateParticles(const ROOT::VecOps::RVec<int>& pdg){
    return pdg.size();
}

ROOT::RDF::RNode RDFUtils::GENIE::AddColumnsFromGENIE(ROOT::RDataFrame& df){
    return df.Define("Interaction_vtxX","EvtVtx[0]")
             .Define("Interaction_vtxY","EvtVtx[1]")
             .Define("Interaction_vtxZ","EvtVtx[2]")
             .Define("Interaction_vtxT","EvtVtx[3]")
             .Define("InteractionTarget", RDFUtils::GENIE::InteractionTarget, {"StdHepPdg"})
             .Define("InteractionTargetFromGEO", RDFUtils::GEO::GetMaterialFromCoordinates, {"Interaction_vtxX",
                                                                                            "Interaction_vtxY",
                                                                                            "Interaction_vtxZ"})
             .Define("StableFinalStateParticles", RDFUtils::GENIE::StableFinalStateParticles,{"StdHepPdg",
                                                                                              "StdHepStatus"})
             .Define("NofFinalStateParticles", RDFUtils::GENIE::NofFinalStateParticles, {"StableFinalStateParticles"});                                                                                             
}

//RDFUtils::GEO___________________________________________________________________

std::string RDFUtils::GEO::GetMaterialFromCoordinates(double x, double y, double z){
    
    if(!geo){
        std::cout<<"GEO not initialized in "<<__FILE__<<__LINE__<<"\n";
        throw "";
    }

    return geo->FindNode(x*100., y*100., z*100.)->GetVolume()->GetMaterial()->GetName();
}