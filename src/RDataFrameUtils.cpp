#include "RDataFrameUtils.h"
#include "GenieUtils.h"
#include "TChain.h"

//RDFUtils______________________________________________________________________

ROOT::RDataFrame RDFUtils::InitDF(const char* production, const char* tree_name){
    
    TChain* chain = new TChain(tree_name, tree_name);

    chain->Add(production);

    ROOT::RDataFrame* df = new ROOT::RDataFrame(*chain);

    if(df) std::cout<<"DF initialized\n";

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

ROOT::RDF::RNode RDFUtils::GENIE::AddColumnsFromGENIE(ROOT::RDataFrame& df){
    return df.Define("InteractionTarget", RDFUtils::GENIE::InteractionTarget, {"StdHepPdg"});
}