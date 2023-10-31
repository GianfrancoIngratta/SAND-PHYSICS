#include <iostream>
#include <string>

#include "RDataFrameUtils.h"
#include "GeoUtils.h"
#include "TFile.h"

int main(int argc, char* argv[]){
    // Analyze output of sandreco 

    if(argc<4)
    {
        std::cout<<"AnalyseGenie <FAST RECO> <GEOMETRY> <FILE OUTPUT NAME>\n";
        throw "";
    }

    // read user inputs

    const char* fInput = argv[1];

    const char* geometry = argv[2];

    const char* fOutput = argv[3];

    // define output files name

    TString fOutput_ = TString::Format("%s.root",fOutput); 

    // if you have multiple files enable multiple thread pocessing

    if(!std::strstr(fInput,"*")) ROOT::EnableImplicitMT();

    // Initialize root DataFrame and add columns
    // auto df = RDFUtils::InitDF(fInput, "tEvent");
    auto df = RDFUtils::InitDF(fInput, "edepsmear_anaevt_tree");

    // RDFUtils::PrintColumns(df);

    auto dfC = RDFUtils::AddConstantsToDF(df); // add some columns with usefull constants

    auto dfSR = RDFUtils::FASTRECO::AddColumnsFromFASTRECO(dfC);

    dfSR.Snapshot("myTree", fOutput_.Data(), {"RecoMuonsVtxX",
                                              "RecoMuonsVtxY",
                                              "RecoMuonsVtxZ",
                                              "RecoMuonsP",
                                              "MuonsMomentumRes"
                                              }
    );
}