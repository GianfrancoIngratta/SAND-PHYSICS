#include <iostream>
#include <string>

#include "RDataFrameUtils.h"
#include "GenieUtils.h"
#include "GeoUtils.h"
#include "TFile.h"

TGeoManager* geo = nullptr;

//___________________________________________________________________
int main(int argc, char* argv[]){
    // Analyze output of genie production

    if(argc<4)
    {
        std::cout<<"AnalyseGenie <GENIE PRODUCTION> <GEOMETRY> <FILE OUTPUT NAME>\n";
        throw "";
    }

    // read user inputs

    const char* fInput = argv[1];

    const char* geometry = argv[2];

    const char* fOutput = argv[3];

    // if you have multiple files enable multiple thread pocessing

    if(!std::strstr(fInput,"*")) ROOT::EnableImplicitMT();
    
    geo = TGeoManager::Import(geometry);

    // define output files name

    TString fOutput_ = TString::Format("%s.root",fOutput); 

    TString fOutput_4KinSelection = TString::Format("%s_4KinSelection.root",fOutput); 
    
    // Initialize root DataFrame and add columns

    auto df = RDFUtils::InitDF(fInput, "gRooTracker");

    // RDFUtils::PrintColumns(df);

    auto dfC = RDFUtils::AddConstantsToDF(df); // add some columns with usefull constants

    auto dfG = RDFUtils::GENIE::AddColumnsFromGENIE(dfC);

    auto dfG_SolidHydrogen = RDFUtils::GENIE::AddColumnsForHydrogenCarbonSampleSelection(dfG); // add columns for Hydrogen sample selection

    dfG.Snapshot("myTree",fOutput_.Data(), {"InteractionTarget",
                                   "InteractionTargetFromGEO",
                                   "EventType",
                                   // initial state particles
                                   "InitialStateNeutrinoP4",
                                   "InitialStateParticlesPDG",
                                   "InitialStateParticlesE",
                                   "InitialStateMomentum",
                                   "InitialStateEnergy",
                                   // stable final state particles
                                   "StableFinalStateParticlesPDG",
                                   "StableFinalStateParticlesE",
                                   "StableFinalStateMomentum",
                                   "StableFinalStateEnergy",
                                   // final state nuclear remnant
                                   "FinalStateNuclearRemnantPDG",
                                   "FinalStateNuclearRemnantE",
                                   "FinalStateNuclearMomentum",
                                   "FinalStateNuclearEnergy",
                                   });
    
    dfG_SolidHydrogen.Snapshot("myTree", fOutput_4KinSelection.Data(), {"InteractionTarget",
                                                                        "EventType",
                                                                        "FinalStateTopologyName",
                                                                        "FinalHadronicSystemP4_TT",
                                                                        "InitialNucleonMomentum",
                                                                        "TransverseBoostingAngle",
                                                                        "Asimmetry_RmH",
    });

    // TString fOutput_1mu_1pr_1pi = TString::Format("test_1mu_1pr_1pi.root");

    // // // example of topology name : "1mu_0pr_1ne_2pi_0em_0ex_0nu"
    // dfG.Filter([](TString s){return s.Contains("1mu_1pr_0ne_1pi_0em_0ex_0nu");}, {"FinalStateTopologyName"})
    //    .Snapshot("selection", fOutput_1mu_1pr_1pi.Data(), {"FinalStateTopologyName",
    //                                                          "InteractionTarget",
    //                                                          "FinalHadronicSystemP4_TT",
    //                                                         });                                   

    return 0;
}
//___________________________________________________________________