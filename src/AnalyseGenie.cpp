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

    const char* fInput = argv[1];

    const char* geometry = argv[2];

    const char* fOutput = argv[3];

    if(!std::strstr(fInput,"*")) ROOT::EnableImplicitMT();

    auto df = RDFUtils::InitDF(fInput, "gRooTracker");

    geo = TGeoManager::Import(geometry);

    // RDFUtils::PrintColumns(df);

    auto dfC = RDFUtils::AddConstantsToDF(df);

    auto dfG = RDFUtils::GENIE::AddColumnsFromGENIE(dfC);

    dfG.Snapshot("myTree",fOutput,{"InteractionTarget",
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
                                   // topology & kinematic sutudies
                                   "FinalStateTopologyName",
                                   "FinalHadronicSystemP4_TT",
                                   "InitialNucleonMomentum",
                                   "TransverseBoostingAngle",
                                   });

    TString fOutput_1mu_1pr_1pi = TString::Format("test_1mu_1pr_1pi.root");

    // // example of topology name : "1mu_0pr_1ne_2pi_0em_0ex_0nu"
    dfG.Filter([](TString s){return s.Contains("1mu_1pr_0ne_1pi_0em_0ex_0nu");}, {"FinalStateTopologyName"})
       .Snapshot("selection", fOutput_1mu_1pr_1pi.Data(), {"FinalStateTopologyName",
                                                             "InteractionTarget",
                                                             "FinalHadronicSystemP4_TT",
                                                            });                                   

    return 0;
}
//___________________________________________________________________