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
        LOG("W", "AnalyseGenie -i <GENIE PRODUCTION> -g <GEOMETRY> -o <FILE OUTPUT>\n");
        throw "";
    }

    // target production
    // read user inputs

    const char* fInput = argv[1];

    const char* geometry = argv[2];

    const char* fOutput = argv[3];

    unsigned int start = 0;

    unsigned int stop = 1;

    int index = 1;

    LOG("I", "Parsing inputs");

    while (index < argc)
    {
        TString opt = argv[index];
        if(opt.CompareTo("-i")==0){
            try
            {
                fInput = argv[++index];
                LOG("ii",TString::Format("Input file : %s", fInput).Data());
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << '\n';
                return 1;
            }
        }else if(opt.CompareTo("-g")==0){
            try
            {
                geometry = argv[++index];
                LOG("ii",TString::Format("Geometry file : %s", geometry).Data());
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << '\n';
                return 1;
            }
        }else if(opt.CompareTo("-o")==0){
            try
            {
                fOutput = argv[++index];
                LOG("ii",TString::Format("Output file : %s", fOutput).Data());
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << '\n';
            }
        }        else{
            auto ui = argv[++index];
            LOG("W", TString::Format("Unknown Input : %s", ui).Data());
            return 1; 
        }
        index++;
    }

    LOG("I", "Reading geometry");        
    geo = TGeoManager::Import(geometry);

    // if you have multiple files enable multiple thread pocessing
    if(TString::Format("%s",fInput).Contains("*")){
        LOG("I","Enabling multiple threading");
        ROOT::EnableImplicitMT();
        geo->SetMaxThreads(100);
        stop = start + 10u;
        };

    LOG("I", "Initialize ROOT DataFrame");
    auto df = RDFUtils::InitDF(fInput, "gRooTracker", start, stop);

    auto dfC = RDFUtils::AddConstantsToDF(df); // add some columns with usefull constants

    auto dfG = RDFUtils::GENIE::AddColumnsFromGENIE(dfC);

    // auto dfG_SolidHydrogen = RDFUtils::GENIE::AddColumnsForHydrogenCarbonSampleSelection(dfG); // add columns for Hydrogen sample selection
    LOG("I", "Writing ouput file");
    dfG.Snapshot("gtrac_extended",fOutput, {
                                            "Interaction_vtxX",
                                            "Interaction_vtxY",
                                            "Interaction_vtxZ",
                                            "EventType",
                                            "NeutrinoFlavor",
                                            "InteractionTarget", // no thread safe
                                            "InteractionTargetFromGEO", // no thread safe
                                            "InteractionVolume", // no thread safe
                                            // incoming neutrinos
                                            "InitialStateNuMu_P4",
                                            "InitialStateAntiNuMu_P4",
                                            // outgoing muons
                                            "FinalStateMuonsP4",
                                            "FinalStateAntiMuonsP4",
                                            });

    // dfG.Snapshot("myTree",fOutput_.Data(), {"InteractionTarget",
    //                                "InteractionTargetFromGEO",
    //                                "EventType",
    //                                // initial state particles
    //                                "InitialStateNuMu_P4",
    //                                "InitialStateParticlesPDG",
    //                                "InitialStateParticlesE",
    //                                "InitialStateMomentum",
    //                                "InitialStateEnergy",
    //                                // stable final state particles
    //                                "StableFinalStateParticlesPDG",
    //                                "StableFinalStateParticlesE",
    //                                "StableFinalStateMomentum",
    //                                "StableFinalStateEnergy",
    //                                // final state nuclear remnant
    //                                "FinalStateNuclearRemnantPDG",
    //                                "FinalStateNuclearRemnantE",
    //                                "FinalStateNuclearMomentum",
    //                                "FinalStateNuclearEnergy",
    //                                });
    
    // dfG_SolidHydrogen.Snapshot("myTree", fOutput_4KinSelection.Data(), {
    //                                                                     "InteractionTarget",
    //                                                                     "BeamDirectionX",
    //                                                                     "BeamDirectionY",
    //                                                                     "BeamDirectionZ",
    //                                                                     "EventType",
    //                                                                     "FinalStateTopologyName",
    //                                                                     "FinalHadronicSystemP4_TT",
    //                                                                     "FinalStateMuonEmissionAngle",
    //                                                                     "FinalStateHadronicSystemEmissionAngle",
    //                                                                     "FinalStateHadronicSystemVSMuonsAngle",
    //                                                                     "InitialNucleonMomentum",
    //                                                                     "TransverseBoostingAngle",
    //                                                                     "Asimmetry_RmH",
    // });

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