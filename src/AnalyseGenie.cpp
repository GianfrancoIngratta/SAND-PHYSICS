#include <iostream>
#include <string>

#include "RDataFrameUtils.h"
#include "TFile.h"

TGeoManager* geo = nullptr;

//___________________________________________________________________
int main(int argc, char* argv[]){
    // Analyze output of genie production

    // if(argc<4)
    // {
    //     LOG("W", "AnalyseGenie -i <GENIE PRODUCTION> -g <GEOMETRY> -o <FILE OUTPUT>\n");
    //     throw "";
    // }

    // // read user inputs

    // const char* fInput = argv[1];

    // const char* geometry = argv[2];

    // const char* fOutput = argv[3];

    unsigned int start = 0;

    // int index = 1;

    // LOG("I", "Parsing inputs");

    // while (index < argc)
    // {
    //     TString opt = argv[index];
    //     if(opt.CompareTo("-i")==0){
    //         try
    //         {
    //             fInput = argv[++index];
    //             LOG("ii",TString::Format("Input file : %s", fInput).Data());
    //         }
    //         catch(const std::exception& e)
    //         {
    //             std::cerr << e.what() << '\n';
    //             return 1;
    //         }
    //     }else if(opt.CompareTo("-g")==0){
    //         try
    //         {
    //             geometry = argv[++index];
    //             LOG("ii",TString::Format("Geometry file : %s", geometry).Data());
    //         }
    //         catch(const std::exception& e)
    //         {
    //             std::cerr << e.what() << '\n';
    //             return 1;
    //         }
    //     }else if(opt.CompareTo("-o")==0){
    //         try
    //         {
    //             fOutput = argv[++index];
    //             LOG("ii",TString::Format("Output file : %s", fOutput).Data());
    //         }
    //         catch(const std::exception& e)
    //         {
    //             std::cerr << e.what() << '\n';
    //         }
    //     }else{
    //         auto ui = argv[++index];
    //         LOG("W", TString::Format("Unknown Input : %s", ui).Data());
    //         return 1; 
    //     }
    //     index++;
    // }

    LOG("I", "Reading geometry");        
    geo = TGeoManager::Import("/storage/gpfs_data/neutrino/users/gi/dunendggd/SAND_opt3_DRIFT1.root");

    auto fInput = "/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc/events-in-SANDtracker.*.gtrac.root";

    auto fOutput = "/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/events-in-SANDtracker.0.gtrac.analysed.root";

    // if you have multiple files enable multiple thread pocessing
    if(TString::Format("%s",fInput).Contains("*")){
        LOG("I","Enabling multiple threading");
        ROOT::EnableImplicitMT();
        geo->SetMaxThreads(100);
        };

    LOG("I", "Initialize ROOT DataFrame");
    auto chain_genie = RDFUtils::InitTChain(fInput, "gRooTracker", start, start + 999u); 
    auto df = RDFUtils::InitDF(chain_genie);

    auto dfC = RDFUtils::AddConstantsToDF(df); // add some columns with usefull constants

    auto dfG = RDFUtils::GENIE::AddColumnsFromGENIE(dfC);

    // auto dfG_SolidHydrogen = RDFUtils::GENIE::AddColumnsForHydrogenCarbonSampleSelection(dfG); // add columns for Hydrogen sample selection
    LOG("I", "Writing ouput file");
    dfG.Snapshot("gtrac_extended",fOutput, {
                                            "EvtNum",
                                            "NeutrinoFlavor",
                                            "isCCEvent",
                                            "EventType",
                                            "InteractionTargetPDG",
                                            "CCQEonHydrogen",
                                            "Interaction_vtxX",
                                            "Interaction_vtxY",
                                            "Interaction_vtxZ",
                                            "isInFiducialVolume",
                                            "InteractionTarget",
                                            "InteractionTargetFromGEO",
                                            "InteractionVolume",
                                            // initial state
                                            "InitialStatePDG",
                                            "InitialStateNames",
                                            "InitialStateMomenta",
                                            "InitialStateTotal4Momentum",
                                            "InitialStateEmissionAngle",
                                            // final state lepton
                                            "FinalStateLeptonPDG",
                                            "FinalStateLeptonNames",
                                            "NofFinalStateChargedParticles",
                                            "FinalStateLepton4Momentum",
                                            "FinalStateLeptonEmissionAngle",
                                            // prmary state hadronic system 
                                            "PrimaryStateHadronicSystemPDG",
                                            "PrimaryStateHadronicSystemNames",
                                            "PrimaryStateHadronicSystemMomenta",
                                            "PrimaryStateHadronicSystemTotal4Momentum",
                                            "PrimaryStateHadronicSystemEmissionAngle",
                                            "PrimaryStateHadronicSystemTotalKinE",
                                            "PrimaryStateHadronicSystemTopology_code",
                                            "PrimaryStateHadronicSystemTopology_name",
                                            // final state hadronic system 
                                            "FinalStateHadronicSystemPDG",
                                            "FinalStateHadronicSystemNames",
                                            "FinalStateHadronicSystemMomenta",
                                            "FinalStateHadronicSystemTotal4Momentum",
                                            "FinalStateHadronicSystemEmissionAngle",
                                            "FinalStateHadronicSystemTotalKinE",
                                            "FinalStateHadronicSystemTopology_code",
                                            "FinalStateHadronicSystemTopology_name",
                                            // nuclear remnant
                                            "NuclearRemnantPDG",
                                            "NuclearRemnantMomenta",
                                            "NuclearTotal4Momentum",
                                            // checks energy momentum conservation
                                            "PrimaryStateTotal4Momantum",
                                            "FinalStateTotal4Momantum",
                                            // for channel selection
                                            "ExpectedNeutrinoP4FromMuon",
                                            "ExpectedHadronSystP3",
                                            "ExpectedHadronSystEnergy",
                                            "FinalStateLeptonTransverseP",
                                            "FinalStateHadronicSystemTransverseP",
                                            "MissingTransverseMomentum",
                                            "RmH",
                                            "DoubleTransverseMomentumImbalance",
                                            "ExpectedNeutronArrivalPositionECAL",
                                            });

    return 0;
}
//___________________________________________________________________