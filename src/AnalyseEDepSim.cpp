#include <iostream>
#include <string>

#include "RDataFrameUtils.h"
// #include "GenieUtils.h"
// #include "GeoUtils.h"
#include "TFile.h"

TGeoManager* geo = nullptr;

//___________________________________________________________________
int main(int argc, char* argv[]){
    // Analyze output of genie production

    // if(argc < 4 || argc > 7)
    // {
    //     LOG("W", "AnalyseEDepSim -i <EDEP PRODUCTION> -g <GEOMETRY> -o <FILE OUTPUT>\n");
    //     throw "";
    // }

    // read user inputs

    // const char* fInput = argv[1];

    // const char* geometry = argv[2];

    // const char* fOutput = argv[3];

    unsigned int start = 0;

    // unsigned int stop = 1;

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
    //     }        else{
    //         auto ui = argv[++index];
    //         LOG("W", TString::Format("Unknown Input : %s", ui).Data());
    //         return 1; 
    //     }
    //     index++;
    // }

    // if you have multiple files enable multiple thread pocessing
    
    LOG("I", "Reading geometry");
    geo = TGeoManager::Import("/storage/gpfs_data/neutrino/users/gi/dunendggd/SAND_opt3_DRIFT1.root");

    auto fInput_genie = "/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc/events-in-SANDtracker.*.gtrac.root";
    auto fInput_edep = "/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc/events-in-SANDtracker.*.edep-sim.root";

    auto fOutput = "/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/events-in-SANDtracker.0.edep-sim.analysed.root";

    // if you have multiple files enable multiple thread pocessing
    if(TString::Format("%s",fInput_genie).Contains("*")){
        LOG("I","Enabling multiple threading");
        ROOT::EnableImplicitMT();
        geo->SetMaxThreads(100);
    };

    LOG("I", "Initialize ROOT DataFrame");
    auto chain_genie = RDFUtils::InitTChain(fInput_genie, "gRooTracker", start, start + 999u); 
    auto chain_edep = RDFUtils::InitTChain(fInput_edep, "EDepSimEvents", start, start + 999u); 
    chain_edep->AddFriend(chain_genie, "genie");
    auto df = RDFUtils::InitDF(chain_edep);

    // RDFUtils::PrintColumns(df);

    // throw "";

    auto dfC = RDFUtils::AddConstantsToDF(df); // add some columns with usefull constants

    auto dfGENIE = RDFUtils::GENIE::AddColumnsFromGENIE(dfC);

    auto dfEDEP = RDFUtils::EDEPSIM::NOSPILL::AddColumnsFromEDEPSIM(dfGENIE);
    
    LOG("I", "Writing ouput file");
    dfEDEP.Snapshot("edep_extended", fOutput, { 
                                                "FileName",
                                                "EventId",
                                                "EventType",
                                                "CCQEonHydrogen",
                                                "NofEvents",
                                                /*
                                                    GENIE INFO
                                                */
                                                "Interaction_vtxX",
                                                "Interaction_vtxY",
                                                "Interaction_vtxZ",
                                                "Interaction_vtxT",
                                                "IncomingNeutrinoP4",
                                                "InteractionVolume",
                                                "FinalStateLeptonPDG",
                                                "FinalStateLeptonNames",
                                                "FinalStateLepton4Momentum",
                                                "FinalStateLeptonEmissionAngle",
                                                "NofPrimaries",
                                                "FinalStateHadronicSystemPDG", 
                                                "FinalStateHadronicSystemNames", 
                                                "FinalStateHadronicSystemMomenta", 
                                                "FinalStateHadronicSystemTotal4Momentum", 
                                                "FinalStateHadronicSystemEmissionAngle", 
                                                "FinalStateHadronicSystemTotalKinE", 
                                                /*
                                                    EDEP INFO
                                                */
                                                "PrimariesPDG",
                                                "PrimariesTrackId",
                                                "PrimariesP4",
                                                "PrimariesFirstHitECAL",
                                                "PrimariesEDepECAL",
                                                "PrimariesEmissionAngle",
                                                /*
                                                    PREDICTIONS FOR CHANNEL antinu on H
                                                */
                                                "ExpectedNeutrinoP4FromMuon",
                                                "ExpectedHadronSystP3",
                                                "ExpectedHadronSystEnergy",
                                                "ExpectedNeutronArrivalPositionECAL",
                                                "MissingTransverseMomentum",
                                                "DoubleTransverseMomentumImbalance",

    });                                                    

    return 0;
}