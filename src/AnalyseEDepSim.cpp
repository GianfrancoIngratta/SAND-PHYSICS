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

    // analyse file from file file_start" to "file_start + nof files 
    unsigned int index;
    if(argc != 2){
        LOG("W", "One input number needed to run the executable");
        throw "";
    }

    index = atoi(argv[1]);
    unsigned int files_per_jobs = 10u;
    unsigned int file_start = index * files_per_jobs;
    unsigned int file_stop = index * files_per_jobs + files_per_jobs;

    LOG("I", TString::Format("Analyze production from file %d to file %d", file_start, file_stop).Data());
    
    LOG("I", "Reading geometry");
    geo = TGeoManager::Import("/storage/gpfs_data/neutrino/users/gi/dunendggd/SAND_opt3_DRIFT1.root");

    auto FOLDER_PRODUCTION = "/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_reverse_current_all_flavor_volSAND/";
    auto FOLDER_ANALYSIS = "/storage/gpfs_data/neutrino/users/gi/sand-physics/production_reverse_current_volSAND/";

    auto fInput_genie =TString::Format("%sevents-in-volSAND.*.gtrac.root", FOLDER_PRODUCTION);
    auto fInput_edep =TString::Format("%sevents-in-volSAND.*.edep-sim.root", FOLDER_PRODUCTION);
    // auto fInput_digit =TString::Format("%sevents-in-volSAND.*.ecal-digit.root", FOLDER_PRODUCTION);
    auto fOutput = TString::Format("%sevents-in-volSAND.%d.to.%d.edep.analysed.root",FOLDER_ANALYSIS, file_start, file_stop);
    
    // if you have multiple files enable multiple thread pocessing
    // if(TString::Format("%s",fInput_digit).Contains("*")){
    if(fInput_edep.Contains("*")){
        LOG("I","Enabling multiple threading");
        ROOT::EnableImplicitMT();
        geo->SetMaxThreads(100);
    };

    LOG("I", "Initialize ROOT DataFrame");
    auto chain_genie = RDFUtils::InitTChain(fInput_genie, "gRooTracker", file_start, file_stop);
    auto chain_edep = RDFUtils::InitTChain(fInput_edep, "EDepSimEvents", file_start, file_stop);
    // auto chain_digit = RDFUtils::InitTChain(fInput_digit, "tDigit", file_start, file_stop); 
    // auto chain_cluster = RDFUtils::InitTChain(fInput_ecal_cluster, "tReco", file_start, file_stop); 
    
    chain_edep->AddFriend(chain_genie, "genie");
    // chain_digit->AddFriend(chain_edep, "edep");
    // chain_digit->AddFriend(chain_cluster, "cluster");
    
    auto df = RDFUtils::InitDF(chain_edep);

    // RDFUtils::PrintColumns(df);

    // throw "";

    auto dfC = RDFUtils::AddConstantsToDF(df); // add some columns with usefull constants

    auto dfGENIE = RDFUtils::GENIE::AddColumnsFromGENIE(dfC);

    auto dfEDEP = RDFUtils::EDEPSIM::NOSPILL::AddColumnsFromEDEPSIM(dfGENIE);
    
    // auto dfDigit = RDFUtils::DIGIT::AddColumnsFromDigit(dfEDEP);

    dfEDEP.Snapshot("edep_extended", fOutput.Data(), { 
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
                                                "IncomingNeutrinoPDG",
                                                "IncomingNeutrinoName",
                                                "InteractionVolume",
                                                "InteractionTarget",
                                                "FinalStateLeptonPDG",
                                                "FinalStateLeptonNames",
                                                "FinalStateLepton4Momentum",
                                                "FinalStateLeptonEmissionAngle",
                                                "NofPrimaries",
                                                "NofFinalStateChargedParticles",
                                                "FinalStateHadronicSystemPDG", 
                                                "FinalStateHadronicSystemNames", 
                                                "FinalStateHadronicSystemMomenta", 
                                                "FinalStateHadronicSystemTotal4Momentum", 
                                                "FinalStateHadronicSystemEmissionAngle", 
                                                "FinalStateHadronicSystemTotalKinE", 
                                                /*
                                                    EDEP INFO
                                                */
                                                // "PrimariesPDG",
                                                // "PrimariesTrackId",
                                                // "PrimariesP4",
                                                // "PrimariesFirstHitECAL",
                                                // "PrimariesEDepECAL",
                                                // "PrimariesEmissionAngle",
                                                /*
                                                    PREDICTIONS FOR CHANNEL antinu on H
                                                */
                                                // "ExpectedNeutrinoP4FromMuon",
                                                // "ExpectedHadronSystP3",
                                                // "ExpectedHadronSystEnergy",
                                                // "MissingTransverseMomentum",
                                                // "DoubleTransverseMomentumImbalance",

    });                                                    

    return 0;
}