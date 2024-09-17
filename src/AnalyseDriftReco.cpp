#include <iostream>
#include <string>

#include "RDataFrameUtils.h"
// #include "GenieUtils.h"
// #include "GeoUtils.h"
#include "TFile.h"

TGeoManager* geo = nullptr;

//___________________________________________________________________
int main(int argc, char* argv[]){
    unsigned int index;
    if(argc != 2){
        LOG("W", "One input number needed to run the executable");
        throw "";
    }

    index = atoi(argv[1]);
    unsigned int files_per_jobs = 1000u;
    unsigned int file_start = index * files_per_jobs;
    unsigned int file_stop = index * files_per_jobs + files_per_jobs;

    LOG("I", TString::Format("Analyze production from file %d to file %d", file_start, file_stop).Data());
    
    LOG("I", "Reading geometry");
    geo = TGeoManager::Import("/storage/gpfs_data/neutrino/users/gi/dunendggd/SAND_opt3_DRIFT1.root");

    auto FOLDER_PRODUCTION = "/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc/";
    auto FOLDER_ANALYSIS = "/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/";

    auto fInput_genie =TString::Format("%sevents-in-SANDtracker.*.gtrac.root", FOLDER_PRODUCTION);
    auto fInput_edep =TString::Format("%sevents-in-SANDtracker.*.edep-sim.root", FOLDER_PRODUCTION);
    auto fInput_drift_reco = TString::Format("%sevents-in-SANDtracker.*.recostruction.NLLmethod.root", FOLDER_PRODUCTION);
    auto fOutput = TString::Format("%sevents-in-SANDtracker.%d.to.%d.drift-reco.analysed.root",FOLDER_ANALYSIS, file_start, file_stop);
    
    // if you have multiple files enable multiple thread pocessing
    if(fInput_drift_reco.Contains("*")){
        LOG("I","Enabling multiple threading");
        ROOT::EnableImplicitMT();
        geo->SetMaxThreads(100);
    };

    LOG("I", "Initialize ROOT DataFrame");
    auto chain_genie = RDFUtils::InitTChain(fInput_genie, "gRooTracker", file_start, file_stop);
    auto chain_edep = RDFUtils::InitTChain(fInput_edep, "EDepSimEvents", file_start, file_stop);
    auto chain_drift_reco = RDFUtils::InitTChain(fInput_drift_reco, "tReco", file_start, file_stop); 
    
    chain_drift_reco->AddFriend(chain_genie, "genie");
    chain_drift_reco->AddFriend(chain_edep, "edep");
    auto df = RDFUtils::InitDF(chain_drift_reco);

    auto dfC = RDFUtils::AddConstantsToDF(df); // add some columns with usefull constants
    auto dfReco = RDFUtils::RECO::AddColumnsFromDriftReco(dfC);
    
    // RDFUtils::PrintColumns(df);
    // throw "";

    LOG("I", "Writing ouput file");
    dfReco.Snapshot("tReco_extended", fOutput.Data(), {
        "edep_file_input",
        "digit_file_input",
        "edep_event_index",
        // "NofFiredWires",
        "pt_true",
        "pt_reco",
        "p_true",
        "p_reco",
        "ptot_true",
        "ptot_reco",
        "true_helix_dip_",
        "reco_helix_dip_",
        });

    return 0;
}