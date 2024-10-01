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
    auto fInput_digit =TString::Format("%sevents-in-SANDtracker.*.ecal-digit.root", FOLDER_PRODUCTION);
    auto fInput_drift_reco = TString::Format("%sevents-in-SANDtracker.*.recostruction.NLLmethod.root", FOLDER_PRODUCTION);
    auto fOutput = TString::Format("%sevents-in-SANDtracker.%d.to.%d.drift-reco.analysed.root",FOLDER_ANALYSIS, file_start, file_stop);
    // auto fOutput = TString::Format("%sevents-in-SANDtracker.%d.to.%d.all.analysed.root",FOLDER_ANALYSIS, file_start, file_stop);
    
    // if you have multiple files enable multiple thread pocessing
    if(fInput_drift_reco.Contains("*")){
        LOG("I","Enabling multiple threading");
        ROOT::EnableImplicitMT();
        geo->SetMaxThreads(100);
    };

    LOG("I", "Initialize ROOT DataFrame");
    auto chain_genie = RDFUtils::InitTChain(fInput_genie, "gRooTracker", file_start, file_stop);
    auto chain_edep = RDFUtils::InitTChain(fInput_edep, "EDepSimEvents", file_start, file_stop);
    auto chain_digit = RDFUtils::InitTChain(fInput_digit, "tDigit", file_start, file_stop);
    auto chain_drift_reco = RDFUtils::InitTChain(fInput_drift_reco, "tReco", file_start, file_stop); 
    
    chain_drift_reco->AddFriend(chain_genie, "genie");
    chain_drift_reco->AddFriend(chain_edep, "edep");
    chain_drift_reco->AddFriend(chain_digit, "ecal-digit");
    
    auto df = RDFUtils::InitDF(chain_drift_reco);
    auto dfC = RDFUtils::AddConstantsToDF(df); // add some columns with usefull constants
    /*
        APPLY FILTERS:
           - KeepThisEvent==1 : vertex is in FV and muon track has enough hits to be reconstructed
    */
    auto df_filtered = RDFUtils::Filter(dfC, "KeepThisEvent==1");

    auto dfGENIE = RDFUtils::GENIE::AddColumnsFromGENIE(df_filtered);
    auto dfEDEP = RDFUtils::EDEPSIM::NOSPILL::AddColumnsFromEDEPSIM(dfGENIE);
    auto dfDigit = RDFUtils::DIGIT::AddColumnsFromDigit(dfEDEP);
    auto dfReco = RDFUtils::RECO::AddColumnsFromDriftReco(dfDigit);

    // RDFUtils::PrintColumns(df);
    // throw "";

    LOG("I", "Writing ouput file");
    dfReco.Snapshot("tReco_extended", fOutput.Data(), {
        "FileName",
        "CCQEonHydrogen",
        "edep_file_input",
        "digit_file_input",
        "edep_event_index",
        "IncomingNeutrinoP4",
        "NuDirection",
        /*
            ECAL DIGIT & RECO
        */
        // "NofEventFiredModules",
        // "EventFiredModules",
        // "Fired_Cells_mod",
        // "Fired_Cells_id",
        // "Fired_Cells_x",
        // "Fired_Cells_y",
        // "Fired_Cells_z",
        // "Fired_Cells_adc1",
        // "Fired_Cells_adc2",
        // "Fired_Cells_tdc1",
        // "Fired_Cells_tdc2",
        // "who_produced_tdc1",
        // "who_produced_tdc2",
        // "Fired_Cell_true_hit1",
        // "Fired_Cell_true_hit2",
        // "isCellComplete",
        // "Cell_Reconstructed_hit",
        // "ExpectedNeutronHit",
        "nof_fired_wires",
        /*
            DRIFT RECO
        */
        "Antimuon_Phi0_true",
        "Antimuon_x0_true",
        "Antimuon_pt_true",
        "Antimuon_p_true",
        "Antimuon_ptot_true",
        "Antimuon_dip_true",
        //
        "Antimuon_Phi0_reco",
        "Antimuon_x0_reco",
        "Antimuon_pt_reco",
        "Antimuon_p_reco",
        "Antimuon_ptot_reco",
        "Antimuon_dip_reco",
        "chi2_fit_zy",
        "chi2_fit_xz",
        //
        "Neutrino_reconstructed_P4_GeV",
        "IncomingNeutrinoP4",
        "PredictedNeutron_P3_GeV",
        "FinalStateHadronicSystemTotal4Momentum",
        });

    return 0;
}