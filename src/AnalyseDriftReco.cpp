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
    unsigned int files_per_jobs = 10u;
    unsigned int file_start = index * files_per_jobs;
    unsigned int file_stop = index * files_per_jobs + files_per_jobs;

    LOG("I", TString::Format("Analyze production from file %d to file %d", file_start, file_stop).Data());
    
    LOG("I", "Reading geometry");
    geo = TGeoManager::Import("/storage/gpfs_data/neutrino/users/gi/dunendggd/SAND_opt3_DRIFT1.root");

    auto FOLDER_PRODUCTION = "/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc/";
    auto FOLDER_ANALYSIS = "/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/DriftReco/Wires_cut/";
    auto FOLDER_ANALYSIS_CELLS = "/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/CompleteCells_neutron_signal/";

    auto fInput_genie =TString::Format("%sevents-in-SANDtracker.*.gtrac.root", FOLDER_PRODUCTION);
    auto fInput_edep =TString::Format("%sevents-in-SANDtracker.*.edep-sim.root", FOLDER_PRODUCTION);
    auto fInput_digit =TString::Format("%sevents-in-SANDtracker.*.ecal-digit.root", FOLDER_PRODUCTION);
    auto fInput_drift_reco = TString::Format("%sevents-in-SANDtracker.*.recostruction.NLLmethod.root", FOLDER_PRODUCTION);
    auto fOutput = TString::Format("%sevents-in-SANDtracker.%d.to.%d.drift-reco.analysed.root",FOLDER_ANALYSIS, file_start, file_stop);
    auto fOutput_cells = TString::Format("%sevents-in-SANDtracker.%d.to.%d.drift-reco.analysed.root",FOLDER_ANALYSIS_CELLS, file_start, file_stop);
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
           - KeepThisEvent==1 : vertex is in FV
           - nof_fired_wires > 70

    */
    LOG("I", "Filter events is FV volume");
    auto df_filtered = RDFUtils::Filter(dfC, "KeepThisEvent");

    auto dfGENIE = RDFUtils::GENIE::AddColumnsFromGENIE(df_filtered);
    auto dfEDEP = RDFUtils::EDEPSIM::NOSPILL::AddColumnsFromEDEPSIM(dfGENIE);
    auto dfDigit = RDFUtils::DIGIT::AddColumnsFromDigit(dfEDEP);
    auto dfReco = RDFUtils::RECO::AddColumnsFromDriftReco(dfDigit);
    
    LOG("I", "Filter events with at least 70 fired wires");
    dfReco = RDFUtils::Filter(dfReco, "pass_nof_wires_cut");

    // RDFUtils::PrintColumns(df);
    // throw "";

    LOG("I", "Writing ouput file");
    // dfReco.Snapshot("tReco_extended", fOutput.Data(), {
    //     "FileName",
    //     "CCQEonHydrogen",
    //     "edep_file_input",
    //     "digit_file_input",
    //     "edep_event_index",
    //     "IncomingNeutrinoP4",
    //     "NuDirection",
    //     /*
    //         ECAL DIGIT & RECO ____________________________
    //     */
    //     "Fired_Cells_mod",
    //     "Fired_Cells_id",
    //     "Fired_Cells_x",
    //     "Fired_Cells_y",
    //     "Fired_Cells_z",
    //     "Fired_Cells_adc1",
    //     "Fired_Cells_adc2",
    //     "Fired_Cells_tdc1",
    //     "Fired_Cells_tdc2",
    //     "Fired_Cell_true_hit1",
    //     "Fired_Cell_true_hit2",
    //     "Fired_by_primary_neutron",
    //     "Fired_by_primary_antimu",
    //     "isCellComplete",
    //     /*
    //         ECAL hit true
    //     */
    //     "Fired_Cell_true_Hit_x",
    //     "Fired_Cell_true_Hit_y",
    //     "Fired_Cell_true_Hit_z",
    //     "Fired_Cell_true_Hit_t",
    //     "True_FlightLength",
    //     /*
    //         ECAL hit reco
    //     */
    //     // "Reconstructed_HitPosition_x",
    //     // "Reconstructed_HitPosition_y",
    //     // "Reconstructed_HitPosition_z",
    //     // "Reconstructed_HitTime",
    //     // "Reconstructed_FlightLength",
    //     /*
    //         DRIFT RECO ___________________________________________
    //     */
    //     "nof_fired_wires",
    //     "Antimuon_Phi0_true",
    //     "Antimuon_x0_true",
    //     "Antimuon_pt_true",
    //     "Antimuon_p_true",
    //     "Antimuon_ptot_true",
    //     "Antimuon_dip_true",
    //     /*
    //         reconstructed antimuon
    //     */
    //     "Antimuon_Phi0_reco",
    //     "Antimuon_x0_reco",
    //     "Antimuon_pt_reco",
    //     "Antimuon_p_reco",
    //     "Antimuon_ptot_reco",
    //     "Antimuon_dip_reco",
    //     "chi2_fit_zy",
    //     "chi2_fit_xz",
    //     /*
    //         predicted neutron
    //     */
    //     "Neutrino_reconstructed_P4_GeV",
    //     "IncomingNeutrinoP4",
    //     "PredictedNeutron_P3_GeV",
    //     "FinalStateHadronicSystemTotal4Momentum",
    //     "PredictedNeutron_E_GeV",
    //     "PredictedNeutron_Beta",
    //     "PredictedNeutron_Angle",
    //     /*
    //         Define expected cell hits
    //     */
    //     "ExpectedNeutron_HitPosition_x_",
    //     "ExpectedNeutron_HitPosition_y_",
    //     "ExpectedNeutron_HitPosition_z_",
    //     "ExpectedNeutron_FlightLength_",
    //     "ExpectedNeutron_TOF_",
    //     "Expected_HitTime_",
    //     });
    // //...................................................................

    auto complete_cells_fired_by_signal_neutron = RDFUtils::DIGIT::GetInfoCellsFromSignal(dfReco);

    complete_cells_fired_by_signal_neutron.Snapshot("tReco_extended",fOutput_cells.Data(),
          {
        "FileName",
        "CCQEonHydrogen",
        //
        "neutrons_cells_mod",
        "neutrons_cells_id",
        "neutrons_cells_x",
        "neutrons_cells_y",
        "neutrons_cells_z",
        /*
            TRUE HIT
        */
        "true_hit_e",
        "true_hit_x",
        "true_hit_y",
        "true_hit_z",
        "true_hit_t",
        "true_hit_FlightLength",
    });
    return 0;
}