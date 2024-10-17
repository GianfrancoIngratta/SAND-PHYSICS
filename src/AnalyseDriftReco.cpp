#include <iostream>
#include <string>

#include "RDataFrameUtils.h"
// #include "GenieUtils.h"
// #include "GeoUtils.h"
#include "TFile.h"

TGeoManager* geo = nullptr;

std::vector<std::string> output_columns_event_selection = {
    /*
        EVENT INFO
    */
    "FileName",
    "CCQEonHydrogen",
    "EventType",
    "NuDirection",
    "Interaction_vtxX",
    "Interaction_vtxY",
    "Interaction_vtxZ",
    "Interaction_vtxT",
    "InteractionVolume_short",
    "NofFinalStateChargedParticles",
    "PrimaryStateHadronicSystemTopology_name",
    "InteractionTarget",
    "candidate_signal_event",
    "nof_fired_wires",
    "IncomingNeutrinoP4",
    "FinalStateHadronicSystemTotal4Momentum", // true neutron
    "Antimuon_p_true", // true antimuon
    /*
        RECONSTRUCTED ANTIMUON
    */
    "Antimuon_reconstructed_P4",
    /*
        PREDICTED NEUTRON
    */
    "Neutrino_reconstructed_P4_GeV",
    "PredictedNeutron_P3_GeV",
    "PredictedNeutron_E_GeV",
    "PredictedNeutron_Beta",
    /*
        EVENT FIRED CELL GENERAL INFO    
    */
    "Fired_Cells_mod",
    "Fired_Cells_id",
    "Fired_Cells_x",
    "Fired_Cells_y",
    "Fired_Cells_z",
    "isCellComplete",
    "Fired_Cells_adc1",
    "Fired_Cells_tdc1",
    "Fired_Cells_adc2",
    "Fired_Cells_tdc2",
    "Fired_Cell_true_hit1",
    "Fired_Cell_true_hit2",
    "Fired_by_primary_neutron",
    "Fired_by_primary_antimu",
    /*
        TRUE NEUTRON HITS
    */
    "Fired_Cell_true_Hit_x",
    "Fired_Cell_true_Hit_y",
    "Fired_Cell_true_Hit_z",
    "Fired_Cell_true_Hit_t",
    "Fired_Cell_true_Hit_e",
    "True_FlightLength",
    /*
        PREDICTED NEUTRON HITS
    */
    "ExpectedNeutron_HitPosition_x_",
    "ExpectedNeutron_HitPosition_y_",
    "ExpectedNeutron_HitPosition_z_",
    "ExpectedNeutron_FlightLength_",
    "ExpectedNeutron_TOF_",
    "Expected_HitTime_",
    /*
        RECONSTRUCTED NEUTRON HITS
    */
    "Reconstructed_HitPosition_x",
    "Reconstructed_HitPosition_y",
    "Reconstructed_HitPosition_z",
    "Reconstructed_HitTime",
    "Reconstructed_Energy",
    "IsEarliestCell_neutron",
    "Reconstructed_FlightLength",
    /*
        CELLS WITH COINCIDENCES
    */
   "Residuals_HitTime_",
   "Residuals_HitSpace_",
   "IsCompatible",
   "earliest_compatible_cell",
   "nof_compatible_cells",
   /*
        RECONSTRUCTED NEUTRON AND NEUTRINO ENERGY
   */
    "reconstructed_neutron_KinE_MeV",
    // "reconstructed_neutrino_Energy",
};

std::vector<std::string> columns_neutron_cells = {
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
    "true_hit_is_earliest",
    "true_neutron_P4",
    "true_neutron_beta",
    "true_neutron_kinE",
    "Neutron_deviation_angle",
    /*
        PREDICTED HIT
    */
    "predicted_hit_x",
    "predicted_hit_y",
    "predicted_hit_z",
    "predicted_hit_t",
    "predicted_hit_FlightLength",
    "predicted_neutron_beta",
    /*
        RECONSTRUCTED HIT
    */
    "reconstructed_hit_x",
    "reconstructed_hit_y",
    "reconstructed_hit_z",
    "reconstructed_hit_t",
    "reconstructed_hit_e",
    "reconstructed_hit_is_earliest",
    "reconstructed_hit_FlightLength",
    "reconstructed_neutron_beta",
    "reconstructed_neutron_kinE",

    };


std::vector<std::string> columns_branch_trj = {
    "FileName",
    "CCQEonHydrogen",
    "NuDirection",
    "Interaction_vtxX",
    "Interaction_vtxY",
    "Interaction_vtxZ",
    "Interaction_vtxT",
    "PrimaryStateHadronicSystemTopology_name",
    "PredictedNeutron_P3_GeV",
    "Antimuon_reconstructed_P4",
    /*
        trj
    */
    "trackid",
    "pdg",
    "point_x",
    "point_y",
    "point_z",
    "point_px",
    "point_py",
    "point_pz",
    "process",
};

//___________________________________________________________________
int main(int argc, char* argv[]){
    unsigned int index;
    if(argc != 2){
        LOG("W", "One input number needed to run the executable");
        throw "";
    }

    index = atoi(argv[1]);
    unsigned int files_per_jobs = 100u;
    unsigned int file_start = index * files_per_jobs;
    unsigned int file_stop = index * files_per_jobs + files_per_jobs;

    LOG("I", TString::Format("Analyze production from file %d to file %d", file_start, file_stop).Data());
    
    LOG("I", "Reading geometry");
    geo = TGeoManager::Import("/storage/gpfs_data/neutrino/users/gi/dunendggd/SAND_opt3_DRIFT1.root");

    auto FOLDER_PRODUCTION = "/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc/";
    auto FOLDER_ANALYSIS = "/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/DriftReco/Wires_cut/";
    auto FOLDER_ANALYSIS_CELLS = "/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/CompleteCells_neutron_signal/";
    auto FOLDER_ANALYSIS_SIGNAL_SELECTION = "/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/Selections/Wires_cut_multi1/";

    auto fInput_genie = TString::Format("%sevents-in-SANDtracker.*.gtrac.root", FOLDER_PRODUCTION);
    auto fInput_edep = TString::Format("%sevents-in-SANDtracker.*.edep-sim.root", FOLDER_PRODUCTION);
    auto fInput_digit = TString::Format("%sevents-in-SANDtracker.*.ecal-digit.root", FOLDER_PRODUCTION);
    auto fInput_drift_reco = TString::Format("%sevents-in-SANDtracker.*.recostruction.NLLmethod.root", FOLDER_PRODUCTION);

    auto fOutput = TString::Format("%sevents-in-SANDtracker.%d.to.%d.drift-reco.analysed.root",FOLDER_ANALYSIS, file_start, file_stop);
    auto fOutput_cells = TString::Format("%sevents-in-SANDtracker.%d.to.%d.cells_neutron.root",FOLDER_ANALYSIS_CELLS, file_start, file_stop);
    auto fOutput_selected_signal = TString::Format("%sevents-in-SANDtracker.%d.to.%d.selected_signal.root",FOLDER_ANALYSIS_SIGNAL_SELECTION, file_start, file_stop);
    auto fOutput_selected_bkg = TString::Format("%sevents-in-SANDtracker.%d.to.%d.selected_bkg.root",FOLDER_ANALYSIS_SIGNAL_SELECTION, file_start, file_stop);
    auto fOutput_trj_signal = TString::Format("%sevents-in-SANDtracker.%d.to.%d.selected_signal.trj.root",FOLDER_ANALYSIS_SIGNAL_SELECTION, file_start, file_stop);
    auto fOutput_trj_bkg = TString::Format("%sevents-in-SANDtracker.%d.to.%d.selected_bkg.trj.root",FOLDER_ANALYSIS_SIGNAL_SELECTION, file_start, file_stop);

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
    LOG("I", "Filter events reconstructed muons with at least 70 fired wires");
    auto df_filtered = RDFUtils::Filter(dfC, "KeepThisEvent", true);

    auto dfGENIE = RDFUtils::GENIE::AddColumnsFromGENIE(df_filtered);
    
    LOG("I", TString::Format("Filter events is FV volume of %f mm", GeoUtils::DRIFT::FIDUCIAL_CUT).Data());
    df_filtered = RDFUtils::Filter(dfGENIE, "isInFiducialVolume==1", true);
    
    auto dfEDEP = RDFUtils::EDEPSIM::NOSPILL::AddColumnsFromEDEPSIM(df_filtered);
    auto dfDigit = RDFUtils::DIGIT::AddColumnsFromDigit(dfEDEP);
    auto dfReco = RDFUtils::RECO::AddColumnsFromDriftReco(dfDigit);

    dfReco = RDFUtils::DIGIT::GetFilteredTrajectories(dfReco, "ECALactive_trajectories");

    // PRE CUT ______________________________________________________________________________________    
    LOG("I", "Filter events with at least 70 fired wires");
    auto dfReco_wires_cut = RDFUtils::Filter(dfReco, "pass_nof_wires_cut", true);
    
    LOG("I", "Filter events with 1 charged particle in final state");
    auto dfReco_wires_cut_1cmulti = dfReco_wires_cut.Filter("NofFinalStateChargedParticles==1");

    // SELECTED SIGNAL _______________________________________________________________________________
    LOG("I", "Writing ouput file: SELECTED SIGNAL");
    LOG("i", fOutput_selected_signal.Data());
    LOG("i", fOutput_trj_signal.Data());
    auto selected_signal = dfReco_wires_cut_1cmulti.Filter("candidate_signal_event == 1");
    selected_signal.Snapshot("selected_signal", fOutput_selected_signal.Data(), output_columns_event_selection);
    selected_signal.Snapshot("selected_signal_traj", fOutput_trj_signal.Data(), columns_branch_trj);
    
    // SELECTED BKG __________________________________________________________________________________
    LOG("I", "Writing ouput file: SELECTED BKG");
    LOG("i", fOutput_selected_bkg.Data());
    LOG("i", fOutput_trj_bkg.Data());
    auto selected_bkg = dfReco_wires_cut_1cmulti.Filter("candidate_signal_event == 0");
    selected_bkg.Snapshot("selected_bkg", fOutput_selected_bkg.Data(), output_columns_event_selection);
    selected_bkg.Snapshot("selected_bkg_traj", fOutput_trj_bkg.Data(), columns_branch_trj);
    
    // NEUTRON CELLS (MC TRUTH) _______________________________________________________________________
    // LOG("I", "Writing ouput file: CELLS FIRED BY NEUTRON");
    // LOG("i", fOutput_cells.Data());
    // auto complete_cells_fired_by_signal_neutron = RDFUtils::DIGIT::GetInfoCellsFromSignal(dfReco);
    // complete_cells_fired_by_signal_neutron.Snapshot("nutron_cells", fOutput_cells.Data(), columns_neutron_cells);

    return 0;
}