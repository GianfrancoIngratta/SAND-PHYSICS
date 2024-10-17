#include <iostream>
#include <string>

#include "RDataFrameUtils.h"
#include "TFile.h"

TGeoManager* geo = nullptr;

std::vector<std::string> columns_df = {
    /*
        EVENT INFO
    */
    "nof_vertices_in_spill",
    "vertex_name",
    "vertex_x",
    "vertex_y",
    "vertex_z",
    "vertex_t",
    "neutrino_flavor",
    "nuclear_target",
    "nucleon_target",
    "event_type",
    "InteractionVolume_short", 
    "vertex_antinumuCCQEonH", 
    "vertex_cherge_multiplicity",
    "Batch_number",
    "Bunch_number",
    "spill_t0",
    "spill_tmax",
    "spill_duration",
};

std::vector<std::string> columns_trj = {
    "trajectories_spill_number",
    "trajectories_id",
    "trajectories_pdg",
    "trajectories_name",
    "trajectories_starting_volume",
    "trajectories_ecal_edep",
    "trajectories_earliest_hit_ecal",
    "trajectories_latest_hit_ecal",
    "trajectories_TOF2ECAL",
};

std::vector<std::string> columns_trj_points = {
    "file_name",
    "trajectories_starting_volume",
    "trajectories_ecal_edep",
    "trajectories_TOF2ECAL",
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

std::vector<std::string> columns_cells = {
    "file_name",
    "spill_number",
    "cell_mod",
    "cell_id",
    "cell_x",
    "cell_y",
    "cell_z",
    "is_complete",
    "track_id_pmt1_hit",
    "track_id_pmt2_hit",
    "track_pdg_pmt1_hit",
    "track_pdg_pmt2_hit",
    "true_hit1",
    "true_hit2",
    "true_edep1",
    "true_edep2",
    "reco_hit",
    "reco_edep",
};

//___________________________________________________________________
int main(int argc, char* argv[]){
        unsigned int index;
    if(argc != 2){
        LOG("W", "One input number needed to run the executable");
        throw "";
    }

    index = atoi(argv[1]);
    // unsigned int files_per_jobs = 100u;
    unsigned int files_per_jobs = 10u;
    unsigned int file_start = index * files_per_jobs;
    unsigned int file_stop = index * files_per_jobs + files_per_jobs;

    LOG("I", TString::Format("Analyze production from file %d to file %d", file_start, file_stop).Data());

    LOG("I", "Reading geometry");
    geo = TGeoManager::Import("/storage/gpfs_data/neutrino/users/gi/dunendggd/SAND_opt3_DRIFT1.root");

    auto FOLDER_PRODUCTION = "/storage/gpfs_data/neutrino/SAND/PRODUCTIONS/PROD/SAND_opt3_DRIFT1_reverse_volSAND/SAND_opt3_DRIFT1_reverse_volSAND_*/";
    auto FOLDER_ECAL_DIGIT = "/storage/gpfs_data/neutrino/users/gi/sand-physics/production_spill_reverse_volSAND/SAND_opt3_DRIFT1_reverse_volSAND/ecal-digit/PROD_0/";
    auto FOLDER_ANALYSIS = "/storage/gpfs_data/neutrino/users/gi/sand-physics/production_spill_reverse_volSAND/SAND_opt3_DRIFT1_reverse_volSAND/";

    // auto fInput_genie = TString::Format("%ssand-drift-spill-revere-volSAND.*.gtrac.root",FOLDER_PRODUCTION);
    auto fInput_edep = TString::Format("%ssand-drift-spill-revere-volSAND.*.edep.root",FOLDER_PRODUCTION);
    auto fInput_ecal_digit = TString::Format("%ssand-drift-spill-revere-volSAND.*.ecal-digit.root",FOLDER_ECAL_DIGIT);

    auto fOutput = TString::Format("%ssand-drift-events.%d.to.%d.analysed.root",FOLDER_ANALYSIS, file_start, file_stop);
    auto fOutput_traj = TString::Format("%ssand-drift-events.%d.to.%d.analysed.trj.root",FOLDER_ANALYSIS, file_start, file_stop);
    auto fOutput_traj_points_spill = TString::Format("%ssand-drift-events.%d.to.%d.analysed.trj_points_spill.root",FOLDER_ANALYSIS, file_start, file_stop);
    auto fOutput_cells = TString::Format("%ssand-drift-events.%d.to.%d.analysed.cells.root",FOLDER_ANALYSIS, file_start, file_stop);

        if(fInput_edep.Contains("*")){
        LOG("I","Enabling multiple threading");
        ROOT::EnableImplicitMT();
        geo->SetMaxThreads(100);
    };

    // auto chain_genie = RDFUtils::InitTChain(fInput_genie, "gRooTracker", file_start, file_stop);
    auto chain_edep = RDFUtils::InitTChain(fInput_edep, "EDepSimEvents", file_start, file_stop);
    auto chain_ecal_digit = RDFUtils::InitTChain(fInput_ecal_digit, "tDigit", file_start, file_stop);

    // chain_ecal_digit->AddFriend(chain_genie, "genie");
    chain_ecal_digit->AddFriend(chain_edep, "edep");

    // auto df = RDFUtils::InitDF(chain_genie);
    auto df = RDFUtils::InitDF(chain_ecal_digit);

    auto df_ = RDFUtils::AddConstantsToDF(df);

    // df_ = RDFUtils::GENIE::AddColumnsFromGENIE(df_);
    df_ = RDFUtils::EDEPSIM::SPILL::AddColumnsFromEDEPSIM(df_);

    df_= RDFUtils::CreateDataFrameCells(df_);
    df_= RDFUtils::CreateDataFrameTrajectories(df_);

    // for plotting
    auto df_trajectory_points = RDFUtils::DIGIT::GetFilteredTrajectories(df_, "Trajectories");

    // RDFUtils::PrintColumns(df);

    df_.Snapshot("digit_extended", fOutput.Data(), columns_df);
    df_.Snapshot("cells", fOutput_cells.Data(), columns_cells);
    df_.Snapshot("trajectories", fOutput_traj.Data(), columns_trj);
    // for plotting
    df_trajectory_points.Snapshot("trajectories_points", fOutput_traj_points_spill.Data(), columns_trj_points);

    return 0;
}