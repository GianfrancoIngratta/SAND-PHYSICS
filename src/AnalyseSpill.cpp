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
    "nof_vertices_in_spill",
    // add vertex time to calculate tof to ecal
    "trajectories_id",
    "trajectories_pdg",
    "trajectories_name",
    "trajectories_ecal_edep",
    "trajectories_earliest_hit_ecal",
    "trajectories_latest_hit_ecal",
};

std::vector<std::string> columns_cells = {
    "nof_vertices_in_spill",
    "cell_x",
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
    unsigned int files_per_jobs = 2u;
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
    df_ = RDFUtils::DIGIT::SPILL::AddColumnsFromDigit(df_);

    df_= RDFUtils::CreateDataFrameTrajectories(df_);
    df_= RDFUtils::CreateDataFrameCells(df_);

    // RDFUtils::PrintColumns(df);

    // df_.Snapshot("digit_extended", fOutput.Data(), columns_df);
    // df_.Snapshot("trajectories", fOutput_traj.Data(), columns_trj);
    df_.Snapshot("cells", fOutput_traj.Data(), columns_cells);

    return 0;
}