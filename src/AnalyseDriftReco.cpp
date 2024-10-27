#include <iostream>
#include <string>
#include <TLegend.h>

#include "RDataFrameUtils.h"
#include "TFile.h"

TGeoManager* geo = nullptr;
/***
 * SCALE: mass ratio carbon_in_graphite / carbon_in_plastic
*/
const double scale_factor = 4.084;

std::vector<std::string> columns_preselection = {
    /*
        EVENT INFO
    */
    "FileName",
    "CCQEonHydrogen",
    "EventType",
    "IncomingNeutrinoP4",
    "InteractionVolume_short",  
    "InteractionTarget",
    "PrimaryStateHadronicSystemTopology_name",
    "Antimuon_p_true", // true antimuon
    "FinalStateHadronicSystemTotal4Momentum", // true neutron
    "Interaction_vtxX",
    "Interaction_vtxY",
    "Interaction_vtxZ",
    "Interaction_vtxT",
    "NofFinalStateChargedParticles",
    "nof_fired_wires",
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
};

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
    "AreTDCsConsistent",
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

const char* hist_names[] = {"total",
                            "CCQE_on_H",
                            "others_on_H",
                            "on_C_plastic",
                            "on_C_graphite",
                            "others",
                            };

TH1I* report_histos[sizeof(hist_names)/sizeof(void*)];

// const char* spectra_names[] = {
//                             "interacting_neutrino_true_energy",
//                             "antimuon_true_energy",
//                             "final_state_hadronic_kinetic_energy"

// }

// TH1D* spectra[sizeof(spectra_names)/sizeof(void*)];

// struct HIST_INFO
// {
//     const char* name;   
//     const char* title;  
//     unsigned int nof_bins;
//     double low_bin;
//     double sup_bin;
//     ROOT::RDF::TH1DModel model  = {name, title, nof_bins, low_bin, sup_bin};
// };

// HIST_INFO interacting_neutrino_true_E = {"interacting neutrino energy", "interacting neutrino energy; [GeV]", 100u, 0., 8.};
// HIST_INFO antimuon_true_E = {"antimuon true E", "antimuon true energy; [GeV]", 100u, 0., 8.};
// HIST_INFO fs_hadron_system_Kin = {"fs_hadron_system_Kin", "final state hadronic system true energy; [GeV]", 100u, 0., 2.};

int nof_report_calls = 1;

void Report(ROOT::RDF::RNode& input_df, std::ofstream& stream, const char* report_stage){
    auto total_in_FV = input_df.Count().GetValue();
    auto nof_CCQE_on_H_in_FV = input_df.Filter("CCQEonHydrogen==1").Count().GetValue();
    auto nof_others_on_H_in_FV = input_df.Filter("CCQEonHydrogen==0").Filter(TString::Format("InteractionTarget == \"%s\"", "proton").Data()).Count().GetValue();
    auto nof_on_C_plastic_in_FV = input_df.Filter(TString::Format("InteractionTarget == \"%s\"", "C12").Data()).Filter(TString::Format("InteractionVolume_short == \"%s\"", "C3H6_Target").Data()).Count().GetValue();
    auto nof_on_C_grphite_in_FV = input_df.Filter(TString::Format("InteractionTarget == \"%s\"", "C12").Data()).Filter(TString::Format("InteractionVolume_short == \"%s\"", "C_Target").Data()).Count().GetValue();
    auto others = total_in_FV - nof_CCQE_on_H_in_FV - nof_others_on_H_in_FV - nof_on_C_plastic_in_FV - nof_on_C_grphite_in_FV;
    stream << "############### " << report_stage << " ###############\n";
    stream << "TOTAL                 : " << total_in_FV << "\n";
    stream << "TOTAL H               : " << nof_CCQE_on_H_in_FV + nof_others_on_H_in_FV << "\n";
    stream << "CCQE ON H             : " << nof_CCQE_on_H_in_FV << "\n";
    stream << "OTHERS ON H           : " << nof_others_on_H_in_FV << "\n";
    stream << "ANY INT ON C PLASTIC  : " << nof_on_C_plastic_in_FV << "\n";
    stream << "ANY INT ON C GRAPHITE : " << nof_on_C_grphite_in_FV << "\n";
    stream << "OTHERS                : " << others << "\n\n";
}

void Report(ROOT::RDF::RNode& input_df, const char* report_stage){
    auto total_in_FV = input_df.Count().GetValue();
    auto nof_CCQE_on_H_in_FV = input_df.Filter("CCQEonHydrogen==1").Count().GetValue();
    auto nof_others_on_H_in_FV = input_df.Filter("CCQEonHydrogen==0").Filter(TString::Format("InteractionTarget == \"%s\"", "proton").Data()).Count().GetValue();
    auto nof_on_C_plastic_in_FV = input_df.Filter(TString::Format("InteractionTarget == \"%s\"", "C12").Data()).Filter(TString::Format("InteractionVolume_short == \"%s\"", "C3H6_Target").Data()).Count().GetValue();
    auto nof_on_C_grphite_in_FV = input_df.Filter(TString::Format("InteractionTarget == \"%s\"", "C12").Data()).Filter(TString::Format("InteractionVolume_short == \"%s\"", "C_Target").Data()).Count().GetValue();
    auto others = total_in_FV - nof_CCQE_on_H_in_FV - nof_others_on_H_in_FV - nof_on_C_plastic_in_FV - nof_on_C_grphite_in_FV;

    report_histos[0]->SetBinContent(nof_report_calls, total_in_FV);
    report_histos[0]->GetXaxis()->SetBinLabel(nof_report_calls, report_stage);
    report_histos[1]->SetBinContent(nof_report_calls, nof_CCQE_on_H_in_FV);
    report_histos[1]->GetXaxis()->SetBinLabel(nof_report_calls, report_stage);
    report_histos[2]->SetBinContent(nof_report_calls, nof_others_on_H_in_FV);
    report_histos[2]->GetXaxis()->SetBinLabel(nof_report_calls, report_stage);
    report_histos[3]->SetBinContent(nof_report_calls, nof_on_C_plastic_in_FV);
    report_histos[3]->GetXaxis()->SetBinLabel(nof_report_calls, report_stage);
    report_histos[4]->SetBinContent(nof_report_calls, nof_on_C_grphite_in_FV);
    report_histos[4]->GetXaxis()->SetBinLabel(nof_report_calls, report_stage);
    report_histos[5]->SetBinContent(nof_report_calls, others);
    report_histos[5]->GetXaxis()->SetBinLabel(nof_report_calls, report_stage);
    nof_report_calls++;
}

//___________________________________________________________________
int main(int argc, char* argv[]){
    
    // unsigned int production; // 0 to 42
    unsigned int run_start; // 0 to 4300
    
    if(argc != 2){
        LOG("W", "One input number needed to run the executable");
        throw "";
    }

    // production = atoi(argv[1]);
    run_start = atoi(argv[1]);
    unsigned int run_per_production = 10u;
    uint production = run_start / 100;
    unsigned int file_start = run_start * run_per_production;
    unsigned int file_stop = file_start + run_per_production;

    LOG("I", TString::Format("Analyze production %d from file %d to file %d", production, file_start, file_stop).Data());
    
    LOG("I", "Reading geometry");
    geo = TGeoManager::Import("/storage/gpfs_data/neutrino/users/gi/dunendggd/SAND_opt3_DRIFT1.root");

    auto FOLDER_PRODUCTION = TString::Format("/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc/production_%04d/", production);
    // auto FOLDER_ECAL_CELLS = "/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/CompleteCells_neutron_signal/";
    // auto FOLDER_PRESELECTION = "/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/antinumu_CCQE_on_H_like/Preselection/";
    // auto FOLDER_SELECTED_SIGNAL = "/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/antinumu_CCQE_on_H_like/";

    /***
     * INPUT:
     * 
    */
    auto fInput_genie = TString::Format("%srun_*/events-in-SANDtracker.*.gtrac.root", FOLDER_PRODUCTION.Data());
    auto fInput_edep = TString::Format("%srun_*/events-in-SANDtracker.*.edep-sim.root", FOLDER_PRODUCTION.Data());
    auto fInput_digit = TString::Format("%srun_*/events-in-SANDtracker.*.ecal-digit.root", FOLDER_PRODUCTION.Data());
    auto fInput_drift_reco = TString::Format("%srun_*/events-in-SANDtracker.*.recostruction.NLLmethod.root", FOLDER_PRODUCTION.Data());
    /***
     * OUTPUT:
     * 
    */
    auto fOutput_preselection = TString::Format("%sevents-in-SANDtracker.%d.to.%d.preselection.root", FOLDER_PRODUCTION.Data(), file_start, file_stop-1);
    auto fOutput_selected_signal = TString::Format("%sevents-in-SANDtracker.%d.to.%d.selected_signal.root", FOLDER_PRODUCTION.Data(), file_start, file_stop-1);
    auto fOutput_selected_bkg = TString::Format("%sevents-in-SANDtracker.%d.to.%d.selected_bkg.root", FOLDER_PRODUCTION.Data(), file_start, file_stop-1);
    auto fOutput_preunfold = TString::Format("%spreunfold/events-in-SANDtracker.%d.to.%d.preunfold.root", FOLDER_PRODUCTION.Data(), file_start, file_stop-1);
    auto fOutput_report = TString::Format("%sevent_selection/events-in-SANDtracker.%d.to.%d.report.root", FOLDER_PRODUCTION.Data(), file_start, file_stop-1);
    auto fOutput_report_txt = TString::Format("%sevent_selection/events-in-SANDtracker.%d.to.%d.report.txt", FOLDER_PRODUCTION.Data(), file_start, file_stop-1);

    // auto fOutput_cells = TString::Format("%sevents-in-SANDtracker.%d.to.%d.cells_neutron.root",FOLDER_ECAL_CELLS, file_start, file_stop);
    // auto fOutput_preselection = TString::Format("%sevents-in-SANDtracker.%d.to.%d.preselection.root", FOLDER_PRESELECTION, file_start, file_stop);
    // auto fOutput_selected_signal = TString::Format("%sevents-in-SANDtracker.%d.to.%d.selected_signal.root", FOLDER_SELECTED_SIGNAL, file_start, file_stop);
    // auto fOutput_selected_bkg = TString::Format("%sevents-in-SANDtracker.%d.to.%d.selected_bkg.root",FOLDER_SELECTED_SIGNAL, file_start, file_stop);
    // auto fOutput_trj_signal = TString::Format("%sevents-in-SANDtracker.%d.to.%d.selected_signal.trj.root",FOLDER_SELECTED_SIGNAL, file_start, file_stop);
    // auto fOutput_trj_bkg = TString::Format("%sevents-in-SANDtracker.%d.to.%d.selected_bkg.trj.root",FOLDER_SELECTED_SIGNAL, file_start, file_stop);
    // auto fOutput_report = TString::Format("%sreports/events-in-SANDtracker.%d.to.%d.report.root", FOLDER_SELECTED_SIGNAL, file_start, file_stop);
    // auto fOutput_report_txt = TString::Format("%sreports/events-in-SANDtracker.%d.to.%d.report.txt", FOLDER_SELECTED_SIGNAL, file_start, file_stop);
    // auto fOutput_preunfold = TString::Format("%sreports/events-in-SANDtracker.%d.to.%d.preunfold.root", FOLDER_SELECTED_SIGNAL, file_start, file_stop);
    
    TFile* report_file = new TFile(fOutput_report.Data(), "RECREATE");
    std::ofstream report_file_txt(fOutput_report_txt.Data());
    
    report_file->cd();
    for (size_t i = 0; i < 6; i++)
    {
        report_histos[i] = new TH1I(hist_names[i], hist_names[i], 6, 0, 6);
    }
    
    // std::ofstream report_file(fOutput_report.Data());
    // report_file.open();
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
        FILTER:
           - KeepThisEvent==1 : vertex is in FV has enough minimum hits to be reconstructed
    */
    auto dfC_initial_filter = RDFUtils::Filter(dfC, "KeepThisEvent==1",true);
    
    LOG("I", "Adding columns ... ");

    LOG("i", "from GENIE");
    auto dfGENIE = RDFUtils::GENIE::AddColumnsFromGENIE(dfC_initial_filter);

    LOG("i", "from edepsim");
    auto dfEDEP = RDFUtils::EDEPSIM::NOSPILL::AddColumnsFromEDEPSIM(dfGENIE);

    LOG("i", "from digit file");
    auto dfDIGIT = RDFUtils::DIGIT::AddColumnsFromDigit(dfEDEP);

    LOG("i", "from drift reco file");
    auto dfReco = RDFUtils::RECO::AddColumnsFromDriftReco(dfDIGIT);
    dfReco = dfReco.Filter("KeepThisEvent==1");
    dfReco = RDFUtils::DIGIT::GetFilteredTrajectories(dfReco, "ECALactive_trajectories");
    
    /***
     * STAGE: CUT FIDUCIAL VOLUME _________________________________________________________________________________________________________________________________________________________________________
     */
    LOG("I", TString::Format("CUT: Filter events is FV volume of %f mm", GeoUtils::DRIFT::FIDUCIAL_CUT).Data());
    
    auto df_filtered = RDFUtils::Filter(dfReco, "isInFiducialVolume==1", true);
    
    Report(df_filtered, "FIDUCIAL VOLUME CUT");
    Report(df_filtered, report_file_txt, "FIDUCIAL VOLUME CUT");
    // df_filtered.Snapshot("preselection", fOutput_preselection.Data(), columns_preselection);

    /***
     * STAGE: CUT 70 FIRED WIRES _________________________________________________________________________________________________________________________________________________________________________
    */
    LOG("I", "CUT: Filter events with at least 70 fired wires"); // 70% events
    
    auto dfReco_wires_cut = RDFUtils::Filter(df_filtered, "pass_nof_wires_cut", true);
    
    Report(dfReco_wires_cut, "MINIMUM NOF WIRES ON RECONSTRUCTED MU+");
    Report(dfReco_wires_cut, report_file_txt, "MINIMUM NOF WIRES ON RECONSTRUCTED MU+");
    
    /***
     * STAGE: 1 FINAL STATE CHARGED PARTICLE _________________________________________________________________________________________________________________________________________________________________________
    */
    LOG("I", "CUT: Filter events with 1 charged particle in final state");
    
    auto dfReco_wires_cut_1cmulti = RDFUtils::Filter(dfReco_wires_cut, "NofFinalStateChargedParticles==1", true);
    
    Report(dfReco_wires_cut_1cmulti, "CHARGE MULTIPLICITY");
    Report(dfReco_wires_cut_1cmulti, report_file_txt, "CHARGE MULTIPLICITY");

    // dfReco_wires_cut_1cmulti.Snapshot("test", "test/test.root", output_columns_event_selection);
    // dfReco_wires_cut_1cmulti.Snapshot("test_trj", "test/test_trj.root", columns_branch_trj);

    /***
     * STAGE: SELECTED SIGNAL _________________________________________________________________________________________________________________________________________________________________________
    */
    LOG("I", "Writing ouput file: SELECTED SIGNAL");
    // LOG("i", fOutput_selected_signal.Data());
    // LOG("i", fOutput_trj_signal.Data());
    
    auto selected_signal = RDFUtils::Filter(dfReco_wires_cut_1cmulti, "candidate_signal_event == 1", true);
    auto selected_signal_on_graphite = RDFUtils::Filter(selected_signal, TString::Format("InteractionVolume_short == \"%s\"", "C_Target").Data(), true);
    auto selected_signal_on_plastic = RDFUtils::Filter(selected_signal, TString::Format("InteractionVolume_short == \"%s\"", "C3H6_Target").Data(), true);
    Report(selected_signal, "SIGNAL SELECTION");
    Report(selected_signal_on_graphite, "SIGNAL SELECTION [ON GRAPHITE]");
    Report(selected_signal_on_plastic, "SIGNAL SELECTION [ON PLASTIC]");

    Report(selected_signal, report_file_txt, "SIGNAL SELECTION");
    Report(selected_signal_on_graphite, report_file_txt, "SIGNAL SELECTION [ON GRAPHITE]");
    Report(selected_signal_on_plastic, report_file_txt, "SIGNAL SELECTION [ON PLASTIC]");
    // selected_signal.Snapshot("selected_signal", fOutput_selected_signal.Data(), output_columns_event_selection);
    // selected_signal.Snapshot("selected_signal_traj", fOutput_trj_signal.Data(), columns_branch_trj);
    
    // // EXCLUDED FROM SELECTION ___________________________________________________________________________________
    // LOG("I", "Writing ouput file: SELECTED BKG");
    // LOG("i", fOutput_selected_bkg.Data());
    // LOG("i", fOutput_trj_bkg.Data());
    // auto selected_bkg = dfReco_wires_cut_1cmulti.Filter("candidate_signal_event == 0");
    // selected_bkg.Snapshot("selected_bkg", fOutput_selected_bkg.Data(), output_columns_event_selection);
    // selected_bkg.Snapshot("selected_bkg_traj", fOutput_trj_bkg.Data(), columns_branch_trj);
    
    // // NEUTRON CELLS (MC TRUTH) _______________________________________________________________________
    // LOG("I", "Writing ouput file: CELLS FIRED BY NEUTRON");
    // LOG("i", fOutput_cells.Data());
    // auto complete_cells_fired_by_signal_neutron = RDFUtils::DIGIT::GetInfoCellsFromSignal(dfReco);
    // complete_cells_fired_by_signal_neutron.Snapshot("neutron_cells", fOutput_cells.Data(), columns_neutron_cells);
    
    /***
     * PLOTS:________________________________________________________________________________________________________________________________________________________________________
     */

    /***
     * STAGE: 1 _____________________________________________________________________________________________________________________________________________________________________________________________
     */
    LOG("I", "Defining HIST STAGE 1");
    // XSEC
    auto total_xsex_vs_neutrino_energy                          = df_filtered.Profile1D({"total_xsex_vs_neutrino_energy", "total_xsex_vs_neutrino_energy", 1000u, 0., 50.},"IncomingNeutrino_energy","EvtXSec");

    auto ccqe_on_H_xsec_vs_neutrino_energy                      = df_filtered.Filter("CCQEonHydrogen==1").Profile1D({"ccqe_on_H_xsec_vs_neutrino_energy", "ccqe_on_H_xsec_vs_neutrino_energy", 1000u, 0., 50.},"IncomingNeutrino_energy","EvtXSec");

    auto ccqe_on_H_xsec_vs_hadronic_K_profile                   = df_filtered.Filter("CCQEonHydrogen==1").Profile1D({"ccqe_on_H_xsec_vs_hadronic_K", "ccqe_on_H_xsec; final state hadronic kinetic E [GeV] ", 100u, 0., 20.},"FinalStateHadronicSystemTotalKinE","EvtXSec");

    auto ccqe_on_H_xsec_vs_hadronic_K                           = df_filtered.Filter("CCQEonHydrogen==1").Graph("FinalStateHadronicSystemTotalKinE","EvtXSec");

    // TRUE
    auto interacting_neutrino_true_E                            = df_filtered.Histo1D({"interacting_nu_stage1", "interacting neutrino true Energy;[GeV]", 100u, 0., 8.}, "IncomingNeutrino_energy");

    auto antimuon_true_E                                        = df_filtered.Histo1D({"antimuon_stage1", "antimuon true Energy;[GeV]", 100u, 0., 8.}, "FinalStateLepton_energy");

    auto fs_hadron_syst_Kin                                     = df_filtered.Histo1D({"fs_hadron_syst_Kin_stage1", "final state hadron system true K;[GeV]", 100u, 0., 2.}, "PrimaryStateHadronicSystemTotalKinE");
    
    // /***
    //  * STAGE: 2 _____________________________________________________________________________________________________________________________________________________________________________________________
    //  */
    LOG("I", "Defining HIST STAGE 2");
    // TRUE
    auto interacting_neutrino_true_E_cut_wires                  = dfReco_wires_cut.Histo1D({"interacting_nu_stage2", "interacting neutrino true Energy;[GeV]", 100u, 0., 8.}, "IncomingNeutrino_energy");

    auto antimuon_true_E_cut_wires                              = dfReco_wires_cut.Histo1D({"antimuon_stage2", "antimuon true Energy;[GeV]", 100u, 0., 8.}, "FinalStateLepton_energy");

    auto fs_hadron_syst_Kin_cut_wires                           = dfReco_wires_cut.Histo1D({"fs_hadron_syst_Kin_stage2", "final state hadron system true K;[GeV]", 100u, 0., 2.}, "PrimaryStateHadronicSystemTotalKinE");
    
    // /***
    //  * STAGE: 3 _____________________________________________________________________________________________________________________________________________________________________________________________
    //  */
    // TRUE
    LOG("I", "Defining HIST STAGE 3");
    auto interacting_neutrino_true_E_charge_multi               = dfReco_wires_cut_1cmulti.Histo1D({"interacting_nu_stage3", "interacting neutrino true Energy;[GeV]", 100u, 0., 8.}, "IncomingNeutrino_energy");

    auto antimuon_true_E_charge_multi                           = dfReco_wires_cut_1cmulti.Histo1D({"antimuon_stage3", "antimuon true Energy ;[GeV]", 100u, 0., 8.}, "FinalStateLepton_energy");

    auto fs_hadron_syst_Kin_charge_multi                        = dfReco_wires_cut_1cmulti.Histo1D({"fs_hadron_syst_Kin_stage3", "final state hadron system true K;[GeV]", 100u, 0., 2.}, "PrimaryStateHadronicSystemTotalKinE");

    // RECO
    auto antimuon_reco_E_charge_multi                           = dfReco_wires_cut_1cmulti.Histo1D({"antimuon_reco_stage3", " #mu+ reco Energy;[GeV]", 100u, 0., 8.}, "Antimuon_reconstructed_energy");

    /***
     * STAGE: 4 _____________________________________________________________________________________________________________________________________________________________________________________________
     */
    // TRUE
    LOG("I", "Defining HIST STAGE 4");
    auto interacting_neutrino_true_E_selection                  = selected_signal.Histo1D({"interacting_nu_stage4", "interacting neutrino true Energy;[GeV]", 100u, 0., 8.}, "IncomingNeutrino_energy");
    
    auto fs_hadron_syst_Kin_selection                           = selected_signal.Histo1D({"fs_hadron_syst_Kin_stage4", "final state hadron system true K;[GeV]", 100u, 0., 2.}, "PrimaryStateHadronicSystemTotalKinE");
    
    auto antimuon_true_E_selection                              = selected_signal.Histo1D({"antimuon_stage4", "antimuon true Energy ;[GeV]", 100u, 0., 8.}, "FinalStateLepton_energy");

    /***
     * STAGE: 4 _____________________________________________________________________________________________________________________________________________________________________________________________
     */

    auto interacting_antimuon_reco_E_selection                  = selected_signal.Histo1D({"antimuon_reco_stage4", " #mu+ reco Energy;[GeV]", 100u, 0., 8.}, "Antimuon_reconstructed_energy");
    
    auto interacting_neutrino_reco_E_selection                  = selected_signal.Histo1D({"interacting_nu_reco_stage4", "interacting neutrino from #mu+ reco Energy;[GeV]", 100u, 0., 8.}, "Neutrino_reconstructed_energy_GeV");
    
    auto interacting_neutrino_reco_E_selected_signal_on_graphite = selected_signal_on_graphite.Histo1D({"interacting_nu_reco_graphite_stage5a", "interacting neutrino reco E (GRAPHITE); [GeV]", 100u, 0., 8.}, "Neutrino_reconstructed_energy_GeV");
    
    auto interacting_neutrino_reco_E_selected_signal_on_plastic = selected_signal_on_plastic.Histo1D({"interacting_nu_reco_plastic_stage5b", "interacting neutrino reco E (PLASTIC); [GeV]", 100u, 0., 8.}, "Neutrino_reconstructed_energy_GeV");
    
    /***
     * ANALYSIS: SOLID HYDROGEN TECNIQUE
     */

    TH1D* interacting_neutrino_reco_E_selected_signal_on_H = interacting_neutrino_reco_E_selected_signal_on_plastic.GetPtr();
    interacting_neutrino_reco_E_selected_signal_on_H->Sumw2();
    interacting_neutrino_reco_E_selected_signal_on_H->Add(interacting_neutrino_reco_E_selected_signal_on_graphite.GetPtr(), - scale_factor);
    interacting_neutrino_reco_E_selected_signal_on_H->SetName("interacting_neutrino_reco_E_selected_signal_on_H");
    interacting_neutrino_reco_E_selected_signal_on_H->SetTitle("interacting_neutrino_reco_E_selected_signal_on_H");
    // auto interacting_neutrino_reco_E_graphite_selection         = selected_signal
    //                                                                     .Filter(TString::Format("InteractionVolume_short == \"%s\"", "C_Target").Data())
    //                                                                     .Histo1D({"interacting_nu_reco_graphite_stage4", "interacting neutrino reco E (GRAPHITE); [GeV]", 100u, 0., 8.}, "Neutrino_reconstructed_energy_GeV");
    
    // auto interacting_neutrino_reco_E_plastic_selection           = selected_signal
    //                                                                     .Filter(TString::Format("InteractionVolume_short == \"%s\"", "C3H6_Target").Data())
    //                                                                     .Histo1D({"interacting_nu_reco_plastic_stage4", "interacting neutrino reco E (PLASTIC); [GeV]", 100u, 0., 8.}, "Neutrino_reconstructed_energy_GeV");
    
    auto interacting_neutrino_reco_E_carbon_selection           = selected_signal
                                                                        .Filter(TString::Format("InteractionTarget == \"%s\"", "C12").Data())
                                                                        .Histo1D({"interacting_nu_reco_int_on_C12_stage4", "interacting neutrino reco E (interaction on Carbon); [GeV]", 100u, 0., 8.}, "Neutrino_reconstructed_energy_GeV");
    
    auto interacting_neutrino_reco_E_proton_selection           = selected_signal
                                                                        .Filter(TString::Format("InteractionTarget == \"%s\"", "proton").Data())
                                                                        .Histo1D({"interacting_nu_reco_int_on_proton_stage4", "interacting neutrino reco E (interaction on proton); [GeV]", 100u, 0., 8.}, "Neutrino_reconstructed_energy_GeV");

    auto interacting_neutrino_reco_E_plastic_proton_selection   = selected_signal
                                                                        .Filter(TString::Format("InteractionVolume_short == \"%s\"", "C3H6_Target").Data())
                                                                        .Filter(TString::Format("InteractionTarget == \"%s\"", "proton").Data())
                                                                        .Histo1D({"interacting_neutrino_reco_E_plastic_proton_selection", "interacting neutrino reco E (interaction on proton of plastic); [GeV]", 100u, 0., 8.}, "Neutrino_reconstructed_energy_GeV");

    auto interacting_neutrino_reco_E_plastic_carbon_selection   = selected_signal
                                                                        .Filter(TString::Format("InteractionVolume_short == \"%s\"", "C_Target").Data())
                                                                        .Filter(TString::Format("InteractionTarget == \"%s\"", "C12").Data())
                                                                        .Histo1D({"interacting_neutrino_reco_E_plastic_carbon_selection", "interacting neutrino reco E (interaction on C12 of plastic); [GeV]", 100u, 0., 8.}, "Neutrino_reconstructed_energy_GeV");

    TProfile* profiles[] = {
        total_xsex_vs_neutrino_energy.GetPtr(),
        ccqe_on_H_xsec_vs_neutrino_energy.GetPtr(),
        ccqe_on_H_xsec_vs_hadronic_K_profile.GetPtr()
    };

    TGraph* graphs[] = {
        ccqe_on_H_xsec_vs_hadronic_K.GetPtr()
    };

    TH1D* histos[] = {
        // TRUE
        interacting_neutrino_true_E.GetPtr(),
        interacting_neutrino_true_E_cut_wires.GetPtr(),
        interacting_neutrino_true_E_charge_multi.GetPtr(),
        interacting_neutrino_true_E_selection.GetPtr(),
        antimuon_true_E.GetPtr(),
        antimuon_true_E_cut_wires.GetPtr(),
        antimuon_true_E_charge_multi.GetPtr(),
        antimuon_true_E_selection.GetPtr(),
        fs_hadron_syst_Kin.GetPtr(),
        fs_hadron_syst_Kin_cut_wires.GetPtr(),
        fs_hadron_syst_Kin_charge_multi.GetPtr(),
        fs_hadron_syst_Kin_selection.GetPtr(),
        // RECO
        antimuon_reco_E_charge_multi.GetPtr(),
        interacting_antimuon_reco_E_selection.GetPtr(),
        interacting_neutrino_reco_E_selection.GetPtr(),
        // STAGE 4
        interacting_neutrino_reco_E_carbon_selection.GetPtr(),
        interacting_neutrino_reco_E_proton_selection.GetPtr(),
        interacting_neutrino_reco_E_plastic_proton_selection.GetPtr(),
        interacting_neutrino_reco_E_plastic_carbon_selection.GetPtr(),
        //  STAGE 5
        interacting_neutrino_reco_E_selected_signal_on_graphite.GetPtr(),
        interacting_neutrino_reco_E_selected_signal_on_plastic.GetPtr(),
        interacting_neutrino_reco_E_selected_signal_on_H,
    };
    
    LOG("I", "WRITING HISTOS ON FILE");
    for (auto& h : histos)
    {
        report_file->cd();
        h->Write();
    }
    
    for (auto& p : profiles)
    {
        report_file->cd();
        p->Write();
    }
    
    for (auto& g : graphs)
    {
        report_file->cd();
        g->Write();
    }

    report_file->cd();
    report_file->Write();

    /***
     * UNFOLD:
     */
    df_filtered.Snapshot("preunfold", fOutput_preunfold.Data(), {
                                                                    "CCQEonHydrogen",
                                                                    "NofFinalStateChargedParticles",
                                                                    "pass_nof_wires_cut",
                                                                    "candidate_signal_event", // 
                                                                    "InteractionTarget",
                                                                    "InteractionVolume_short",
                                                                    "IncomingNeutrino_energy",
                                                                    "Neutrino_reconstructed_energy_GeV",
                                                                    });
                        
    return 0;
}


    // // NEUTRINOS
    // TCanvas* c_neutrino_E = new TCanvas("neutrino_true_E", "interacting neutrino true energy", 800, 600);
    
    // interacting_neutrino_true_E->SetLineColor(kRed);     
    // interacting_neutrino_true_E_cut_wires->SetLineColor(kBlue);    
    // interacting_neutrino_true_E_charge_multi->SetLineColor(kGreen);
    // interacting_neutrino_true_E_selection->SetLineColor(kMagenta); 
    
    // interacting_neutrino_true_E->Scale(1./interacting_neutrino_true_E->Integral());     
    // interacting_neutrino_true_E_cut_wires->Scale(1./interacting_neutrino_true_E_cut_wires->Integral());   
    // interacting_neutrino_true_E_charge_multi->Scale(1./interacting_neutrino_true_E_charge_multi->Integral()); 
    // interacting_neutrino_true_E_selection->Scale(1./interacting_neutrino_true_E_selection->Integral());

    // interacting_neutrino_true_E->Draw("");
    // interacting_neutrino_true_E_cut_wires->Draw("SAME");
    // interacting_neutrino_true_E_charge_multi->Draw("SAME");
    // interacting_neutrino_true_E_selection->Draw("SAME");

    // TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Posizione della legenda
    // legend->AddEntry(interacting_neutrino_true_E.GetPtr(), "[STAGE 1: FIDUCIAL VOLUME]", "l");
    // legend->AddEntry(interacting_neutrino_true_E_cut_wires.GetPtr(), "[STAGE 2: CUT NOF WIRES]", "l");
    // legend->AddEntry(interacting_neutrino_true_E_charge_multi.GetPtr(), "[STAGE 3: CHARGE MULTIPLICITY]", "l");
    // legend->AddEntry(interacting_neutrino_true_E_selection.GetPtr(), "[STAGE 4: SIGNAL SELECTION]", "l");
    // legend->Draw();
    // c_neutrino_E->Update();
    // report_file->cd();
    // c_neutrino_E->Write();
    // //

    // // NEUTRINOS
    // TCanvas* c_antimuon_E = new TCanvas("antimuon", "antimuon true energy", 800, 600);
    
    // antimuon_true_E->SetLineColor(kRed);   
    // antimuon_true_E_cut_wires->SetLineColor(kBlue);    
    // antimuon_true_E_charge_multi->SetLineColor(kGreen);   
    // antimuon_true_E_selection->SetLineColor(kMagenta); 
    
    // antimuon_true_E->Scale(1./antimuon_true_E->Integral());   
    // antimuon_true_E_cut_wires->Scale(1./antimuon_true_E_cut_wires->Integral());    
    // antimuon_true_E_charge_multi->Scale(1./antimuon_true_E_charge_multi->Integral());   
    // antimuon_true_E_selection->Scale(1./antimuon_true_E_selection->Integral()); 

    // antimuon_true_E->Draw("");
    // antimuon_true_E_cut_wires->Draw("SAME");
    // antimuon_true_E_charge_multi->Draw("SAME");
    // antimuon_true_E_selection->Draw("SAME");

    // TLegend* legend2 = new TLegend(0.7, 0.7, 0.9, 0.9); 
    // legend2->AddEntry(antimuon_true_E.GetPtr(), "[STAGE 1: FIDUCIAL VOLUME]", "l");
    // legend2->AddEntry(antimuon_true_E_cut_wires.GetPtr(), "[STAGE 2: CUT NOF WIRES]", "l");
    // legend2->AddEntry(antimuon_true_E_charge_multi.GetPtr(), "[STAGE 3: CHARGE MULTIPLICITY]", "l");
    // legend2->AddEntry(antimuon_true_E_selection.GetPtr(), "[STAGE 4: SIGNAL SELECTION]", "l");
    // legend2->Draw();
    // c_antimuon_E->Update();
    // report_file->cd();
    // c_antimuon_E->Write();
    // //

    // // HADRONIC SYSTEM
    // TCanvas* c_hadronic_E = new TCanvas("hadronic_system", "final state hadronic system true Kin E", 800, 600);

    // fs_hadron_syst_Kin->SetLineColor(kRed);   
    // fs_hadron_syst_Kin_cut_wires->SetLineColor(kBlue);    
    // fs_hadron_syst_Kin_charge_multi->SetLineColor(kGreen);   
    // fs_hadron_syst_Kin_selection->SetLineColor(kMagenta); 

    // fs_hadron_syst_Kin->Scale(1.0 / fs_hadron_syst_Kin->Integral());
    // fs_hadron_syst_Kin_cut_wires->Scale(1.0 / fs_hadron_syst_Kin_cut_wires->Integral());
    // fs_hadron_syst_Kin_charge_multi->Scale(1.0 / fs_hadron_syst_Kin_charge_multi->Integral());
    // fs_hadron_syst_Kin_selection->Scale(1.0 / fs_hadron_syst_Kin_selection->Integral());

    // c_hadronic_E->SetLogy();

    // fs_hadron_syst_Kin->Draw("");
    // fs_hadron_syst_Kin_cut_wires->Draw("SAME");
    // fs_hadron_syst_Kin_charge_multi->Draw("SAME");
    // fs_hadron_syst_Kin_selection->Draw("SAME");

    // // Creazione della legenda
    // TLegend* legend3 = new TLegend(0.7, 0.7, 0.9, 0.9); 
    // legend3->AddEntry(fs_hadron_syst_Kin.GetPtr(), "[STAGE 1: FIDUCIAL VOLUME]", "l");
    // legend3->AddEntry(fs_hadron_syst_Kin_cut_wires.GetPtr(), "[STAGE 2: CUT NOF WIRES]", "l");
    // legend3->AddEntry(fs_hadron_syst_Kin_charge_multi.GetPtr(), "[STAGE 3: CHARGE MULTIPLICITY]", "l");
    // legend3->AddEntry(fs_hadron_syst_Kin_selection.GetPtr(), "[STAGE 4: SIGNAL SELECTION]", "l");
    // legend3->Draw();

    // // Aggiorna il canvas e scrivi nel file
    // c_hadronic_E->Update();
    // report_file->cd();
    // c_hadronic_E->Write();
    //