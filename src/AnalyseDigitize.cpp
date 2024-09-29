#include <iostream>
#include <string>

#include "RDataFrameUtils.h"
// #include "GenieUtils.h"
// #include "GeoUtils.h"
#include "TFile.h"

TGeoManager* geo = nullptr;

//___________________________________________________________________
int main(int argc, char* argv[]){

    // analyse file from file file_start" to "file_start + nof files 
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
    auto FOLDER_ANALYSIS = "/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/";

    auto fInput_genie =TString::Format("%sevents-in-SANDtracker.*.gtrac.root", FOLDER_PRODUCTION);
    auto fInput_edep =TString::Format("%sevents-in-SANDtracker.*.edep-sim.root", FOLDER_PRODUCTION);
    auto fInput_digit =TString::Format("%sevents-in-SANDtracker.*.ecal-digit.root", FOLDER_PRODUCTION);
    // auto fInput_ecal_cluster = "events-in-SANDtracker.*.ecal-cluster.root";
    auto fOutput = TString::Format("%sevents-in-SANDtracker.%d.to.%d.ecal-digit.analysed.root",FOLDER_ANALYSIS, file_start, file_stop);

    
    // if you have multiple files enable multiple thread pocessing
    // if(TString::Format("%s",fInput_digit).Contains("*")){
    if(fInput_digit.Contains("*")){
        LOG("I","Enabling multiple threading");
        ROOT::EnableImplicitMT();
        geo->SetMaxThreads(100);
    };

    LOG("I", "Initialize ROOT DataFrame");
    auto chain_genie = RDFUtils::InitTChain(fInput_genie, "gRooTracker", file_start, file_stop);
    auto chain_edep = RDFUtils::InitTChain(fInput_edep, "EDepSimEvents", file_start, file_stop);
    auto chain_digit = RDFUtils::InitTChain(fInput_digit, "tDigit", file_start, file_stop); 
    // auto chain_cluster = RDFUtils::InitTChain(fInput_ecal_cluster, "tReco", file_start, file_stop); 
    
    chain_digit->AddFriend(chain_genie, "genie");
    chain_digit->AddFriend(chain_edep, "edep");
    // chain_digit->AddFriend(chain_cluster, "cluster");
    
    auto df = RDFUtils::InitDF(chain_digit);

    // RDFUtils::PrintColumns(df);

    // throw "";

    auto dfC = RDFUtils::AddConstantsToDF(df); // add some columns with usefull constants

    auto dfGENIE = RDFUtils::GENIE::AddColumnsFromGENIE(dfC);

    auto dfEDEP = RDFUtils::EDEPSIM::NOSPILL::AddColumnsFromEDEPSIM(dfGENIE);
    
    auto dfDigit = RDFUtils::DIGIT::AddColumnsFromDigit(dfEDEP);

    /*
        Filter
        - Fiducial Volume
        - Signal events (MC truth)
        - Complete Cells (reco)
    */

    LOG("I", "Filter events in fiducial volume (MC truth)");
    dfDigit = dfDigit.Filter("isInFiducialVolume"); // genie

    LOG("I", "Filter signal events (MC truth)");
    dfDigit = dfDigit.Filter("CCQEonHydrogen==1"); // genie

    auto filter_description = "FV_Signal";
    // auto filter_description = "FV";

    auto fOutput_filtered = TString::Format("%sevents-in-SANDtracker.%d.to.%d.ecal-digit.analysed.%s.root",FOLDER_ANALYSIS, file_start, file_stop, filter_description);
    auto fOutput_filtered_trj = TString::Format("%sevents-in-SANDtracker.%d.to.%d.ecal-digit.analysed.%s_trj.root",FOLDER_ANALYSIS, file_start, file_stop, filter_description);
    
    LOG("I", "Writing ouput file");
    dfDigit.Snapshot("digit_extended", fOutput_filtered.Data(), {
                                                /*
                                                    GENIE INFO
                                                */
                                                "FileName",
                                                "EventId",
                                                "EventType",
                                                "CCQEonHydrogen",
                                                "NuDirection",
                                                "NofEvents",
                                                "Interaction_vtxX",
                                                "Interaction_vtxY",
                                                "Interaction_vtxZ",
                                                "Interaction_vtxT",
                                                "InteractionVolume",
                                                "NofFinalStateChargedParticles",
                                                "FinalStateLeptonEmissionAngle",
                                                "PrimaryStateHadronicSystemTotalKinE",
                                                "PrimaryStateHadronicSystemTopology_name",
                                                "InteractionTarget",
                                                /*
                                                    EDEP INFO
                                                */
                                                "PrimariesPDG",
                                                "PrimariesTrackId",
                                                "PrimariesP4",
                                                "PrimariesBeta",
                                                "PrimariesFirstHitECAL",
                                                "PrimariesEDepECAL",
                                                "PrimariesEmissionAngle",
                                                "IsECALHitMissing",
                                                "DeviationAngle",
                                                /*
                                                    PREDICTIONS FOR CHANNEL antinu on H
                                                */
                                                "ExpectedNeutrinoP4FromMuon",
                                                "ExpectedHadronSystP3",
                                                "ExpectedHadronSystEnergy",
                                                "MissingTransverseMomentum",
                                                 /*
                                                      DIGIT INFO
                                                 */
                                                "NofEventFiredModules",
                                                "EventFiredModules",
                                                "Fired_Cells_mod",
                                                "Fired_Cells_id",
                                                "Fired_Cells_x",
                                                "Fired_Cells_y",
                                                "Fired_Cells_z",
                                                "Fired_Cells_adc1",
                                                "Fired_Cells_adc2",
                                                "Fired_Cells_tdc1",
                                                "Fired_Cells_tdc2",
                                                "Fired_Cell_true_hit1",
                                                "Fired_Cell_true_hit2",
                                                "Fired_by_primary_neutron",
                                                "Fired_by_primary_antimu",
                                                // "who_produced_tdc2",
                                                "isCellComplete",
                                                /*
                                                    true
                                                */
                                                "Fired_Cell_true_Hit_x",
                                                "Fired_Cell_true_Hit_y",
                                                "Fired_Cell_true_Hit_z",
                                                "Fired_Cell_true_Hit_t",
                                                "True_FlightLength",
                                                /*
                                                    expected
                                                */
                                                "ExpectedNeutron_Beta",
                                                "ExpectedNeutron_HitPosition_x",
                                                "ExpectedNeutron_HitPosition_y",
                                                "ExpectedNeutron_HitPosition_z",
                                                "ExpectedNeutron_TOF",
                                                "ExpectedNeutron_FlightLength",
                                                /*
                                                    reco
                                                */
                                                "Reconstructed_HitPosition_x",
                                                "Reconstructed_HitPosition_y",
                                                "Reconstructed_HitPosition_z",
                                                "Reconstructed_HitTime",
                                                "Reconstructed_FlightLength",
                                                //
                                                "Residuals_HitTime",
                                                "Residuals_HitSpace",
                                                "IsSpaceCompatible",
                                                // "isCandidate",
    });                                                    

    LOG("I", "Writing trajectory file");
    auto dfTraj = RDFUtils::DIGIT::GetFilteredTrajectories(dfDigit);
    
    dfTraj.Snapshot("trj_extended", fOutput_filtered_trj.Data(), {
                                                 "FileName",
                                                 "trackid",
                                                 "pdg",
                                                 "point_x",
                                                 "point_y",
                                                 "point_z",
                                                 "point_px",
                                                 "point_py",
                                                 "point_pz",
                                                 "process",

    });
    
    return 0;
}