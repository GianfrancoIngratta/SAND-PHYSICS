#include <iostream>
#include <string>

#include "RDataFrameUtils.h"
// #include "GenieUtils.h"
// #include "GeoUtils.h"
#include "TFile.h"

TGeoManager* geo = nullptr;

struct SIGNAL
{ // antinumu + H -> n + mu+
    static const int charge_multuplicity = 1;
};


// FUNCTIONS

ROOT::RDF::RNode Filter_ChargeMultiplicity(ROOT::RDF::RNode& df){
    /*
        CUT : Require a number of charget particles
              associated with the vertex to be 1.
              This cut is temporary done on MC truth
    */
    auto condition_multiplicity = TString::Format("NofFinalStateChargedParticles==%d", SIGNAL::charge_multuplicity);
    std::cout << "Filtering only events with 1 final state charged particle\n";
    return df.Filter(condition_multiplicity.Data());
}

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
    auto fOutput = TString::Format("%sevents-in-SANDtracker.%d.to.%d.event_selection.root",FOLDER_ANALYSIS, file_start, file_stop);
    
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

    auto df_cut1 = Filter_ChargeMultiplicity(dfDigit);
    
    LOG("I", "Writing ouput file");
    df_cut1.Snapshot("digit_extended", fOutput.Data(), {
                                                /*
                                                    GENIE INFO
                                                */
                                                "FileName",
                                                "EventId",
                                                "EventType",
                                                "CCQEonHydrogen",
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
                                                "PrimariesFirstHitECAL",
                                                "PrimariesEDepECAL",
                                                "PrimariesEmissionAngle",
                                                "PrimaryHasNoECALHit",
                                                "PrimaryHasChangedDirection",
                                                /*
                                                    PREDICTIONS FOR CHANNEL antinu on H
                                                */
                                                "ExpectedNeutrinoP4FromMuon",
                                                "ExpectedHadronSystP3",
                                                "ExpectedHadronSystEnergy",
                                                "ExpectedNeutronArrivalPositionECAL",
                                                "ExpectedNeutronTOF",
                                                "ExpectedFiredModuleByNeutron",
                                                "MissingTransverseMomentum",
    //                                             "DoubleTransverseMomentumImbalance",
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
                                                "Fired_Cells_tdc1",
                                                "Fired_Cells_tdc2",
                                                "who_produced_tdc1",
                                                "who_produced_tdc2",
                                                "isCellComplete",
                                                "Cell_Reconstructed_hits",
                                                "ExpectedNeutronHit",
                                                 /*
                                                    ECAL CLUSTER INFO
                                                 */
                                                //  "NofEventClusters",
                                                //  "ClusterX4",
                                                //  "Cluster2Vertex4Distance",
    });                                                    

    return 0;
}