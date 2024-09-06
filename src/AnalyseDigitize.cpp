#include <iostream>
#include <string>

#include "RDataFrameUtils.h"
// #include "GenieUtils.h"
// #include "GeoUtils.h"
#include "TFile.h"

TGeoManager* geo = nullptr;

//___________________________________________________________________
int main(int argc, char* argv[]){

    unsigned int start = 2;
    
    LOG("I", "Reading geometry");
    // geo = TGeoManager::Import("/storage/gpfs_data/neutrino/users/gi/dunendggd/SAND_opt3_DRIFT1.root");
    // geo = TGeoManager::Import("/storage/gpfs_data/neutrino/users/gi/dunendggd/SAND_opt3_DRIFT1.gdml");
    TFile f("/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc/events-in-SANDtracker.1.edep-sim.root", "READ");
    geo = (TGeoManager*)f.Get("EDepSimGeometry");

    auto fInput_genie = "/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc/events-in-SANDtracker.*.gtrac.root";
    auto fInput_edep = "/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc/events-in-SANDtracker.*.edep-sim.root";
    auto fInput_digit = "/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc/events-in-SANDtracker.*.ecal-digit.root";
    auto fInput_ecal_cluster = "/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc/events-in-SANDtracker.*.ecal-cluster.root";

    auto fOutput = "/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/events-in-SANDtracker.0.ecal-digit.analysed.root";
    
    // if you have multiple files enable multiple thread pocessing
    if(TString::Format("%s",fInput_digit).Contains("*")){
        LOG("I","Enabling multiple threading");
        ROOT::EnableImplicitMT();
        geo->SetMaxThreads(100);
    };

    LOG("I", "Initialize ROOT DataFrame");
    auto chain_genie = RDFUtils::InitTChain(fInput_genie, "gRooTracker", start, start + 999u);
    auto chain_edep = RDFUtils::InitTChain(fInput_edep, "EDepSimEvents", start, start + 999u);
    auto chain_digit = RDFUtils::InitTChain(fInput_digit, "tDigit", start, start + 999u); 
    // auto chain_cluster = RDFUtils::InitTChain(fInput_ecal_cluster, "tReco", start, start + 999u); 
    
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
    
    LOG("I", "Writing ouput file");
    dfDigit.Snapshot("digit_extended", fOutput, {
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
                                                /*
                                                    GENIE INFO
                                                */
                                                "FinalStateLeptonEmissionAngle",
                                                "PrimaryStateHadronicSystemTotalKinE",
                                                /*
                                                    EDEP INFO
                                                */
                                                "PrimariesPDG",
    //                                             "PrimariesTrackId",
                                                "PrimariesP4",
                                                "PrimariesFirstHitECAL",
                                                "PrimariesEDepECAL",
                                                "PrimariesEmissionAngle",
                                                /*
                                                    PREDICTIONS FOR CHANNEL antinu on H
                                                */
                                                "ExpectedNeutrinoP4FromMuon",
                                                "ExpectedHadronSystP3",
                                                "ExpectedHadronSystEnergy",
                                                "ExpectedNeutronArrivalPositionECAL",
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
                                                 "isCellComplete",
                                                 "Cell_Reconstructed_hits",
                                                 "STDistToNeutronExpectedHit",
                                                 /*
                                                    ECAL CLUSTER INFO
                                                 */
                                                //  "NofEventClusters",
                                                //  "ClusterX4",
                                                //  "Cluster2Vertex4Distance",
    });                                                    

    return 0;
}