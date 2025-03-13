#include <iostream>
#include <string>

#include "RDataFrameUtils.h"
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
    unsigned int files_per_jobs = 5000u;
    unsigned int file_start = index * files_per_jobs;
    unsigned int file_stop = index * files_per_jobs + files_per_jobs;

    LOG("I", TString::Format("Analyze production from file %d to file %d", file_start, file_stop).Data());
    
    LOG("I", "Reading geometry");
    geo = TGeoManager::Import("/storage/gpfs_data/neutrino/users/gi/dunendggd/SAND_opt3_DRIFT1.root");

    auto FOLDER_PRODUCTION = "/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc/";
    auto FOLDER_ANALYSIS = "/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/antinumu_CCQE_on_H_like/genie/";

    auto fInput_genie = TString::Format("%sevents-in-SANDtracker.*.gtrac.root", FOLDER_PRODUCTION);
    auto fOutput = TString::Format("%sevents-in-SANDtracker.from.%d.to.%d.genie_extended.root", FOLDER_ANALYSIS, file_start, file_stop);

    if(fInput_genie.Contains("*")){
        LOG("I","Enabling multiple threading");
        ROOT::EnableImplicitMT();
        geo->SetMaxThreads(100);
    };

    LOG("I", "Initialize ROOT DataFrame");
    auto chain_genie = RDFUtils::InitTChain(fInput_genie, "gRooTracker", file_start, file_stop); 
    auto df = RDFUtils::InitDF(chain_genie);

    auto dfC = RDFUtils::AddConstantsToDF(df); // add some columns with usefull constants

    auto dfG = RDFUtils::GENIE::AddColumnsFromGENIE(dfC);
    LOG("I", "Filtere Fiducial Volume");
    auto df_filtered = RDFUtils::Filter(dfG, "isInFiducialVolume==1", true);
    // auto dfG_SolidHydrogen = RDFUtils::GENIE::AddColumnsForHydrogenCarbonSampleSelection(dfG); // add columns for Hydrogen sample selection
    LOG("I", "Writing ouput file");
    dfG.Snapshot("gtrac_extended",fOutput.Data(), {
                                            "CCQEonHydrogen",
                                            "EvtXSec",
                                            "EvtDXSec",
                                            "EvtProb",
                                            "InteractionVolume",
                                            "InteractionVolume_short",
                                            "isInFiducialVolume",
                                            "NeutrinoFlavor",
                                            "isCCEvent",
                                            "EventType",
                                            "Interaction_vtxX",
                                            "Interaction_vtxY",
                                            "Interaction_vtxZ",
                                            "Interaction_vtxT",
                                            "InteractionTarget",
                                            "IncomingNeutrinoP4",
                                            "NucleonTargetP4",
                                            "FinalStateLepton4Momentum",
                                            "FinalStateLeptonPDG",
                                            "FinalStateHadronicSystemTopology_name",
                                            "FinalStateHadronicSystemTotal4Momentum",
                                            });

    return 0;
}
//___________________________________________________________________