#include <iostream>
#include <string>

#include "RDataFrameUtils.h"
#include "GenieUtils.h"
#include "GeoUtils.h"
#include "TFile.h"

TGeoManager* geo = nullptr;

//___________________________________________________________________
int main(int argc, char* argv[]){
    // Analyze output of genie production

    if(argc<4)
    {
        std::cout<<"AnalyseGenie <GENIE PRODUCTION> <GEOMETRY> <FILE OUTPUT NAME>\n";
        throw "";
    }

    const char* fInput = argv[1];

    const char* geometry = argv[2];

    const char* fOutput = argv[3];

    if(!std::strstr(fInput,"*")) ROOT::EnableImplicitMT();

    auto df = RDFUtils::InitDF(fInput, "gRooTracker");

    // RDFUtils::PrintColumns(df);

    auto dfG = RDFUtils::GENIE::AddColumnsFromGENIE(df);

    geo = TGeoManager::Import(geometry);

    dfG.Snapshot("myTree",fOutput,{"InteractionTarget",
                                   "InteractionTargetFromGEO",
                                   "EventType",
                                   // initial state particles
                                   "InitialStateParticlesPDG",
                                   "InitialStateParticlesE",
                                   "InitialStateEnergy",
                                   // stable final state particles
                                   "StableFinalStateParticlesPDG",
                                   "StableFinalStateParticlesE",
                                   "StableFinalStateEnergy",
                                   // final state nuclear remnant
                                   "FinalStateNuclearRemnantPDG",
                                   "FinalStateNuclearRemnantE",
                                   "FinalStateNuclearRemnantEnergy"
                                   });

    return 0;
}
//___________________________________________________________________