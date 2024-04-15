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

    if(argc < 4 || argc > 7)
    {
        LOG("W", "AnalyseEDepSim -i <EDEP PRODUCTION> -g <GEOMETRY> -o <FILE OUTPUT>\n");
        throw "";
    }

    // read user inputs

    const char* fInput = argv[1];

    const char* geometry = argv[2];

    const char* fOutput = argv[3];

    unsigned int start = 0;

    unsigned int stop = 1;

    int index = 1;

    LOG("I", "Parsing inputs");

    while (index < argc)
    {
        TString opt = argv[index];
        if(opt.CompareTo("-i")==0){
            try
            {
                fInput = argv[++index];
                LOG("ii",TString::Format("Input file : %s", fInput).Data());
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << '\n';
                return 1;
            }
        }else if(opt.CompareTo("-g")==0){
            try
            {
                geometry = argv[++index];
                LOG("ii",TString::Format("Geometry file : %s", geometry).Data());
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << '\n';
                return 1;
            }
        }else if(opt.CompareTo("-o")==0){
            try
            {
                fOutput = argv[++index];
                LOG("ii",TString::Format("Output file : %s", fOutput).Data());
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << '\n';
            }
        }        else{
            auto ui = argv[++index];
            LOG("W", TString::Format("Unknown Input : %s", ui).Data());
            return 1; 
        }
        index++;
    }

    // if you have multiple files enable multiple thread pocessing

    if(TString::Format("%s",fInput).Contains("*")){
        LOG("I","Enabling multiple threading");
        ROOT::EnableImplicitMT();
        stop = start + 10u;
        };
    
    LOG("I", "Reading geometry");
    geo = TGeoManager::Import(geometry);

    LOG("I", "Initialize ROOT DataFrame");
    auto df = RDFUtils::InitDF(fInput, "EDepSimEvents", start, stop);

    auto dfC = RDFUtils::AddConstantsToDF(df); // add some columns with usefull constants

    auto dfEDEP = RDFUtils::EDEPSIM::AddColumnsFromEDEPSIM(dfC);
    
    LOG("I", "Writing ouput file");
    dfEDEP.Snapshot("edep_extended", fOutput, {
                                                "NofEventsPerSpill",
                                                "PrimariesVertexX",
                                                "PrimariesVertexY",
                                                "PrimariesVertexZ",
                                                "PrimariesVertexT",
                                                "PrimariesVertexTimeDiff",

    });                                                    

    return 0;
}