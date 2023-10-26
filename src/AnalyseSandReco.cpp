#include <iostream>
#include <string>

#include "RDataFrameUtils.h"
#include "GenieUtils.h"
#include "GeoUtils.h"
#include "TFile.h"

int main(int argc, char* argv[]){
    // Analyze output of sandreco 

    if(argc<4)
    {
        std::cout<<"AnalyseGenie <SAND RECO> <GEOMETRY> <FILE OUTPUT NAME>\n";
        throw "";
    }

    // read user inputs

    const char* fInput = argv[1];

    const char* geometry = argv[2];

    const char* fOutput = argv[3];

    // if you have multiple files enable multiple thread pocessing

    if(!std::strstr(fInput,"*")) ROOT::EnableImplicitMT();

    // Initialize root DataFrame and add columns

    auto df = RDFUtils::InitDF(fInput, "tEvent");

    // RDFUtils::PrintColumns(df);

    
}