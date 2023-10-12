#include <iostream>

#include "RDataFrameUtils.h"
#include "GenieUtils.h"
#include "TFile.h"


//___________________________________________________________________
int main(int argc, char* argv[]){
    // Analyze output of 1 genie production

    if(argc!=2)
    {
        std::cout<<"AnalyseGenie <GENIE PRODUCTION>\n";
        throw "";
    }

    const char* fInput = argv[1];

    auto df = RDFUtils::InitDF(fInput, "gRooTracker");

    // RDFUtils::PrintColumns(df);

    auto dfG = RDFUtils::GENIE::AddColumnsFromGENIE(df);

    dfG.Snapshot("myTree","prova.root","InteractionTarget");

    return 0;
}
//___________________________________________________________________