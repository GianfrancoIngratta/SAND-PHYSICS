#ifndef RDF_UTILS_H
#define RDF_UTILS_H

#include <string>

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

namespace RDFUtils{//RDFUtils

ROOT::RDataFrame InitDF(const char* production, const char* tree_name);

void PrintColumns(ROOT::RDataFrame& df);

namespace GENIE{//GENIE

std::string InteractionTarget(const ROOT::VecOps::RVec<int>& pdg);

ROOT::RDF::RNode AddColumnsFromGENIE(ROOT::RDataFrame& df);

}//GENIE

}//RDFUtils

#endif