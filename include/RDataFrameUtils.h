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

ROOT::VecOps::RVec<int> StableFinalStateParticles(const ROOT::VecOps::RVec<int>& pdg,
                                           const ROOT::VecOps::RVec<int>& status);

int NofFinalStateParticles(const ROOT::VecOps::RVec<int>& pdg);

ROOT::RDF::RNode AddColumnsFromGENIE(ROOT::RDataFrame& df);

}//GENIE

namespace GEO{//GEO

std::string GetMaterialFromCoordinates(double x, double y, double z);

}

}//RDFUtils

#endif