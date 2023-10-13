#ifndef RDF_UTILS_H
#define RDF_UTILS_H

#include <string>

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TString.h"
#include "TObject.h"
#include "TObjString.h"
#include "TLorentzVector.h"

namespace RDFUtils{//RDFUtils

ROOT::RDataFrame InitDF(const char* production, const char* tree_name);

void PrintColumns(ROOT::RDataFrame& df);

template<int coord>
ROOT::VecOps::RVec<double> GetComponent(const ROOT::VecOps::RVec<TLorentzVector>& vTL);

namespace GENIE{//GENIE

std::string InteractionTarget(const ROOT::VecOps::RVec<int>& pdg);

ROOT::VecOps::RVec<int> StableFinalStateParticles(const ROOT::VecOps::RVec<int>& pdg,
                                                  const ROOT::VecOps::RVec<int>& status);

int NofFinalStateParticles(const ROOT::VecOps::RVec<int>& pdg);

std::string EventType(TObjString& t);

template<int PDG>
ROOT::VecOps::RVec<TLorentzVector> GetMomentum(const ROOT::VecOps::RVec<int>& pdg,
                                               const ROOT::VecOps::RVec<int>& status,
                                               const ROOT::VecOps::RVec<double>& P4);

ROOT::RDF::RNode AddColumnsFromGENIE(ROOT::RDataFrame& df);

}//GENIE

namespace GEO{//GEO

std::string GetMaterialFromCoordinates(double x, double y, double z);

}

}//RDFUtils

#endif