#ifndef RDF_UTILS_H
#define RDF_UTILS_H

#include <string>

#include "GenieUtils.h"
#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TString.h"
#include "TObject.h"
#include "TObjString.h"
#include "TLorentzVector.h"
#include "TVector2.h"

namespace RDFUtils{//RDFUtils

ROOT::RDataFrame InitDF(const char* production, const char* tree_name);

void PrintColumns(ROOT::RDataFrame& df);

template<int coord>
ROOT::VecOps::RVec<double> GetComponent(const ROOT::VecOps::RVec<TLorentzVector>& vTL);

double GetColumnSum(const ROOT::VecOps::RVec<double>& v);

ROOT::VecOps::RVec<double> VectorDifference(const ROOT::VecOps::RVec<double>& v1,
                                            const ROOT::VecOps::RVec<double>& v2);

ROOT::VecOps::RVec<double> VectorSubtractConst(const ROOT::VecOps::RVec<double>& v1, double c);

TLorentzVector VectorFilterByHighest(const ROOT::VecOps::RVec<double>& filter,
                                     const ROOT::VecOps::RVec<TLorentzVector>& v);                                            

namespace GENIE{//GENIE

std::string InteractionTarget(const ROOT::VecOps::RVec<int>& pdg);

std::string EventType(TObjString& t);

ROOT::VecOps::RVec<genie::GHepParticle> AllGenieParticles(const ROOT::VecOps::RVec<int>& pdg,
                                                          const ROOT::VecOps::RVec<int>& status,
                                                          const ROOT::VecOps::RVec<double>& P4);                                                

template<genie::GHepStatus_t STATUS>
ROOT::VecOps::RVec<genie::GHepParticle> GetParticlesWithStatus(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

template<int component>
ROOT::VecOps::RVec<double> GetMomentum(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<int> GetPDG(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);  

int NofFinalStateParticles(const ROOT::VecOps::RVec<int>& pdg);

ROOT::RDF::RNode AddColumnsFromGENIE(ROOT::RDataFrame& df);

}//GENIE

namespace GEO{//GEO

std::string GetMaterialFromCoordinates(double x, double y, double z);

}

}//RDFUtils

#endif