#ifndef RDF_UTILS_H
#define RDF_UTILS_H

#include <string>

#include "GenieUtils.h"
#include "RecoUtils.h"
#include "struct.h" // struct sandreco
#include "evtinfo.h" // struct fastreco

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TString.h"
#include "TObject.h"
#include "TObjString.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TVector3.h"

namespace RDFUtils{//RDFUtils

ROOT::RDataFrame InitDF(const char* production, const char* tree_name);

void PrintColumns(ROOT::RDataFrame& df);

ROOT::RDF::RNode AddConstantsToDF(ROOT::RDataFrame& df);

template<int coord>
ROOT::VecOps::RVec<double> GetComponent(const ROOT::VecOps::RVec<TLorentzVector>& vTL);

double GetColumnSum(const ROOT::VecOps::RVec<double>& v);

ROOT::VecOps::RVec<double> VectorDifference(const ROOT::VecOps::RVec<double>& v1,
                                            const ROOT::VecOps::RVec<double>& v2);

TVector3 TLVectorCrossProduct(const TLorentzVector& v1,
                              const TLorentzVector& v2);

ROOT::VecOps::RVec<double> VectorSubtractConst(const ROOT::VecOps::RVec<double>& v1, double c);

TLorentzVector VectorFilterByHighest(const ROOT::VecOps::RVec<double>& filter,
                                     const ROOT::VecOps::RVec<TLorentzVector>& v);

// ROOT::VecOps::RVec<double> GetResolution(const ROOT::VecOps::RVec<double>& reco,
//                                          const ROOT::VecOps::RVec<double>& true);

namespace GENIE{//GENIE

std::string InteractionTarget(const ROOT::VecOps::RVec<int>& pdg);

std::string EventType(TObjString& t);

ROOT::VecOps::RVec<genie::GHepParticle> AllGenieParticles(const ROOT::VecOps::RVec<int>& pdg,
                                                          const ROOT::VecOps::RVec<int>& status,
                                                          const ROOT::VecOps::RVec<double>& P4);                                                

template<genie::GHepStatus_t STATUS>
ROOT::VecOps::RVec<genie::GHepParticle> GetParticlesWithStatus(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

template<int PDG>
ROOT::VecOps::RVec<genie::GHepParticle> GetParticlesWithPDG(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<genie::GHepParticle> GetExhoticMesons(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<genie::GHepParticle> GetExhoticHadrons(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<genie::GHepParticle> GetRecoiledNuclei(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<genie::GHepParticle> GetFinalHadronicSystem(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

GenieUtils::event_topology GetFinalStateTopology(const ROOT::VecOps::RVec<int>& pdgs);

ROOT::VecOps::RVec<int> GetPDG(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<TLorentzVector> GetMomentum(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

TLorentzVector SumLorentzVectors(const ROOT::VecOps::RVec<TLorentzVector>& VTL);

int NofFinalStateParticles(const ROOT::VecOps::RVec<int>& pdg);

double GetInitialNucleonMomentum(const TLorentzVector& FinalStateMomentum,
                                 const TVector3& FinalStateLongitudinalMomentum,
                                 const TVector3& FinalStateDeltaPT,
                                 int InteractionTarget);

ROOT::RDF::RNode AddColumnsFromGENIE(ROOT::RDF::RNode& df);

ROOT::RDF::RNode AddColumnsForHydrogenCarbonSampleSelection(ROOT::RDF::RNode& df); // these is expected to have columns from AddColumnsFromGENIE output


}//GENIE

namespace SANDRECO{//SANDRECO

ROOT::VecOps::RVec<particle> ParticleGoodTrackFit(const ROOT::VecOps::RVec<particle>& particles);

ROOT::VecOps::RVec<particle> FilterPrimaries(const ROOT::VecOps::RVec<particle>& particles);

template<int PDG>
ROOT::VecOps::RVec<particle> GetParticlesWithPDG(const ROOT::VecOps::RVec<particle>& particles);

ROOT::VecOps::RVec<TLorentzVector> GetMomentum(const ROOT::VecOps::RVec<particle>& particles);

ROOT::VecOps::RVec<TLorentzVector> GetTrackVertex(const ROOT::VecOps::RVec<particle>& particles);

ROOT::RDF::RNode AddColumnsFromSANDRECO(ROOT::RDF::RNode& df);

}//SANDRECO

namespace FASTRECO{//FASTRECO

template<int PDG>
ROOT::VecOps::RVec<fast::particle> GetParticlesWithPDG(const std::map<int, fast::particle>& particle_map);

ROOT::VecOps::RVec<TVector3> GetTrackVertex(const ROOT::VecOps::RVec<fast::particle>& particles);

template<int coord>
double GetEventVertex(const ROOT::VecOps::RVec<TVector3>& vertices);

ROOT::VecOps::RVec<TLorentzVector> GetTrue4Momentum(const ROOT::VecOps::RVec<fast::particle>& particles);

ROOT::VecOps::RVec<TLorentzVector> GetReco4Momentum(const ROOT::VecOps::RVec<fast::particle>& particles);

ROOT::RDF::RNode AddColumnsFromFASTRECO(ROOT::RDF::RNode& df);

}//FASTRECO

namespace GEO{//GEO

std::string GetMaterialFromCoordinates(double x, double y, double z);

}

}//RDFUtils

#endif