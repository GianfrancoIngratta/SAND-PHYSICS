#ifndef RDF_UTILS_H
#define RDF_UTILS_H

#include <string>

#include "GenieUtils.h"
#include "TSystem.h"
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
#include "TG4Event.h"
#include "TG4HitSegment.h"

namespace Color {
    enum Code {
        FG_RED      = 31,
        FG_GREEN    = 32,
        FG_BLUE     = 34,
        FG_DEFAULT  = 39,
        BG_RED      = 41,
        BG_GREEN    = 42,
        BG_BLUE     = 44,
        BG_DEFAULT  = 49
    };
    class Modifier {
        Code code;
    public:
        Modifier(Code pCode) : code(pCode) {}
        friend std::ostream&
        operator<<(std::ostream& os, const Modifier& mod) {
            return os << "\033[" << mod.code << "m";
        }
    };
}

Color::Modifier def(Color::FG_DEFAULT);
Color::Modifier red(Color::FG_RED); // warnings
Color::Modifier green(Color::FG_GREEN); // opeartions
Color::Modifier blue(Color::FG_BLUE); // results

void LOG(TString i, const char* out){
    if(i.Contains("I")){ // info
        std::cout << green << "[INFO] " << out << def << std::endl;
        std::cout << "\n";
    }else if(i.Contains("i")){
        std::cout << green << " |___ " << out << def << std::endl;
        std::cout << "\n";
    }else if(i.Contains("R")){ // result
        std::cout << blue << "[RESULT] " << out << def << std::endl;
        std::cout << "\n";
    }else if(i.Contains("W")){ // waring
        std::cout << red << "[WARNING] " << out << def << std::endl;
        std::cout << "\n";
    }else{
        std::cout << out;
    }
}

namespace RDFUtils{//RDFUtils

ROOT::RDataFrame InitDF(TString production, const char* tree_name, unsigned int file_index_start = 0, unsigned int file_index_stop = 10);

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

std::string GetVolumeFromCoordinates(double x, double y, double z);

}//GEO

namespace EDEPSIM{//EDEPSIM

ROOT::RDF::RNode AddColumnsFromEDEPSIM(ROOT::RDF::RNode& df);

template<int coordinate>
ROOT::VecOps::RVec<double> PrimariesVertex(const ROOT::VecOps::RVec<TG4PrimaryVertex>& vertices);

ROOT::VecOps::RVec<double> PrimariesVertexTimeDiff(ROOT::VecOps::RVec<double>& primaries_time);

int NofEventsPerSpill(const ROOT::VecOps::RVec<TG4PrimaryVertex>& vertices);

}//EDEPSIM

}//RDFUtils

#endif