#ifndef RDF_UTILS_H
#define RDF_UTILS_H

#include <string>

#include "TSystem.h"

#include "GenieUtils.h"
#include "GeoUtils.h"
// #include "SANDRecoUtils.h" // helix fitting

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

namespace TextColor {
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

static TextColor::Modifier def_(TextColor::FG_DEFAULT);
static TextColor::Modifier red_(TextColor::FG_RED); // warnings
static TextColor::Modifier green_(TextColor::FG_GREEN); // opeartions
static TextColor::Modifier blue_(TextColor::FG_BLUE); // results

inline void LOG(TString i, const char* out){
    if(i.Contains("I")){ // info
        std::cout << green_ << "[INFO] " << out << def_ << std::endl;
        std::cout << "\n";
    }else if(i.Contains("i")){
        std::cout << green_ << " |___ " << out << def_ << std::endl;
        std::cout << "\n";
    }else if(i.Contains("R")){ // result
        std::cout << blue_ << "[RESULT] " << out << def_ << std::endl;
        std::cout << "\n";
    }else if(i.Contains("W")){ // waring
        std::cout << red_ << "[WARNING] " << out << def_ << std::endl;
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

// bool IsUnphysical(genie::NtpMCEventRecord& m);

std::string InteractionTarget(const ROOT::VecOps::RVec<int>& pdg);

int NeutrinoFlavor(const ROOT::VecOps::RVec<int>& pdg);

std::string EventType(TObjString& t);

bool IsCC(TObjString& s);

ROOT::VecOps::RVec<genie::GHepParticle> AllGenieParticles(const ROOT::VecOps::RVec<int>& pdg,
                                                          const ROOT::VecOps::RVec<int>& status,
                                                          const ROOT::VecOps::RVec<double>& P4,
                                                          const ROOT::VecOps::RVec<double>& X4,
                                                          const ROOT::VecOps::RVec<int>& FirstMother,
                                                          const ROOT::VecOps::RVec<int>& LastMother,
                                                          const ROOT::VecOps::RVec<int>& FirstDaugther,
                                                          const ROOT::VecOps::RVec<int>& Lastdaugther);                                               

template<genie::GHepStatus_t STATUS>
ROOT::VecOps::RVec<genie::GHepParticle> GetParticlesWithStatus(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

template<int PDG>
ROOT::VecOps::RVec<genie::GHepParticle> GetParticlesWithPDG(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

template<int PDG>
ROOT::VecOps::RVec<genie::GHepParticle> ExcludePDG(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<std::string> PDG2Name(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<genie::GHepParticle> GetExhoticMesons(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<genie::GHepParticle> GetExhoticHadrons(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<genie::GHepParticle> GetRecoiledNuclei(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<genie::GHepParticle> FinalStateHadronicSystem(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

template<genie::GHepStatus_t STATUS>
ROOT::VecOps::RVec<genie::GHepParticle> PrimaryStateHadronicSystem_status(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<genie::GHepParticle> PrimaryStateHadronicSystem(const ROOT::VecOps::RVec<genie::GHepParticle>& particles, 
                                                                   std::string event_type);

ROOT::VecOps::RVec<genie::GHepParticle> GetFinalHadronicSystem(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

GenieUtils::event_topology GetInteractionTopology(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

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

// namespace SANDRECO{//SANDRECO

// ROOT::VecOps::RVec<particle> ParticleGoodTrackFit(const ROOT::VecOps::RVec<particle>& particles);

// ROOT::VecOps::RVec<particle> FilterPrimaries(const ROOT::VecOps::RVec<particle>& particles);

// template<int PDG>
// ROOT::VecOps::RVec<particle> GetParticlesWithPDG(const ROOT::VecOps::RVec<particle>& particles);

// ROOT::VecOps::RVec<TLorentzVector> GetMomentum(const ROOT::VecOps::RVec<particle>& particles);

// ROOT::VecOps::RVec<TLorentzVector> GetTrackVertex(const ROOT::VecOps::RVec<particle>& particles);

// ROOT::RDF::RNode AddColumnsFromSANDRECO(ROOT::RDF::RNode& df);

// }//SANDRECO

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