#ifndef RDF_UTILS_H
#define RDF_UTILS_H

#include <string>

#include "TSystem.h"

#include "GenieUtils.h"
#include "GeoUtils.h"
#include "EDepSimUtils.h"
#include "struct.h"
// #include "SANDRecoUtils.h" // helix fitting

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TString.h"
#include "TChain.h"
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

// ROOT::RDataFrame InitDF(TString production, const char* tree_name, unsigned int file_index_start = 0, unsigned int file_index_stop = 10);
TChain* InitTChain(TString production,
                            const char* tree_name,
                            unsigned int file_index_start,
                            unsigned int file_index_stop);

ROOT::RDataFrame InitDF(TChain* input_chain);

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

double SumDuble(const ROOT::VecOps::RVec<double>& v);

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

double GetNeutronTOF(std::string units_l,
                     std::string units_E,
                     const double vtx,
                     const double vty,
                     const double vtz,
                     TVector3 X3,
                     TVector3 P3
                     );

int CountChargedParticles(const genie::GHepParticle& fs_lepton, 
                          const ROOT::VecOps::RVec<genie::GHepParticle>& fs_HadronicSystem);

template<genie::GHepStatus_t STATUS>
ROOT::VecOps::RVec<genie::GHepParticle> GetParticlesWithStatus(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

template<int PDG>
ROOT::VecOps::RVec<genie::GHepParticle> GetParticlesWithPDG(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

template<int PDG>
ROOT::VecOps::RVec<genie::GHepParticle> ExcludePDG(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<std::string> PDG2Name(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

genie::GHepParticle GetIncomingNeutrino(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

genie::GHepParticle GetFinalStateLepton(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<genie::GHepParticle> GetInitialState(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<genie::GHepParticle> FinalStateHadronicSystem(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

genie::GHepParticle GetNucleonTarget(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

template<genie::GHepStatus_t STATUS>
ROOT::VecOps::RVec<genie::GHepParticle> PrimaryStateHadronicSystem_status(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<genie::GHepParticle> PrimaryStateHadronicSystem(const ROOT::VecOps::RVec<genie::GHepParticle>& particles, 
                                                                   std::string event_type);

GenieUtils::event_topology GetInteractionTopology(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<int> GetPDG(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<TLorentzVector> GetMomentum(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

ROOT::VecOps::RVec<double> GetKinE(const ROOT::VecOps::RVec<genie::GHepParticle>& particles);

TLorentzVector SumLorentzVectors(const ROOT::VecOps::RVec<TLorentzVector>& VTL);

int NofFinalStateParticles(const ROOT::VecOps::RVec<int>& pdg);

TVector3 GetTransverseComponent(const TLorentzVector& v1, const TLorentzVector v2);

ROOT::VecOps::RVec<double> DotProductWithAxis(const TVector3& axis, const ROOT::VecOps::RVec<TLorentzVector>& V);

TLorentzVector GetNup4FromMu(double proton_mass, double neutron_mass, double muon_mass,
                             const TLorentzVector& lepton, const TVector3 nu_direction);

TVector3 NeutronArrivalPosECAL(std::string units, double vtx_x, double vtx_y, double vtx_z, const TVector3& neutron_momentum);                             

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

std::string GetVolumeFromCoordinates(std::string, double x, double y, double z);

}//GEO

namespace EDEPSIM{//EDEPSIM

// TVector3 GetExpectedNeutronFromMu(const TLorentzVector& lepton);

// double GetBetaFromMomentum(const TVector3 p, double mass);

double NeutrinoEnergyFromCCQEonH(const TLorentzVector& muon, double muon_angle);

namespace SPILL{// EDEPSIM::SPILL
ROOT::RDF::RNode AddColumnsFromEDEPSIM(ROOT::RDF::RNode& df);

template<int coordinate>
ROOT::VecOps::RVec<double> PrimariesVertex(const ROOT::VecOps::RVec<TG4PrimaryVertex>& vertices);

ROOT::VecOps::RVec<double> PrimariesVertexTimeDiff(ROOT::VecOps::RVec<double>& primaries_time);

ROOT::VecOps::RVec<std::string> EventType(const ROOT::VecOps::RVec<TG4PrimaryVertex>& V);

ROOT::VecOps::RVec<std::string> FileName(const ROOT::VecOps::RVec<TG4PrimaryVertex>& V);

int NofEventsPerSpill(const ROOT::VecOps::RVec<TG4PrimaryVertex>& vertices);
} // SPILL

namespace NOSPILL{ // EDEPSIM::NOSPILL

ROOT::RDF::RNode AddColumnsFromEDEPSIM(ROOT::RDF::RNode& df);

std::string EventType(const ROOT::VecOps::RVec<TG4PrimaryVertex>& V);

bool IsCCQEonH(const ROOT::VecOps::RVec<TG4PrimaryVertex>& V);

std::string FileName(const ROOT::VecOps::RVec<TG4PrimaryVertex>& V);

int NofPrimaries(const ROOT::VecOps::RVec<TG4PrimaryVertex>& v);

ROOT::VecOps::RVec<int> GetPrimariesPDG(const ROOT::VecOps::RVec<TG4PrimaryVertex>& v);

ROOT::VecOps::RVec<TLorentzVector> GetPrimariesP4(const ROOT::VecOps::RVec<TG4PrimaryVertex>& v);

TLorentzVector GetPrimaryHadronicP4(const ROOT::VecOps::RVec<TLorentzVector>& momenta, const ROOT::VecOps::RVec<int> pdgs);

TLorentzVector GetPrimaryLeptonP4(const ROOT::VecOps::RVec<TLorentzVector>& momenta, const ROOT::VecOps::RVec<int> pdgs);

ROOT::VecOps::RVec<double> GetEmissionAngle(const TVector3 nuDirection, const ROOT::VecOps::RVec<TLorentzVector>& primariesP4);

ROOT::VecOps::RVec<int> GetPrimariesTrackId(const ROOT::VecOps::RVec<TG4PrimaryVertex>& v);

ROOT::VecOps::RVec<std::string> PDG2Name(const ROOT::VecOps::RVec<int>& pdgs);

ROOT::VecOps::RVec<EDepUtils::track_hits> GetPrimariesHits(TG4Event& ev, const ROOT::VecOps::RVec<int>& ids, const ROOT::VecOps::RVec<int>& pdgs);

ROOT::VecOps::RVec<ROOT::VecOps::RVec<TG4HitSegment>> GroupHitsByPrimary(TG4Event& ev, const ROOT::VecOps::RVec<int>& ids);

ROOT::VecOps::RVec<double> GetPrimariesEDepECAL(const ROOT::VecOps::RVec<EDepUtils::track_hits>& vector_track_hits);

ROOT::VecOps::RVec<TLorentzVector> GetPrimariesFirstHitECAL(const ROOT::VecOps::RVec<EDepUtils::track_hits>& vector_track_hits);

template<int coordinate>
ROOT::VecOps::RVec<double> GetECALHitPos(const ROOT::VecOps::RVec<EDepUtils::track_hits>& vector_primary_hits);

} // NOSPILL

}//EDEPSIM

namespace DIGIT{// DIGIT

ROOT::RDF::RNode AddColumnsFromDigit(ROOT::RDF::RNode& df);

int NofFiredECALMods(const ROOT::VecOps::RVec<int>& fired_cells_modules);

ROOT::VecOps::RVec<TLorentzVector> ReconstructHitFromCell(const ROOT::VecOps::RVec<dg_cell>& cells);

ROOT::VecOps::RVec<int> IsCellComplete(const ROOT::VecOps::RVec<dg_cell>& cells);

std::vector<int> FiredECALMods(const ROOT::VecOps::RVec<int>& fired_cells_modules);

template<int side>
ROOT::VecOps::RVec<double> FiredECALGetTDC(const ROOT::VecOps::RVec<dg_cell>& cells);

int NofClusters(const ROOT::VecOps::RVec<cluster>& clusters);

ROOT::VecOps::RVec<TLorentzVector> GetClusterX4(const ROOT::VecOps::RVec<cluster>& clusters);

ROOT::VecOps::RVec<TLorentzVector> Cluster2Vertex4Distance(double x,
                                                           double y,
                                                           double z,
                                                           double t,
                                                           const ROOT::VecOps::RVec<cluster>& clusters);

ROOT::VecOps::RVec<double> STDistance(const TVector3& expected_neutron_hit3,
                                      const double exp_neutron_tof,
                                      const ROOT::VecOps::RVec<TLorentzVector>& cells_reco_hits);                                                       

}// DIGIT

}//RDFUtils

#endif