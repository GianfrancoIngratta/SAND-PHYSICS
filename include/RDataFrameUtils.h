#ifndef RDF_UTILS_H
#define RDF_UTILS_H

#include <string>

#include "TSystem.h"

#include "GenieUtils.h"
#include "GeoUtils.h"
#include "EDepSimUtils.h"
#include "struct.h"
#include "SANDRecoUtils.h" // helix fitting

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

/** ==================
 * RDF:
 */ 
namespace RDFUtils{

const int DEFAULT_NO_DATA = -999;

const TLorentzVector DEFAULT_VECTOR_NO_DATA = {DEFAULT_NO_DATA, DEFAULT_NO_DATA, DEFAULT_NO_DATA, DEFAULT_NO_DATA}; 

struct trajectory
{   
    // /*** see : TG4PrimaryVertex
    //  * @attribute: GeneratorName
    //  * The name of the generator that created this trajecotry's vertex.
    //  */
    // std::string GeneratorName;
    std::string file_edep;

    int spill_number = DEFAULT_NO_DATA;
    
    int vertex_interaction_number = DEFAULT_NO_DATA;

    int id = DEFAULT_NO_DATA;

    int pdg = DEFAULT_NO_DATA;

    std::string name;
    
    /***
     * @attribute: start_point
     * trajectory starting point
     */
    TLorentzVector start_point;
    
    /***
     * @attribute: start_point
     * trajectory starting point
     */
    std::string start_volume;
    
    /***
     * @attribute: initial_momentum
     * trajectory initial momentum
     */
    TLorentzVector initial_momentum;
    
    /***
     * @attribute: hindex_ecal
     * list of all track hits in ECAL
     */
    std::vector<int> hindex_ecal;

    /**
     * @attribute:earliest_hit
     * earliest hit in ecal
     */
    TLorentzVector earliest_hit = DEFAULT_VECTOR_NO_DATA;

    /**
     * @attribute:earliest_hit
     * latest hit in ecal
     */
    TLorentzVector latest_hit = DEFAULT_VECTOR_NO_DATA;

    /***
     * @attribute: edep_in_ecal
     * energy deposited in ecal
     */
    double edep_in_ecal = 0.;

    /***
     * @attribute: cell_id_ecal
     * list of all ecal cells whose TDC
     * is produced by this trajectory
     */
    std::vector<int> cell_id_ecal;

    /***
     * @attribute: tof_to_ecal
     * true trajectory tof to ecal from
     * vertex to first hit in ecal
     */
    double tof_to_ecal;

};

struct cell
{
    std::string file_edep;

    int spill_number = DEFAULT_NO_DATA;
    
    int vertex_interaction_number = DEFAULT_NO_DATA;

    int id = DEFAULT_NO_DATA;

    /**
     * @attribute: cell has both tdc with signal
     * -1 uninitialized
     * 0 broken, 1 complete
     */
    int is_complete = -1;
    
    /**
     * @attribute: track_id of the trajectory
     * that whose hit crated the observed tdc on pmt1 (pmt2)
     */
    int track_id_pmt1_hit;
    int track_id_pmt2_hit;

    /**
     * @attribute: pdg of the trajectory
     * that whose hit crated the observed tdc on pmt1 (pmt2)
     */
    int track_pdg_pmt1_hit;
    int track_pdg_pmt2_hit;

    /**
     * @attribute: edepsim hit that produced the 
     * observed TDC on pmt1 and pmt2
     */ 
    TLorentzVector true_hit1 = DEFAULT_VECTOR_NO_DATA;
    TLorentzVector true_hit2 = DEFAULT_VECTOR_NO_DATA;

    /**
     * @attribute: true_energy deposited in the cell
     */
    double true_edep1 = DEFAULT_NO_DATA;
    double true_edep2 = DEFAULT_NO_DATA;

    /**
     * @attribute: reco_hit reconstructed hit from obseved
     * TDCs
     */
    TLorentzVector reco_hit = DEFAULT_VECTOR_NO_DATA;

    /**
     * @attribute: reconstructed energy deposited in the cell
     */
    double reco_edep = DEFAULT_NO_DATA;
};

/**
 * TRAJECTORY FUNCTIONS: 
 */
ROOT::RDF::RNode CreateDataFrameTrajectories(ROOT::RDF::RNode& df);
ROOT::VecOps::RVec<RDFUtils::trajectory> CreateTrajectories(TG4Event& ev, 
                                                            const ROOT::VecOps::RVec<TG4Trajectory>& TG4_trajectories,
                                                            const ROOT::VecOps::RVec<dg_cell>& cells,
                                                            const ROOT::VecOps::RVec<int>& track_id_tdc1,
                                                            const ROOT::VecOps::RVec<int>& track_id_tdc2);
ROOT::VecOps::RVec<std::string> GetTrajectoryFileName(const ROOT::VecOps::RVec<RDFUtils::trajectory>& trajectoires);
ROOT::VecOps::RVec<std::string> GetTrajectoryName(const ROOT::VecOps::RVec<RDFUtils::trajectory>& trajectoires);
ROOT::VecOps::RVec<int> GetTrajectoryId(const ROOT::VecOps::RVec<RDFUtils::trajectory>& trajectoires);
ROOT::VecOps::RVec<int> GetTrajectorySpillNumber(const ROOT::VecOps::RVec<RDFUtils::trajectory>& trajectoires);
ROOT::VecOps::RVec<int> GetTrajectoryPDG(const ROOT::VecOps::RVec<RDFUtils::trajectory>& trajectoires);
ROOT::VecOps::RVec<std::string> GetTrajectoryName(const ROOT::VecOps::RVec<RDFUtils::trajectory>& trajectoires);
ROOT::VecOps::RVec<std::string> GetTrajectoryVolume(const ROOT::VecOps::RVec<RDFUtils::trajectory>& trajectoires);
ROOT::VecOps::RVec<double> GetTrajectoryECALedep(const ROOT::VecOps::RVec<RDFUtils::trajectory>& trajectoires);
ROOT::VecOps::RVec<TLorentzVector> GetTrajectoryECALearliest(const ROOT::VecOps::RVec<RDFUtils::trajectory>& trajectoires);
ROOT::VecOps::RVec<TLorentzVector> GetTrajectoryECALlatest(const ROOT::VecOps::RVec<RDFUtils::trajectory>& trajectoires);
ROOT::VecOps::RVec<double> GetTrajectoryTOF2ECAL(const ROOT::VecOps::RVec<RDFUtils::trajectory>& trajectoires);

/**
 * CELLS FUNCTIONS: 
 */
ROOT::RDF::RNode CreateDataFrameCells(ROOT::RDF::RNode& df);
ROOT::VecOps::RVec<RDFUtils::cell> CreateCells(TG4Event& ev, 
                                               const ROOT::VecOps::RVec<TG4Trajectory>& TG4_trajectories,
                                               const ROOT::VecOps::RVec<dg_cell>& input_cells,
                                               double calibration_x,
                                               double calibration_y,
                                               double calibration_t);
ROOT::VecOps::RVec<int> GetCellId(const ROOT::VecOps::RVec<RDFUtils::cell>& cells);
ROOT::VecOps::RVec<std::string> GetCellFileName(const ROOT::VecOps::RVec<RDFUtils::cell>& cells);
ROOT::VecOps::RVec<int> GetCellSpillNumber(const ROOT::VecOps::RVec<RDFUtils::cell>& cells);
ROOT::VecOps::RVec<int> IsCellComplete(const ROOT::VecOps::RVec<RDFUtils::cell>& cells);
ROOT::VecOps::RVec<int> GetCellTrackId1(const ROOT::VecOps::RVec<RDFUtils::cell>& cells);
ROOT::VecOps::RVec<int> GetCellTrackId2(const ROOT::VecOps::RVec<RDFUtils::cell>& cells);
ROOT::VecOps::RVec<int> GetCellTrackPDG1(const ROOT::VecOps::RVec<RDFUtils::cell>& cells);
ROOT::VecOps::RVec<int> GetCellTrackPDG2(const ROOT::VecOps::RVec<RDFUtils::cell>& cells);
ROOT::VecOps::RVec<TLorentzVector> GetCellTrueHit1(const ROOT::VecOps::RVec<RDFUtils::cell>& cells);
ROOT::VecOps::RVec<TLorentzVector> GetCellTrueHit2(const ROOT::VecOps::RVec<RDFUtils::cell>& cells);
ROOT::VecOps::RVec<TLorentzVector> GetCellRecoHit(const ROOT::VecOps::RVec<RDFUtils::cell>& cells);
ROOT::VecOps::RVec<double> GetCellTrueEdep1(const ROOT::VecOps::RVec<RDFUtils::cell>& cells);
ROOT::VecOps::RVec<double> GetCellTrueEdep2(const ROOT::VecOps::RVec<RDFUtils::cell>& cells);
ROOT::VecOps::RVec<double> GetCellRecoEdep(const ROOT::VecOps::RVec<RDFUtils::cell>& cells);

TChain* InitTChain(TString production,
                            const char* tree_name,
                            unsigned int file_index_start,
                            unsigned int file_index_stop);

ROOT::RDataFrame InitDF(TChain* input_chain);

ROOT::RDF::RNode Filter(ROOT::RDF::RNode& df, const char* condition, bool report);

void PrintColumns(ROOT::RDataFrame& df);

ROOT::RDF::RNode AddConstantsToDF(ROOT::RDataFrame& df);

template<int coord>
ROOT::VecOps::RVec<double> GetComponent(const ROOT::VecOps::RVec<TLorentzVector>& vTL);

template<int coord>
ROOT::VecOps::RVec<double> GetComponent3(const ROOT::VecOps::RVec<TVector3>& vTV);

ROOT::VecOps::RVec<TVector3> GetVect(const ROOT::VecOps::RVec<TLorentzVector>& vTL);

ROOT::VecOps::RVec<int> FindMinimum(const ROOT::VecOps::RVec<double> v);

ROOT::VecOps::RVec<int> FindMinimum2(const ROOT::VecOps::RVec<double> vector,
                                               const ROOT::VecOps::RVec<int> condition);

double FindMinimum3(const ROOT::VecOps::RVec<double> values, const ROOT::VecOps::RVec<int> condition);

double slice_w_condition(const ROOT::VecOps::RVec<double>& values,
                                   const ROOT::VecOps::RVec<double>& mask,
                                   const double value);

double GetColumnSum(const ROOT::VecOps::RVec<double>& v);

ROOT::VecOps::RVec<double> VectorDifference(const ROOT::VecOps::RVec<double>& v1,
                                            const ROOT::VecOps::RVec<double>& v2);

TVector3 TLVectorCrossProduct(const TLorentzVector& v1,
                              const TLorentzVector& v2);

ROOT::VecOps::RVec<double> VectorSubtractConst(const ROOT::VecOps::RVec<double>& v1, double c);

TLorentzVector VectorFilterByHighest(const ROOT::VecOps::RVec<double>& filter,
                                     const ROOT::VecOps::RVec<TLorentzVector>& v);

double SumDuble(const ROOT::VecOps::RVec<double>& v);

ROOT::VecOps::RVec<int> CompareVectors(const ROOT::VecOps::RVec<double>& V1, 
                                        const ROOT::VecOps::RVec<double>& V2);

ROOT::VecOps::RVec<int> IsNULLTLV(const ROOT::VecOps::RVec<TLorentzVector>& V);

ROOT::VecOps::RVec<int> IsValidVector(const ROOT::VecOps::RVec<double>& v);

/** ==================
 * GENIE:
 */ 

namespace GENIE{

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

ROOT::VecOps::RVec<double> GetBeta(const ROOT::VecOps::RVec<TLorentzVector>& P4);

TLorentzVector SumLorentzVectors(const ROOT::VecOps::RVec<TLorentzVector>& VTL);

int NofFinalStateParticles(const ROOT::VecOps::RVec<int>& pdg);

TVector3 GetTransverseComponent(const TLorentzVector& v1, const TLorentzVector v2);

ROOT::VecOps::RVec<double> DotProductWithAxis(const TVector3& axis, const ROOT::VecOps::RVec<TLorentzVector>& V);

TLorentzVector GetNup4FromMu(double proton_mass, double neutron_mass, double muon_mass,
                             const TLorentzVector& lepton, const TVector3 nu_direction);

TVector3 NeutronArrivalPosECAL(std::string units, double vtx_x, double vtx_y, double vtx_z, const TVector3& neutron_momentum);                             

ROOT::RDF::RNode AddColumnsFromGENIE(ROOT::RDF::RNode& df);

ROOT::RDF::RNode AddColumnsForHydrogenCarbonSampleSelection(ROOT::RDF::RNode& df); // these is expected to have columns from AddColumnsFromGENIE output

}
/**
 * GENIE:
 ==================*/ 

/** =================
 * GEO:
 * */ 

namespace GEO{

std::string GetMaterialFromCoordinates(double x, double y, double z);

std::string GetVolumeFromCoordinates(std::string, double x, double y, double z);

ROOT::VecOps::RVec<std::string> GetVolumeFromCoordinates_v(std::string units,
                                                           const ROOT::VecOps::RVec<double>& vertices_x, 
                                                           const ROOT::VecOps::RVec<double>& vertices_y, 
                                                           const ROOT::VecOps::RVec<double>& vertices_z);

}
/**
 * GEO:
 ==================*/ 

/** =================
 * EDEPSIM:
*/ 
namespace EDEPSIM{

// TVector3 GetExpectedNeutronFromMu(const TLorentzVector& lepton);

// double GetBetaFromMomentum(const TVector3 p, double mass);

double NeutrinoEnergyFromCCQEonH(const TLorentzVector& muon, double muon_angle);

ROOT::VecOps::RVec<int> GetHitTrajectoryId(const ROOT::VecOps::RVec<int>& pmts_hindex, 
                                 TG4Event& ev);

template<int PDG>
ROOT::VecOps::RVec<int> CheckTrajId(const ROOT::VecOps::RVec<int>& track_ids_pmt1,
                                    const ROOT::VecOps::RVec<int>& track_ids_pmt2,
                                    TG4Event& ev);

ROOT::VecOps::RVec<double> GetNeutronEKin(const ROOT::VecOps::RVec<int>& track_ids,
                                         const ROOT::VecOps::RVec<int>& is_primary_neutron, 
                                         TG4Event& ev,
                                         const double neutron_mass);

ROOT::VecOps::RVec<TG4Trajectory> GetTrajectories(const ROOT::VecOps::RVec<int>& track_ids_pmt1, 
                                                  const ROOT::VecOps::RVec<int>& track_ids_pmt2,
                                                  TG4Event& ev);

ROOT::VecOps::RVec<int> ExpandIndex(const ROOT::VecOps::RVec<TG4Trajectory>& trajectories);

ROOT::VecOps::RVec<int> ExpandPDG(const ROOT::VecOps::RVec<TG4Trajectory>& trajectories);

ROOT::VecOps::RVec<TG4TrajectoryPoint> GetTrjPointsAll(const ROOT::VecOps::RVec<TG4Trajectory>& trajectories);

template<int PDG>
ROOT::VecOps::RVec<double> GetPointComponent(const ROOT::VecOps::RVec<TG4TrajectoryPoint>& points);

template<int PDG>
ROOT::VecOps::RVec<double> GetPointMomentum(const ROOT::VecOps::RVec<TG4TrajectoryPoint>& points);

ROOT::VecOps::RVec<double> GetPointProcess(const ROOT::VecOps::RVec<TG4TrajectoryPoint>& points);

ROOT::VecOps::RVec<TLorentzVector> GetHitFromIndex(const ROOT::VecOps::RVec<int>& h_index, TG4Event& ev);

ROOT::VecOps::RVec<double> GetHitEdepFromIndex(const ROOT::VecOps::RVec<int>& h_index, TG4Event& ev);

/** =================
 * EDEPSIM::SPILL:
*/ 
namespace SPILL{
ROOT::RDF::RNode AddColumnsFromEDEPSIM(ROOT::RDF::RNode& df);

template<int coordinate>
ROOT::VecOps::RVec<double> PrimariesVertex(const ROOT::VecOps::RVec<TG4PrimaryVertex>& vertices);

// ROOT::VecOps::RVec<double> PrimariesVertexTimeDiff(ROOT::VecOps::RVec<double>& primaries_time);

// ROOT::VecOps::RVec<std::string> EventType(const ROOT::VecOps::RVec<TG4PrimaryVertex>& V);

template<int token_number>
ROOT::VecOps::RVec<std::string> EventInfo(const ROOT::VecOps::RVec<TG4PrimaryVertex>& vertices);

ROOT::VecOps::RVec<std::string> FileName(const ROOT::VecOps::RVec<TG4PrimaryVertex>& vertex);

int NofEventsPerSpill(const ROOT::VecOps::RVec<TG4PrimaryVertex>& vertices);

ROOT::VecOps::RVec<int> GetBatchNumber(const ROOT::VecOps::RVec<double>& vertex_time, 
                                                                 const double batch_duration,
                                                                 const int nof_batches_in_spill,
                                                                 const double spill_t0);

ROOT::VecOps::RVec<int> GetBunchNumber(const ROOT::VecOps::RVec<double>& vertex_time,
                                                                const double batch_duration,
                                                                const double bunch_separation,
                                                                const int nof_bunches_in_batch,
                                                                const double spill_t0);

ROOT::VecOps::RVec<std::string> GetEventType(const ROOT::VecOps::RVec<TG4PrimaryVertex>& vertices);

ROOT::VecOps::RVec<int> IsSignal(const ROOT::VecOps::RVec<TG4PrimaryVertex>& vertices); 

ROOT::VecOps::RVec<int> CountChargedPrimaries(const ROOT::VecOps::RVec<TG4PrimaryVertex>& vertices);

}
/**
 * EDEPSIM::SPILL:
 ==================*/ 
/** =================
 * EDEPSIM::NOSPILL:
 */ 

namespace NOSPILL{

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

ROOT::VecOps::RVec<double> GetDirectionInECAL(const ROOT::VecOps::RVec<TLorentzVector>& P4,
                                              const double vtxX,
                                              const double vtxY,
                                              const double vtxZ,
                                              const ROOT::VecOps::RVec<TLorentzVector>& primaries_first_hit);

}
/** 
 * EDEPSIM::NOSPILL:
 ==================*/ 

}
/** 
 * EDEPSIM:
 ==================*/ 

/** =================
 * DIGIT:
*/ 
namespace DIGIT{// DIGIT

ROOT::RDF::RNode AddColumnsFromDigit(ROOT::RDF::RNode& df);

ROOT::RDF::RNode GetInfoCellsFromSignal(ROOT::RDF::RNode& df);

ROOT::RDF::RNode GetFilteredTrajectories(ROOT::RDF::RNode& df, std::string column_name_TG4_trajectories);

int NofFiredECALMods(const ROOT::VecOps::RVec<int>& fired_cells_modules);

ROOT::VecOps::RVec<TVector3> XfromTDC(const ROOT::VecOps::RVec<dg_cell>& cells,
                                      const double calibration_x,
                                      const double calibration_y);

ROOT::VecOps::RVec<double> TfromTDC(const ROOT::VecOps::RVec<dg_cell>& cells, 
                                    const double calibration_const);

ROOT::VecOps::RVec<int> CheckTDCsConsistency(const ROOT::VecOps::RVec<dg_cell>& cells,
                                                              const ROOT::VecOps::RVec<TVector3> x_hit_reco, 
                                                              const ROOT::VecOps::RVec<double> t_hit_reco);

ROOT::VecOps::RVec<double> EfromTDC(const ROOT::VecOps::RVec<dg_cell>& cells,
                                                     const double calibration_const,
                                                     const ROOT::VecOps::RVec<double>& e_true);

TVector3 WeightedCell(const ROOT::VecOps::RVec<TVector3> x_reco,
                                             const ROOT::VecOps::RVec<double> e_reco,
                                             const ROOT::VecOps::RVec<int> fired_by_primary);

double WeightedCellTime(const ROOT::VecOps::RVec<double>& t_hit_reco,
                        const ROOT::VecOps::RVec<double>& e_reco,
                        const ROOT::VecOps::RVec<int>& isCompatible);

ROOT::VecOps::RVec<double> GetFlightLength(const TVector3& vtx, ROOT::VecOps::RVec<TVector3>& hits);


ROOT::VecOps::RVec<double> GetTOF(const ROOT::VecOps::RVec<double>& flight_length, const double beta);

ROOT::VecOps::RVec<int> IsCellComplete(const ROOT::VecOps::RVec<dg_cell>& cells);

std::vector<int> FiredECALMods(const ROOT::VecOps::RVec<int>& fired_cells_modules);

int Get_TDC_hindex(const dg_ps& photo_sensor);

template<int side>
ROOT::VecOps::RVec<double> FiredECALGetTDC(const ROOT::VecOps::RVec<dg_cell>& cells);

template<int side>
ROOT::VecOps::RVec<double> FiredECALGetADC(const ROOT::VecOps::RVec<dg_cell>& cells);

template<int side>
ROOT::VecOps::RVec<int> GetHindexOfTDC(const ROOT::VecOps::RVec<dg_cell>& cells);

int NofClusters(const ROOT::VecOps::RVec<cluster>& clusters);

ROOT::VecOps::RVec<TLorentzVector> GetClusterX4(const ROOT::VecOps::RVec<cluster>& clusters);

ROOT::VecOps::RVec<TLorentzVector> Cluster2Vertex4Distance(double x,
                                                           double y,
                                                           double z,
                                                           double t,
                                                           const ROOT::VecOps::RVec<cluster>& clusters);

ROOT::VecOps::RVec<TVector3> GetExpectedHitPosition(TVector3 vertex,
                                                    const TVector3& momentum_vector,
                                                    const ROOT::VecOps::RVec<TVector3>& hits_reco);

ROOT::VecOps::RVec<double> TimeResiduals(const ROOT::VecOps::RVec<double>& expected_t_hit,
                                                          const ROOT::VecOps::RVec<double>& reconstruced_t_hit);

ROOT::VecOps::RVec<double> SpaceResiduals(const ROOT::VecOps::RVec<TVector3>& expected_x_hit,
                                                           const ROOT::VecOps::RVec<TVector3>& reconstruced_x_hit);


ROOT::VecOps::RVec<int> IsCompatible(const ROOT::VecOps::RVec<double>& space_residuals,
                                     const ROOT::VecOps::RVec<double>& time_residuals,
                                     const ROOT::VecOps::RVec<int>& are_tdc_cell_consistent);

ROOT::VecOps::RVec<int> IsCompatible2(const ROOT::VecOps::RVec<double>& space_residuals,
                                      const ROOT::VecOps::RVec<double>& time_residuals,
                                      const ROOT::VecOps::RVec<double>& reco_energy,
                                      const ROOT::VecOps::RVec<int>& are_tdc_cell_consistent);

ROOT::VecOps::RVec<int> IsCandidateCell(const ROOT::VecOps::RVec<int>& isSpaceCompatible,
                                        const ROOT::VecOps::RVec<double>& t_hit_reco);

int isCandidateSignal(ROOT::VecOps::RVec<int>& cell_has_coincidence);


ROOT::VecOps::RVec<double> GetBeta(const ROOT::VecOps::RVec<double>& flight_length,
                                    const ROOT::VecOps::RVec<double>& t_hit_reco,
                                    double vertex_time);

ROOT::VecOps::RVec<double> GetPFromBeta(const ROOT::VecOps::RVec<double>& beta,
                                                         const double neutron_mass);

ROOT::VecOps::RVec<double> GetKinEFromBeta(const ROOT::VecOps::RVec<double>& beta,
                                                         const double neutron_mass);

ROOT::VecOps::RVec<TVector3> GetRecoP3(const ROOT::VecOps::RVec<double>& ptot_reco,
                                       const TVector3& expected_direction);

                                       
ROOT::VecOps::RVec<double> GetMissingPT(const TLorentzVector& lepton_P4,
                                        const ROOT::VecOps::RVec<TVector3>& neutron_reco_p3);


ROOT::VecOps::RVec<dg_cell> FitlerCells(const ROOT::VecOps::RVec<dg_cell>& cells,
                                        const ROOT::VecOps::RVec<int>& is_complete,
                                        const ROOT::VecOps::RVec<int>& is_fired_by_neutron);

ROOT::VecOps::RVec<int>  GetCellMod(const ROOT::VecOps::RVec<dg_cell>& cells);                                        

ROOT::VecOps::RVec<int>  GetCellId(const ROOT::VecOps::RVec<dg_cell>& cells);

ROOT::VecOps::RVec<double>  GetCellX(const ROOT::VecOps::RVec<dg_cell>& cells);

ROOT::VecOps::RVec<double>  GetCellY(const ROOT::VecOps::RVec<dg_cell>& cells);

ROOT::VecOps::RVec<double>  GetCellZ(const ROOT::VecOps::RVec<dg_cell>& cells);

/** =================
 * DIGIT::SPILL:
*/ 

namespace SPILL
{

}
/**
 * DIGIT::SPILL:
 ==================*/

}
/**
 * DIGIT:
 ==================*/

/** =================
 * RECO:
*/
namespace RECO{

int GetNofFiredWires(const ROOT::VecOps::RVec<int>& wires_id);

int Test(bool b, const MinuitFitInfos& i);

ROOT::RDF::RNode AddColumnsFromDriftReco(ROOT::RDF::RNode& df);

int isCandidateSignal(int charge_multiplicity, int nof_cell_candidates);

}
/**
 * RECO:
 ==================*/
}

/**
 * RDF:
 ==================*/

#endif