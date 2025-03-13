#include <string>
#include <vector>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

// source /opt/exp_software/neutrino/ROOUNFOLD/RooUnfold/setup.sh
// .L /opt/exp_software/neutrino/ROOUNFOLD/RooUnfold/src/RooUnfold.h

#include <TEfficiency.h>
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

enum UNFOLD_COMPONTENTS{
    SIGNAL_SAMPLE1 = 0,
    SIGNAL_BKG_SAMPLE1 = 1,
    BKG_SAMPLE1 = 2,
    //
    SIGNAL_SAMPLE2 = 3,
    SIGNAL_BKG_SAMPLE2 = 4,
    BKG_SAMPLE2 = 5,
    //
    SIGNAL_SAMPLE3 = 6,
    SIGNAL_BKG_SAMPLE3 = 7,
    BKG_SAMPLE3 = 8,
    UNFOLD_NONE
};


/**
UNFOLD:
 */
TFile* unfold_file = new TFile("unfold_file.root", "RECREATE");

RooUnfoldResponse* response_matrix[UNFOLD_NONE];

TEfficiency* efficiency[UNFOLD_NONE];

TH1D* hist4unfold[UNFOLD_NONE][3];

uint nof_bin_nu_E = 12;
double low_bin_nu_E = 0.;
double up_bin_nu_E = 6.;

void Init_Unfold(){
    // ALL EFFICIENCIES
    efficiency[SIGNAL_SAMPLE1] = new TEfficiency("eff_signal_sample1", "Efficiency Signal Sample 1", nof_bin_nu_E, 0., 6.);
    efficiency[SIGNAL_BKG_SAMPLE1] = new TEfficiency("eff_signal_bkg_sample1", "Efficiency Signal + Background Sample 1", nof_bin_nu_E, 0., 6.);
    efficiency[BKG_SAMPLE1] = new TEfficiency("eff_bkg_sample1", "Efficiency Background Sample 1", nof_bin_nu_E, 0., 6.);
    // EFFICIECNY_SAMPLE2
    efficiency[SIGNAL_SAMPLE2] = new TEfficiency("eff_signal_sample2", "Efficiency Signal Sample 2", nof_bin_nu_E, 0., 6.);
    efficiency[SIGNAL_BKG_SAMPLE2] = new TEfficiency("eff_signal_bkg_sample2", "Efficiency Signal + Background Sample 2", nof_bin_nu_E, 0., 6.);
    efficiency[BKG_SAMPLE2] = new TEfficiency("eff_bkg_sample2", "Efficiency Background Sample 2", nof_bin_nu_E, 0., 6.);

    // ALL RESPONSES
    // MARTRIX RESPONSE SAMPLE 1
    response_matrix[SIGNAL_SAMPLE1] = new RooUnfoldResponse(nof_bin_nu_E,0.,6.,"SIGNAL_SAMPLE1","SIGNAL_SAMPLE1");
    response_matrix[SIGNAL_BKG_SAMPLE1] = new RooUnfoldResponse(nof_bin_nu_E,0.,6.,"SIGNAL_BKG_SAMPLE1","SIGNAL_BKG_SAMPLE1");
    response_matrix[BKG_SAMPLE1] = new RooUnfoldResponse(nof_bin_nu_E,0.,6.,"BKG_SAMPLE1","BKG_SAMPLE1");
    // MARTRIX RESPONSE SAMPLE 2
    response_matrix[SIGNAL_SAMPLE2] = new RooUnfoldResponse(nof_bin_nu_E,0.,6.,"SIGNAL_SAMPLE2","SIGNAL_SAMPLE2");
    response_matrix[SIGNAL_BKG_SAMPLE2] = new RooUnfoldResponse(nof_bin_nu_E,0.,6.,"SIGNAL_BKG_SAMPLE2","SIGNAL_BKG_SAMPLE2");
    response_matrix[BKG_SAMPLE2] = new RooUnfoldResponse(nof_bin_nu_E,0.,6.,"BKG_SAMPLE2","BKG_SAMPLE2");

    // ALL HISTOS
    hist4unfold[SIGNAL_SAMPLE1][0] = new TH1D("h_SIGNAL_SAMPLE1_true", "h_SIGNAL_SAMPLE1_true", nof_bin_nu_E, 0, 6);
    hist4unfold[SIGNAL_SAMPLE1][1] = new TH1D("h_SIGNAL_SAMPLE1_reco", "h_SIGNAL_SAMPLE1_reco", nof_bin_nu_E, 0, 6);
    hist4unfold[SIGNAL_BKG_SAMPLE1][0] = new TH1D("h_SIGNAL_BKG_SAMPLE1_true", "h_SIGNAL_BKG_SAMPLE1_true", nof_bin_nu_E, 0, 6);
    hist4unfold[SIGNAL_BKG_SAMPLE1][1] = new TH1D("h_SIGNAL_BKG_SAMPLE1_reco", "h_SIGNAL_BKG_SAMPLE1_reco", nof_bin_nu_E, 0, 6);
    hist4unfold[BKG_SAMPLE1][0] = new TH1D("h_BKG_SAMPLE1_true", "h_BKG_SAMPLE1_true", nof_bin_nu_E, 0, 6);
    hist4unfold[BKG_SAMPLE1][1] = new TH1D("h_BKG_SAMPLE1_reco", "h_BKG_SAMPLE1_reco", nof_bin_nu_E, 0, 6);
    
    hist4unfold[SIGNAL_SAMPLE2][0] = new TH1D("h_SIGNAL_SAMPLE2_true", "h_SIGNAL_SAMPLE2_true", nof_bin_nu_E, 0, 6);
    hist4unfold[SIGNAL_SAMPLE2][1] = new TH1D("h_SIGNAL_SAMPLE2_reco", "h_SIGNAL_SAMPLE2_reco", nof_bin_nu_E, 0, 6);
    hist4unfold[SIGNAL_BKG_SAMPLE2][0] = new TH1D("h_SIGNAL_BKG_SAMPLE2_true", "h_SIGNAL_BKG_SAMPLE2_true", nof_bin_nu_E, 0, 6);
    hist4unfold[SIGNAL_BKG_SAMPLE2][1] = new TH1D("h_SIGNAL_BKG_SAMPLE2_reco", "h_SIGNAL_BKG_SAMPLE2_reco", nof_bin_nu_E, 0, 6);
    hist4unfold[BKG_SAMPLE2][0] = new TH1D("h_BKG_SAMPLE2_true", "h_BKG_SAMPLE2_true", nof_bin_nu_E, 0, 6);
    hist4unfold[BKG_SAMPLE2][1] = new TH1D("h_BKG_SAMPLE2_reco", "h_BKG_SAMPLE2_reco", nof_bin_nu_E, 0, 6);
    
    hist4unfold[SIGNAL_SAMPLE3][0] = new TH1D("h_SIGNAL_SAMPLE3_true", "h_SIGNAL_SAMPLE3_true", nof_bin_nu_E, 0, 6);
    hist4unfold[SIGNAL_SAMPLE3][1] = new TH1D("h_SIGNAL_SAMPLE3_reco", "h_SIGNAL_SAMPLE3_reco", nof_bin_nu_E, 0, 6);
    hist4unfold[SIGNAL_BKG_SAMPLE3][0] = new TH1D("h_SIGNAL_BKG_SAMPLE3_true", "h_SIGNAL_BKG_SAMPLE3_true", nof_bin_nu_E, 0, 6);
    hist4unfold[SIGNAL_BKG_SAMPLE3][1] = new TH1D("h_SIGNAL_BKG_SAMPLE3_reco", "h_SIGNAL_BKG3_reco", nof_bin_nu_E, 0, 6);
    hist4unfold[BKG_SAMPLE3][0] = new TH1D("h_BKG_SAMPLE3_true", "h_BKG_SAMPLE3_true", nof_bin_nu_E, 0, 6);
    hist4unfold[BKG_SAMPLE3][1] = new TH1D("h_BKG_SAMPLE3_reco", "h_BKG_SAMPLE3_reco", nof_bin_nu_E, 0, 6);
}

void PrintAllUnfoldInfos()
{
    unfold_file->cd();
    for(size_t i=0; i< SIGNAL_SAMPLE3; i++){
        efficiency[i] -> Write();
        response_matrix[i] -> Write();
    }
    for(size_t i=0; i< UNFOLD_NONE; i++){
        hist4unfold[i][0] -> Write(); 
        hist4unfold[i][1] -> Write();
    }
    unfold_file->Close();
}

struct CellData {
    std::vector<int>* Fired_Cells_mod = nullptr;
    std::vector<int>* Fired_Cells_id = nullptr;
    std::vector<double>* Fired_Cells_x = nullptr;
    std::vector<double>* Fired_Cells_y = nullptr;
    std::vector<double>* Fired_Cells_z = nullptr;
    std::vector<int>* isCellComplete = nullptr;
    std::vector<int>* Fired_Cells_adc1 = nullptr;
    std::vector<int>* Fired_Cells_adc2 = nullptr;
    std::vector<int>* Fired_Cells_tdc1 = nullptr;
    std::vector<int>* Fired_Cells_tdc2 = nullptr;
    std::vector<double>* Fired_Cell_true_Hit_x = nullptr;
    std::vector<double>* Fired_Cell_true_Hit_y = nullptr;
    std::vector<double>* Fired_Cell_true_Hit_z = nullptr;
    std::vector<double>* Fired_Cell_true_Hit_t = nullptr;
    std::vector<double>* Fired_Cell_true_Hit_e = nullptr;
    std::vector<int>* Fired_by_primary_antimu = nullptr;
    std::vector<int>* Fired_by_primary_neutron = nullptr;
    std::vector<int>* IsEarliestCell_neutron = nullptr;
    std::vector<double>* Reconstructed_HitPosition_x = nullptr;
    std::vector<double>* Reconstructed_HitPosition_y = nullptr;
    std::vector<double>* Reconstructed_HitPosition_z = nullptr;
    std::vector<double>* Reconstructed_HitTime = nullptr;
    std::vector<double>* Reconstructed_Energy = nullptr;
    std::vector<double>* ExpectedNeutron_HitPosition_x_ = nullptr;
    std::vector<double>* ExpectedNeutron_HitPosition_y_ = nullptr;
    std::vector<double>* ExpectedNeutron_HitPosition_z_ = nullptr;
    std::vector<double>* Expected_HitTime_ = nullptr;
    std::vector<double>* True_FlightLength = nullptr;
    std::vector<double>* Reconstructed_FlightLength = nullptr;
    std::vector<double>* Residuals_HitTime_ = nullptr;
    std::vector<double>* Residuals_HitSpace_ = nullptr;
};


/***
STRUCT: Cell Manager
*/

enum CELL_TYPE{
    COMPLETE_CELLS_SIGNAL_NEUTRONS = 0,
    COMPLETE_CELLS_SIGNAL_MUONS = 1,
    COMPLETE_CELLS_BKG_NEUTRONS = 2,
    CELLS_NONE,
};

struct cell {
    int mod, id;
    bool is_complete, is_fired_by_antimu, is_fired_by_neutron;
    double adc1, adc2, tdc1, tdc2;
    double x, y, z;
    double true_x, true_y, true_z, true_t, true_e;
    double reco_x, reco_y, reco_z, reco_t, reco_e;
    double pred_x, pred_y, pred_z, pred_t;
    double true_flight, reco_flight, hit_time_res, hit_space_res;
    CELL_TYPE ctype;

    // Costruttore di default
    cell()
        : mod(-999), id(-999),
          is_complete(false), is_fired_by_antimu(false), is_fired_by_neutron(false),
          adc1(-999.), adc2(-999.), tdc1(-999.), tdc2(-999.),
          x(-999.), y(-999.), z(-999.),
          true_x(-999.), true_y(-999.), true_z(-999.), true_t(-999.), true_e(-999.),
          reco_x(-999.), reco_y(-999.), reco_z(-999.), reco_t(-999.), reco_e(-999.),
          pred_x(-999.), pred_y(-999.), pred_z(-999.), pred_t(-999.),
          true_flight(-999.), reco_flight(-999.), hit_time_res(-999.), hit_space_res(-999.),
          ctype(CELLS_NONE) {};
};

CELL_TYPE DetermineCellType(bool is_signal,
                            const cell& c){
    if(is_signal)
    {
        if(c.is_fired_by_neutron)
        {
            return COMPLETE_CELLS_SIGNAL_NEUTRONS;
        }else
        {
            return COMPLETE_CELLS_SIGNAL_MUONS;
        }
    }else
    {
        if(c.is_fired_by_neutron)
        {
            return COMPLETE_CELLS_BKG_NEUTRONS;
        }
        
    }

    return CELLS_NONE;                             
}

/***
 * LEAVES: leaves of the output files
*/

int CCQEonHydrogen;
int NofFinalStateChargedParticles;
int nof_fired_wires;
int candidate_signal_event; // 
std::string* InteractionTarget = new std::string();
std::string* InteractionVolume_short = new std::string();
TLorentzVector* IncomingNeutrinoP4 = nullptr; 
TLorentzVector* FinalStateHadronicSystemTotal4Momentum = nullptr; 
double Interaction_vtxX;
double Interaction_vtxY;
double Interaction_vtxZ;
TVector3* Antimuon_p_true = nullptr;
TLorentzVector* Antimuon_reconstructed_P4 = nullptr;
TLorentzVector* Neutrino_reconstructed_P4_GeV = nullptr;
TVector3* PredictedNeutron_P3_GeV = nullptr; 
double PredictedNeutron_E_GeV;
// list of cells fired by neutron in ecal
// vector<int,ROOT::Detail::VecOps::RAdoptAllocator<int>>* is_cell_fired_by_neutron = nullptr;
std::vector<int>* is_cell_fired_by_neutron = nullptr;
std::vector<int>* is_cell_fired_by_antimu = nullptr;
std::vector<int>* IsEarliestCell_neutron = nullptr;
std::vector<int>* nof_compatible_cells = nullptr;
std::vector<int>* nof_compatible_cells2 = nullptr;
std::vector<int>* Fired_Cells_mod = nullptr;
std::vector<int>* isCellComplete = nullptr;
std::vector<int>* IsCompatible = nullptr;
std::vector<int>* Fired_Cells_adc1 = nullptr;
std::vector<int>* Fired_Cells_adc2 = nullptr;
std::vector<int>* Fired_Cells_tdc1 = nullptr;
std::vector<int>* Fired_Cells_tdc2 = nullptr;
std::vector<double>* Fired_Cell_true_Hit_x = nullptr;
std::vector<double>* Fired_Cell_true_Hit_y = nullptr;
std::vector<double>* Fired_Cell_true_Hit_z = nullptr;
std::vector<double>* Fired_Cell_true_Hit_t = nullptr;
std::vector<double>* Reconstructed_HitPosition_x = nullptr;
std::vector<double>* Reconstructed_HitPosition_y = nullptr;
std::vector<double>* Reconstructed_HitPosition_z = nullptr;
std::vector<double>* Reconstructed_HitTime = nullptr;
std::vector<double>* True_FlightLength = nullptr;
std::vector<double>* Reconstructed_FlightLength = nullptr;
std::vector<double>* Residuals_HitTime_ = nullptr;
std::vector<double>* Residuals_HitSpace_ = nullptr;
std::vector<double>* Reconstructed_Energy = nullptr;

/***
 * SCALE: mass ratio carbon_in_graphite / carbon_in_plastic
*/
const double scale_factor = 4.049;

enum STAGE{
    FIDUCIAL_VOLUME = 0,
    WIRES_CUT = 1,
    CHARGE_MULTIPLICITY = 2,
    ECAL_COINCIDENCE = 3,
    STAGE_NONE = 4
};

std::string getStageAsString(int stage) {
    switch (stage) {
        case FIDUCIAL_VOLUME: return "FIDUCIAL_VOLUME";
        case WIRES_CUT: return "WIRES_CUT";
        case CHARGE_MULTIPLICITY: return "CHARGE_MULTIPLICITY";
        case ECAL_COINCIDENCE: return "ECAL_COINCIDENCE";
        case STAGE_NONE: return "STAGE_NONE";
        default: return "UNKNOWN_STAGE";
    }
}

enum SELECTION{
    // SELECTION
    SELECTED_TRUE_POSITIVE = 0,
    FALSE_NEGATIVE = 1,
    // FALSE POSITIVES
    SELECTED_FALSE_POSITIVE_GRAPHITE  = 2,
    SELECTED_FALSE_POSITIVE_PLASTIC_C  = 3,
    SELECTED_FALSE_POSITIVE_PLASTIC_H  = 4,
    SELECTED_FALSE_POSITIVE_OTHER  = 5,
    // TRUE NEGATIVES
    TRUE_NEGATIVE_GRAPHITE  = 6,
    TRUE_NEGATIVE_PLASTIC_C  = 7,
    TRUE_NEGATIVE_PLASTIC_H  = 8,
    TRUE_NEGATIVE_OTHERS  = 9,
    // 
    SELECTED_POSITIVE_PLASTIC = 10,
    SELECTION_NONE
};


std::string getSelectionAsString(int selection) {
    switch (selection) {
        case SELECTED_TRUE_POSITIVE: return "SELECTED_TRUE_POSITIVE";
        case FALSE_NEGATIVE: return "FALSE_NEGATIVE";
        case SELECTED_FALSE_POSITIVE_GRAPHITE: return "SELECTED_FALSE_POSITIVE_GRAPHITE";
        case SELECTED_FALSE_POSITIVE_PLASTIC_C: return "SELECTED_FALSE_POSITIVE_PLASTIC_C";
        case SELECTED_FALSE_POSITIVE_PLASTIC_H: return "SELECTED_FALSE_POSITIVE_PLASTIC_H";
        case SELECTED_FALSE_POSITIVE_OTHER: return "SELECTED_FALSE_POSITIVE_OTHER";
        case TRUE_NEGATIVE_GRAPHITE: return "TRUE_NEGATIVE_GRAPHITE";
        case TRUE_NEGATIVE_PLASTIC_C: return "TRUE_NEGATIVE_PLASTIC_C";
        case TRUE_NEGATIVE_PLASTIC_H: return "TRUE_NEGATIVE_PLASTIC_H";
        case TRUE_NEGATIVE_OTHERS: return "TRUE_NEGATIVE_OTHERS";
        case SELECTED_POSITIVE_PLASTIC: return "SELECTED_POSITIVE_PLASTIC";
        case SELECTION_NONE: return "SELECTION_NONE";
        default: return "unknown";
    }
}

enum PARTICLE{
    ANTIMUON = 0,
    NEUTRON = 1,
    NEUTRINO = 2,
    PARTICLE_NONE = 3
};

std::string getTypeAsString(PARTICLE particle) {
    switch (particle) {
        case ANTIMUON: return "antimuon";
        case NEUTRON: return "neutron";
        case NEUTRINO: return "neutrino";
        case PARTICLE_NONE: return "none";
        default: return "unknown";
    }
}

enum KINETIC_VARIABLE{
    ENERGY = 0,
    KINETIC_ENERGY = 1,
    MOMENTUM = 2,
    TRANSVERSE_MOMENTUM = 3,
    LONGITUDINAL_MOMENTUM = 4,
    ANGLE_OF_EMISSION = 5,
    KINETIC_VARIABLE_NONE = 6
};

enum TARGET {
    GRAPHITE = 0,
    PLASTIC = 1,
    OTHER = 2,
    ANY_TARGET = 3,
    TARGET_NONE = 4
};

enum REACTION {
    CCQE_ON_H = 0,
    CC_ON_CARBON = 1,
    CCRES_ON_H = 2,
    REACTION_ANY = 3,
    REACTION_NONE = 4
};

std::string getReactionAsString(int reaction)
{
    switch (reaction) {
        case CCQE_ON_H: return "CCQE_ON_H";
        case CC_ON_CARBON: return "CC_ON_CARBON";
        case CCRES_ON_H: return "CCRES_ON_H";
        case REACTION_NONE: return "REACTION_NONE";
        default: return "unknown";
    }
}


SELECTION DetermineSelection(bool is_event_selected, bool is_signal, bool is_in_graphite, bool is_in_plastic, bool is_on_C, bool is_on_H){
    if((is_event_selected))
    // selected CCQE like on H
    {
        if(is_signal) return SELECTED_TRUE_POSITIVE;
        if(is_in_graphite) return SELECTED_FALSE_POSITIVE_GRAPHITE;
        if(is_in_plastic)
        {
            if(is_on_C) return SELECTED_FALSE_POSITIVE_PLASTIC_C;
            if(is_on_H) return SELECTED_FALSE_POSITIVE_PLASTIC_H;
        }
        return SELECTED_FALSE_POSITIVE_OTHER;
    }

    // non selected
    if(is_signal) return FALSE_NEGATIVE;
    if(is_in_graphite) return TRUE_NEGATIVE_GRAPHITE;
    if(is_in_plastic)
    {
        if(is_on_C) return TRUE_NEGATIVE_PLASTIC_C;
        if(is_on_H) return TRUE_NEGATIVE_PLASTIC_H;
    }
    return TRUE_NEGATIVE_OTHERS;
}

void plot_histograms2(TH2D* h2, std::string canvas_name){
    TCanvas* c_h2 = new TCanvas(canvas_name.c_str(), "", 900, 700);
    h2->Draw("COLZ");
    c_h2->Draw();
}

/***
RELATIONS: TRUE KINEMATIC VARIABLES, TH2D to relate quantities
 */

TH2D* true_nu_E_vs_true_n_kin_E[STAGE_NONE] = {
    new TH2D("true_nu_E_vs_true_n_kin_E_FIDUCIAL_VOLUME", "", 20u, 0., 6., 20u, 0., 1.),
    new TH2D("true_nu_E_vs_true_n_kin_E_WIRES_CUT", "", 20u, 0., 6., 20u, 0., 1.),
    new TH2D("true_nu_E_vs_true_n_kin_E_CHARGE_MULTIPLICITY", "", 20u, 0., 6., 20u, 0., 1.),
    new TH2D("true_nu_E_vs_true_n_kin_E_ECAL_COINCIDENCE", "", 20u, 0., 6., 20u, 0., 1.)
};

TH2D* true_nu_E_vs_true_mu_E[STAGE_NONE] = {
    new TH2D("true_nu_E_vs_true_mu_E_FIDUCIAL_VOLUME", "", 20u, 0., 6., 20u, 0., 6000.),
    new TH2D("true_nu_E_vs_true_mu_E_WIRES_CUT", "", 20u, 0., 6., 20u, 0., 6000.),
    new TH2D("true_nu_E_vs_true_mu_E_CHARGE_MULTIPLICITY", "", 20u, 0., 6., 20u, 0., 6000.),
    new TH2D("true_nu_E_vs_true_mu_E_ECAL_COINCIDENCE", "", 20u, 0., 6., 20u, 0., 6000.)
};

TH2D* true_mu_E_vs_true_n_kin_E[STAGE_NONE] = {
    new TH2D("true_mu_E_vs_true_n_kin_E_FIDUCIAL_VOLUME", "", 20u, 0., 6000., 20u, 0., 6.),
    new TH2D("true_mu_E_vs_true_n_kin_E_WIRES_CUT", "", 20u, 0., 6000., 20u, 0., 6.),
    new TH2D("true_mu_E_vs_true_n_kin_E_CHARGE_MULTIPLICITY", "", 20u, 0., 6000., 20u, 0., 6.),
    new TH2D("true_mu_E_vs_true_n_kin_E_ECAL_COINCIDENCE", "", 20u, 0., 6000., 20u, 0., 6.)
};

TH2D* neutron_space_time_residulas = new TH2D("neutron_space_time_residulas", ";neutron hit predicted - reco [mm]; neutron hit predicted - reco [ns]", 250u, 0, 500., 1000u, -5, 5.);

bool IsInsideSpot(const double vtxX, const double vtxY, const double vtxZ, double radius_m = 1.){
    TVector3 spot_center = {0., -2.5, 23.5};
    TVector3 vtx = {vtxX, vtxY, vtxZ};
    return ((spot_center - vtx).Mag() < radius_m);
}

struct h_bins
{
    uint nof_bins = nof_bin_nu_E;
    double low = 0.;
    double up = 6.;
    h_bins(uint bins, double l, double u) : nof_bins(bins), low(l), up(u) {};
};

h_bins Get_h_bins(PARTICLE p)
{
    if(p == ANTIMUON)
    {
        return {nof_bin_nu_E, 0., 6000.};
    }else if(p == NEUTRINO)
    {
        return {nof_bin_nu_E, 0., 6.};
    }else if(p == NEUTRON){
        return {20, 0., 1.};
    }else
    {
        return {nof_bin_nu_E, 0., 6.};
    }
};

/***
 CLASS: Hist Manager
 */
class Hist_Manager{
    private:
        PARTICLE particle;

        TH1D* energy_spectrum[STAGE_NONE+1][SELECTION_NONE+1][REACTION_NONE+1][2];

        void InitArrayHist() {
            // Get the name of the particle for histogram naming
            std::string particleName = getTypeAsString(particle);

            auto bins = Get_h_bins(particle);

            // Loop through all stages
            for (int stage = 0; stage < STAGE_NONE + 1; ++stage) {
                std::string stageName = getStageAsString(stage);

                // Loop through all selections
                for (int selection = 0; selection < SELECTION_NONE + 1; ++selection) {
                    std::string selectionName = getSelectionAsString(selection);

                    // Loop through all reactions
                    for (int reaction = 0; reaction < REACTION_NONE + 1; ++reaction) {
                        std::string reactionName = getReactionAsString(reaction);

                        // Generate histogram names
                        std::string hNameTrue = particleName + "_true_" + stageName + "_" + selectionName + "_" + reactionName;
                        std::string hNameReco = particleName + "_reco_" + stageName + "_" + selectionName + "_" + reactionName;

                        // Initialize histograms for true and reconstructed energy spectra
                        energy_spectrum[stage][selection][reaction][0] = new TH1D(hNameTrue.c_str(), "", bins.nof_bins, bins.low, bins.up);
                        energy_spectrum[stage][selection][reaction][1] = new TH1D(hNameReco.c_str(), "", bins.nof_bins, bins.low, bins.up);
                    }
                }
            }
        }

    public:
        // Constructor
        // explicit avoid unwanted conversions of type HistManager (if any possible)
        explicit Hist_Manager(PARTICLE p) : particle(p) {
            InitArrayHist(); // Call InitHist in the constructor
        }

        std::string GetParticleName(){
            return getTypeAsString(particle);
        }

        TH1D* GetHistogram(STAGE stage, SELECTION selection, REACTION reaction, int i) const {
            if (stage >= 0 && stage <= STAGE_NONE &&
                selection >= 0 && selection <= SELECTION_NONE &&
                reaction >= 0 && reaction <= REACTION_NONE &&
                (i == 0 || i == 1)) {
                return energy_spectrum[stage][selection][reaction][i];
            } else {
                throw std::out_of_range("Invalid indices in Hist_Manager::GetHistogram");
            }
    }

        void Fill(STAGE stage, SELECTION selection, REACTION reaction, int i, double E)
        {
            // Check for valid indices to avoid out-of-bounds errors
            if (stage >= 0 && stage <= STAGE_NONE &&
                selection >= 0 && selection <= SELECTION_NONE &&
                reaction >= 0 && reaction <= REACTION_NONE &&
                (i == 0 || i == 1)) { // i must be 0 (true) or 1 (reco)

                // Fill the corresponding histogram
                energy_spectrum[stage][selection][reaction][i]->Fill(E);
            } else {
                throw std::out_of_range("Invalid indices in Hist_Manager::Fill");
            }
        }

        void CompareStages(std::string canvas_name)
        {
            TCanvas *c1 = new TCanvas(canvas_name.c_str(), canvas_name.c_str(), 900, 700);
            
            TPad* pad1 = new TPad("pad1", "Main Plot", 0.0, 0.3, 1.0, 1.0); // Top pad
            TPad* pad2 = new TPad("pad2", "Ratio Plot", 0.0, 0.0, 1.0, 0.3); // Bottom pad

            pad1->SetBottomMargin(0.02); // Reduce bottom margin for top pad
            pad2->SetTopMargin(0.02);   // Reduce top margin for bottom pad
            pad2->SetBottomMargin(0.3); // Increase bottom margin for labels
            pad1->Draw();
            pad2->Draw();
            
            pad1->cd();

            auto h1 = energy_spectrum[FIDUCIAL_VOLUME][SELECTION_NONE][CCQE_ON_H][0];
            auto h2 = energy_spectrum[WIRES_CUT][SELECTION_NONE][CCQE_ON_H][0];
            // auto h3 = energy_spectrum[CHARGE_MULTIPLICITY][SELECTION_NONE][CCQE_ON_H][0];
            auto h4 = energy_spectrum[ECAL_COINCIDENCE][SELECTION_NONE][CCQE_ON_H][0];
        
            h1 -> SetLineColor(kRed);
            h2 -> SetLineColor(kBlue);
            // h3 -> SetLineColor(kBlue);
            h4 -> SetLineColor(kBlack);
        
            h1 -> Draw("E HIST");
            h2 -> Draw("E HIST SAME");
            // h3 -> Draw("E HIST SAME");
            h4 -> Draw("E HIST SAME");

            TLegend* legend = new TLegend(0.5, 0.5, 0.9, 0.9);
            std::string particle_name = GetParticleName();
            legend->AddEntry(h1, TString::Format("%s (STAGE FIDUCIAL)", particle_name.c_str()).Data(), "l");
            legend->AddEntry(h2, TString::Format("%s (STAGE WIRES CUT)", particle_name.c_str()).Data(), "l");
            // legend->AddEntry(h3, TString::Format("%s (STAGE CHARGE MULTI)", particle_name.c_str()).Data(), "l");
            legend->AddEntry(h4, TString::Format("%s (STAGE ECAL COINC.)", particle_name.c_str()).Data(), "l");
            
            legend->Draw();
            
            // Calculate the ratio
            pad2->cd();
            h1->Sumw2();
            h2->Sumw2();
            // h3->Sumw2();
            h4->Sumw2();
            TH1D* ratio21 = (TH1D*)h2->Clone("ratio21");
            // TH1D* ratio32 = (TH1D*)h3->Clone("ratio32");
            TH1D* ratio43 = (TH1D*)h4->Clone("ratio43");
            
            ratio21->Divide(h1);
            // ratio32->Divide(h2);
            ratio43->Divide(h2);
            
            ratio21->SetTitle(""); // Remove title for ratio plot
            
            ratio21->SetLineColor(kBlue);
            // ratio32->SetLineColor(kBlue);
            ratio43->SetLineColor(kBlack);
            
            ratio21->GetYaxis()->SetTitle("Cut Efficiency");
            std::string x_label = (particle == ANTIMUON) ? "True Energy [MeV]" : "True Energy [GeV]";
            ratio21->GetXaxis()->SetTitle(x_label.c_str());
            ratio21->GetYaxis()->SetNdivisions(505);
            ratio21->GetYaxis()->SetTitleSize(0.1);
            ratio21->GetYaxis()->SetTitleOffset(0.4);
            ratio21->GetYaxis()->SetLabelSize(0.08);
            ratio21->GetXaxis()->SetTitleSize(0.1);
            ratio21->GetXaxis()->SetTitleOffset(0.8);
            ratio21->GetXaxis()->SetLabelSize(0.08);
            
            ratio21->Draw("E1");
            // ratio32->Draw("E1 SAME");
            ratio43->Draw("E1 SAME");
        }

};

Hist_Manager positive_mu_hist(ANTIMUON);
Hist_Manager antinu_hist(NEUTRINO);
Hist_Manager neutron_hist(NEUTRON);
Hist_Manager Q2(PARTICLE_NONE);


class Cells_Manager {
private:
    
    std::vector<cell> incomplete_cells; 
    std::vector<cell> complete_cells;
    std::vector<cell> fired_by_signal_neutrons;
    std::vector<cell> fired_by_bkg_neutrons;
    std::vector<cell> fired_by_signal_muons;
    std::vector<cell> compatible_cells;
    
    double best_time_residual;
    double best_space_residual; 

public:

    void add_cell(const cell& c) {

        switch (c.ctype) {
            case COMPLETE_CELLS_SIGNAL_NEUTRONS:
            
                complete_cells.push_back(c);
                fired_by_signal_neutrons.push_back(c);
                break;
            
            case COMPLETE_CELLS_BKG_NEUTRONS:
            
                complete_cells.push_back(c);
                fired_by_bkg_neutrons.push_back(c);
                break;
            
            case COMPLETE_CELLS_SIGNAL_MUONS:
                complete_cells.push_back(c);
                fired_by_signal_muons.push_back(c);
                break;
            
            default:
                break;
        }
    }

    // Ottieni il numero di celle per ogni tipo
    size_t get_cell_count(CELL_TYPE type) const {
        switch (type) {
            
            case COMPLETE_CELLS_SIGNAL_NEUTRONS:
                return fired_by_signal_neutrons.size();
            
            case COMPLETE_CELLS_BKG_NEUTRONS:
                return fired_by_bkg_neutrons.size();
            
            case COMPLETE_CELLS_SIGNAL_MUONS:
                return fired_by_signal_muons.size();
            
            default:
                return 0; // Tipo non valido
        }
    }

    // Resetta tutte le liste
    void clear_all_cells() {
        fired_by_signal_neutrons.clear();
        fired_by_bkg_neutrons.clear();
        fired_by_signal_muons.clear();
    }

    void init_event_cells(const CellData& cell_data, bool is_signal_event) 
    {
        size_t n_cells = cell_data.isCellComplete->size();
        for (size_t i = 0; i < n_cells; ++i) 
        {
            cell this_cell;
            this_cell.mod = cell_data.Fired_Cells_mod->at(i);
            this_cell.id = cell_data.Fired_Cells_id->at(i);
            this_cell.is_complete = cell_data.isCellComplete->at(i);
            this_cell.is_fired_by_antimu = cell_data.Fired_by_primary_antimu->at(i);
            this_cell.is_fired_by_neutron = cell_data.Fired_by_primary_neutron->at(i);
            this_cell.tdc1 = cell_data.Fired_Cells_tdc1->at(i);
            this_cell.tdc2 = cell_data.Fired_Cells_tdc2->at(i);
            this_cell.adc1 = cell_data.Fired_Cells_adc1->at(i);
            this_cell.adc2 = cell_data.Fired_Cells_adc2->at(i);
            this_cell.x = cell_data.Fired_Cells_x->at(i);
            this_cell.y = cell_data.Fired_Cells_y->at(i);
            this_cell.z = cell_data.Fired_Cells_z->at(i);
            this_cell.true_x = cell_data.Fired_Cell_true_Hit_x->at(i);
            this_cell.true_y = cell_data.Fired_Cell_true_Hit_y->at(i);
            this_cell.true_z = cell_data.Fired_Cell_true_Hit_z->at(i);
            this_cell.true_t = cell_data.Fired_Cell_true_Hit_t->at(i);
            this_cell.true_e = cell_data.Fired_Cell_true_Hit_e->at(i);
            this_cell.reco_x = cell_data.Reconstructed_HitPosition_x->at(i);
            this_cell.reco_y = cell_data.Reconstructed_HitPosition_y->at(i);
            this_cell.reco_z = cell_data.Reconstructed_HitPosition_z->at(i);
            this_cell.reco_t = cell_data.Reconstructed_HitTime->at(i);
            this_cell.reco_e = cell_data.Reconstructed_Energy->at(i);
            this_cell.pred_x = cell_data.ExpectedNeutron_HitPosition_x_->at(i);
            this_cell.pred_y = cell_data.ExpectedNeutron_HitPosition_y_->at(i);
            this_cell.pred_z = cell_data.ExpectedNeutron_HitPosition_z_->at(i);
            this_cell.pred_t = cell_data.Expected_HitTime_->at(i);
            this_cell.true_flight = cell_data.True_FlightLength->at(i);
            this_cell.reco_flight = cell_data.Reconstructed_FlightLength->at(i);
            this_cell.hit_time_res = cell_data.Residuals_HitTime_->at(i);
            this_cell.hit_space_res = cell_data.Residuals_HitSpace_->at(i);

            this_cell.ctype = DetermineCellType(is_signal_event, this_cell);

            add_cell(this_cell);
        }
    }



};


TH1D* flight_length_neutron_signal = new TH1D("flight_length_neutron_signal", ";reconstructed flight length [mm];count", 300, 0, 4000);
TH1D* flight_length_neutron_bkg = new TH1D("flight_length_neutron_bkg", ";reconstructed flight length [mm];count", 300, 0, 4000);
TH2D* flight_length_neutron_signal_vs_res_time = new TH2D("fl_neutron_signal_vs_res1", ";reconstructed flight length [mm];predicted time - reco time [ns]", 300, 0, 4000, 100, -4, 4 );
TH2D* flight_length_neutron_bkg_vs_res_time = new TH2D("fl_neutron_signal_vs_res2", ";reconstructed flight length [mm];predicted time - reco time [ns]", 300, 0, 4000, 100, -4, 4 );


TH2D* flight_length_neutron_signal_vs_res_space = new TH2D("fl_neutron_signal_vs_res11", ";reconstructed flight length [mm];space_residuals [mm]", 300, 0, 4000, 250, 0, 1000 );
TH2D* flight_length_neutron_bkg_vs_res_space = new TH2D("fl_neutron_signal_vs_res22", ";reconstructed flight length [mm];space_residuals [mm]", 300, 0, 4000, 250, 0, 1000);

void LoopOverCells1D(const CellData& cell_data, TH1D* hist) {
    for (size_t i = 0; i < cell_data.Fired_Cells_mod->size(); i++) {
        if (cell_data.isCellComplete->at(i) && cell_data.Fired_by_primary_neutron->at(i) == 1) 
        {
            // if(cell_data.IsEarliestCell_neutron -> at(i))
            if(1)
            {
                hist->Fill(cell_data.Reconstructed_FlightLength -> at(i));
            }
        }
    }
}

void LoopOverCells2D(const CellData& cell_data, TH2D* hist1, TH2* hist2) {
    for (size_t i = 0; i < cell_data.Fired_Cells_mod->size(); i++) {
        if (cell_data.isCellComplete->at(i) && cell_data.Fired_by_primary_neutron->at(i) == 1) 
        {
            // if(cell_data.IsEarliestCell_neutron -> at(i))
            if(fabs(cell_data.Residuals_HitSpace_ -> at(i)) < 150.)
            {
                hist1->Fill(cell_data.Reconstructed_FlightLength -> at(i), cell_data.Residuals_HitTime_ -> at(i));
            }else if(fabs(cell_data.Residuals_HitTime_ -> at(i) / cell_data.Reconstructed_HitTime-> at(i)) < 0.05)
            {
                hist2->Fill(cell_data.Reconstructed_FlightLength -> at(i), cell_data.Residuals_HitSpace_ -> at(i));
            }
        }
    }
}

// void LoopOverCells2D(const CellData& cell_data, TH2D* hist1, TH2* hist2) {
//     for (size_t i = 0; i < cell_data.Fired_Cells_mod->size(); i++) {
//         if (cell_data.isCellComplete->at(i) && cell_data.Fired_by_primary_neutron->at(i) == 1) 
//         {
//             // if(cell_data.IsEarliestCell_neutron -> at(i))
//             if(fabs(cell_data.Residuals_HitSpace_ -> at(i)) < 150. && cell_data.IsEarliestCell_neutron -> at(i))
//             {
//                 hist1->Fill(cell_data.Reconstructed_FlightLength -> at(i), cell_data.Residuals_HitTime_ -> at(i));
//             }else if(fabs(cell_data.Residuals_HitTime_ -> at(i)) < 1. && cell_data.IsEarliestCell_neutron -> at(i))
//             {
//                 hist2->Fill(cell_data.Reconstructed_FlightLength -> at(i), cell_data.Residuals_HitSpace_ -> at(i));
//             }
//         }
//     }
// }

struct Counter
{
    int ccqe_on_H[4] = {0, 0, 0, 0};
    int cc_on_C_graphite[4] = {0, 0, 0, 0};
    int cc_on_C_plastic[4] = {0, 0, 0, 0};
    int ccres_on_H[4] = {0, 0, 0, 0};
    int cc_on_others[4] = {0, 0, 0, 0};
};


void PrintCounter(const Counter& counter) {
    // Larghezza delle colonne
    const int col_width = 25;

    // Intestazione della tabella
    std::cout << std::left
              << std::setw(col_width) << "Category"
              << std::setw(col_width) << "Entry 0"
              << std::setw(col_width) << "Entry 1"
              << std::setw(col_width) << "Entry 2"
              << std::setw(col_width) << "Entry 3"
              << std::endl;

    std::cout << std::string(5 * col_width, '-') << std::endl;

    // Righe della tabella
    std::cout << std::setw(col_width) << "CCQE on H"
              << std::setw(col_width) << counter.ccqe_on_H[0]
              << std::setw(col_width) << counter.ccqe_on_H[1]
              << std::setw(col_width) << counter.ccqe_on_H[2]
              << std::setw(col_width) << counter.ccqe_on_H[3]
              << std::endl;

    std::cout << std::setw(col_width) << "CC on C (Graphite)"
              << std::setw(col_width) << counter.cc_on_C_graphite[0]
              << std::setw(col_width) << counter.cc_on_C_graphite[1]
              << std::setw(col_width) << counter.cc_on_C_graphite[2]
              << std::setw(col_width) << counter.cc_on_C_graphite[3]
              << std::endl;

    std::cout << std::setw(col_width) << "CC on C (Plastic)"
              << std::setw(col_width) << counter.cc_on_C_plastic[0]
              << std::setw(col_width) << counter.cc_on_C_plastic[1]
              << std::setw(col_width) << counter.cc_on_C_plastic[2]
              << std::setw(col_width) << counter.cc_on_C_plastic[3]
              << std::endl;

    std::cout << std::setw(col_width) << "CCRes on H"
              << std::setw(col_width) << counter.ccres_on_H[0]
              << std::setw(col_width) << counter.ccres_on_H[1]
              << std::setw(col_width) << counter.ccres_on_H[2]
              << std::setw(col_width) << counter.ccres_on_H[3]
              << std::endl;

    std::cout << std::setw(col_width) << "CC on Others"
              << std::setw(col_width) << counter.cc_on_others[0]
              << std::setw(col_width) << counter.cc_on_others[1]
              << std::setw(col_width) << counter.cc_on_others[2]
              << std::setw(col_width) << counter.cc_on_others[3]
              << std::endl;
}

void Increment_counter(Counter& counter, STAGE stage, 
                        bool is_signal, bool is_on_C, bool is_in_graphite, bool is_in_plastic, bool is_on_H)
{
    if (is_signal) {
        counter.ccqe_on_H[stage]++;
    } else if (is_on_C && is_in_graphite) {
        counter.cc_on_C_graphite[stage]++;
    } else if (is_on_C && is_in_plastic) {
        counter.cc_on_C_plastic[stage]++;
    } else if (is_on_H && !is_signal) {
        counter.ccres_on_H[stage]++;
    }else{
        counter.cc_on_others[stage]++;
    }
}

void Fill_Final_Histos(SELECTION selection, 
                      double true_nu_E, double reco_nu_E, 
                      const CellData& cell_data,
                      Counter& counter,
                      int event_number
                    //   double true_mu_E, double reco_mu_E,
                    //   double true_n_K, double reco_n_K,
                      ){
    switch (selection)
    {
        case SELECTED_TRUE_POSITIVE:
            // efficiency & matrix unfold_____________________________
            if(event_number%2){
                response_matrix[SIGNAL_SAMPLE2]->Fill(reco_nu_E, true_nu_E);
                efficiency[SIGNAL_SAMPLE2] -> Fill(true, true_nu_E);
                hist4unfold[SIGNAL_SAMPLE2][0] -> Fill(true_nu_E); 
                hist4unfold[SIGNAL_SAMPLE2][1] -> Fill(reco_nu_E); 

                response_matrix[SIGNAL_BKG_SAMPLE2]->Fill(reco_nu_E, true_nu_E);
                efficiency[SIGNAL_BKG_SAMPLE2] -> Fill(true, true_nu_E);
                hist4unfold[SIGNAL_BKG_SAMPLE2][0] -> Fill(true_nu_E); 
                hist4unfold[SIGNAL_BKG_SAMPLE2][1] -> Fill(reco_nu_E); 

            }else{
                hist4unfold[SIGNAL_SAMPLE3][0] -> Fill(true_nu_E); 
                hist4unfold[SIGNAL_SAMPLE3][1] -> Fill(reco_nu_E); 
                hist4unfold[SIGNAL_BKG_SAMPLE3][0] -> Fill(true_nu_E); 
                hist4unfold[SIGNAL_BKG_SAMPLE3][1] -> Fill(reco_nu_E); 
            }
            response_matrix[SIGNAL_SAMPLE1]->Fill(reco_nu_E, true_nu_E);
            efficiency[SIGNAL_SAMPLE1] -> Fill(true, true_nu_E);
            hist4unfold[SIGNAL_SAMPLE1][0] -> Fill(true_nu_E); 
            hist4unfold[SIGNAL_SAMPLE1][1] -> Fill(reco_nu_E);

            response_matrix[SIGNAL_BKG_SAMPLE1]->Fill(reco_nu_E, true_nu_E);
            efficiency[SIGNAL_BKG_SAMPLE1] -> Fill(true, true_nu_E);
            hist4unfold[SIGNAL_BKG_SAMPLE1][0] -> Fill(true_nu_E); 
            hist4unfold[SIGNAL_BKG_SAMPLE1][1] -> Fill(reco_nu_E);

            // _________________________________________________

            counter.ccqe_on_H[FIDUCIAL_VOLUME]++;
            counter.ccqe_on_H[ECAL_COINCIDENCE]++;

            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_POSITIVE_PLASTIC, REACTION_ANY, 0, true_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_POSITIVE_PLASTIC, REACTION_ANY, 1, reco_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_TRUE_POSITIVE, CCQE_ON_H, 0, true_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_TRUE_POSITIVE, CCQE_ON_H, 1, reco_nu_E);

            LoopOverCells1D(cell_data, flight_length_neutron_signal);
            LoopOverCells2D(cell_data, flight_length_neutron_signal_vs_res_time, flight_length_neutron_signal_vs_res_space);
            
            break;

        case FALSE_NEGATIVE:
            // efficiency & matrix unfold_____________________________
            if(event_number%2){
                response_matrix[SIGNAL_SAMPLE2]->Miss(true_nu_E);
                efficiency[SIGNAL_SAMPLE2] -> Fill(false, true_nu_E);
                hist4unfold[SIGNAL_SAMPLE2][0] -> Fill(true_nu_E); 

                response_matrix[SIGNAL_BKG_SAMPLE2]->Miss(true_nu_E);
                efficiency[SIGNAL_BKG_SAMPLE2] -> Fill(false, true_nu_E);
                hist4unfold[SIGNAL_BKG_SAMPLE2][0] -> Fill(true_nu_E); 
            } else {
                hist4unfold[SIGNAL_SAMPLE3][0] -> Fill(true_nu_E); 
                hist4unfold[SIGNAL_BKG_SAMPLE3][0] -> Fill(true_nu_E); 
            }
            response_matrix[SIGNAL_SAMPLE1]->Miss(true_nu_E);
            efficiency[SIGNAL_SAMPLE1] -> Fill(false, true_nu_E);
            hist4unfold[SIGNAL_SAMPLE1][0] -> Fill(true_nu_E); 
            
            response_matrix[SIGNAL_BKG_SAMPLE1]->Miss(true_nu_E);
            efficiency[SIGNAL_BKG_SAMPLE1] -> Fill(false, true_nu_E);
            hist4unfold[SIGNAL_BKG_SAMPLE1][0] -> Fill(true_nu_E); 
            
            // _________________________________________________

            counter.ccqe_on_H[FIDUCIAL_VOLUME]++;

            LoopOverCells1D(cell_data, flight_length_neutron_signal);
            LoopOverCells2D(cell_data, flight_length_neutron_signal_vs_res_time, flight_length_neutron_signal_vs_res_space);

            break;

        case SELECTED_FALSE_POSITIVE_GRAPHITE:
            // efficiency & matrix unfold_____________________________
            if(event_number%2){
                response_matrix[BKG_SAMPLE2]->Fill(reco_nu_E, true_nu_E);
                efficiency[BKG_SAMPLE2] -> Fill(true, true_nu_E);
                hist4unfold[BKG_SAMPLE2][0] -> Fill(true_nu_E); 
                hist4unfold[BKG_SAMPLE2][1] -> Fill(reco_nu_E); 
            } else {
                hist4unfold[BKG_SAMPLE3][0] -> Fill(true_nu_E); 
                hist4unfold[BKG_SAMPLE3][1] -> Fill(reco_nu_E); 
            }
            response_matrix[BKG_SAMPLE1]->Fill(reco_nu_E, true_nu_E);
            efficiency[BKG_SAMPLE1] -> Fill(true, true_nu_E);
            hist4unfold[BKG_SAMPLE1][0] -> Fill(true_nu_E); 
            hist4unfold[BKG_SAMPLE1][1] -> Fill(reco_nu_E); 
            
            // _________________________________________________
            
            counter.cc_on_C_graphite[FIDUCIAL_VOLUME]++;
            counter.cc_on_C_graphite[ECAL_COINCIDENCE]++;

            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_GRAPHITE, CC_ON_CARBON, 0, true_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_GRAPHITE, CC_ON_CARBON, 1, reco_nu_E);

            LoopOverCells1D(cell_data, flight_length_neutron_bkg);
            LoopOverCells2D(cell_data, flight_length_neutron_bkg_vs_res_time, flight_length_neutron_bkg_vs_res_space);

            break;

        case SELECTED_FALSE_POSITIVE_PLASTIC_C:
            // efficiency & matrix unfold_____________________________
            if(event_number%2){
                response_matrix[SIGNAL_BKG_SAMPLE2]->Fill(reco_nu_E, true_nu_E);
                efficiency[SIGNAL_BKG_SAMPLE2] -> Fill(true, true_nu_E);
                hist4unfold[SIGNAL_BKG_SAMPLE2][0] -> Fill(true_nu_E);
                hist4unfold[SIGNAL_BKG_SAMPLE2][1] -> Fill(reco_nu_E);
            } else {
                hist4unfold[SIGNAL_BKG_SAMPLE3][0] -> Fill(true_nu_E);
                hist4unfold[SIGNAL_BKG_SAMPLE3][1] -> Fill(reco_nu_E);
            }
            response_matrix[SIGNAL_BKG_SAMPLE1]->Fill(reco_nu_E, true_nu_E);
            efficiency[SIGNAL_BKG_SAMPLE1] -> Fill(true, true_nu_E);
            hist4unfold[SIGNAL_BKG_SAMPLE1][0] -> Fill(true_nu_E);
            hist4unfold[SIGNAL_BKG_SAMPLE1][1] -> Fill(reco_nu_E);
            // _________________________________________________

            counter.cc_on_C_plastic[FIDUCIAL_VOLUME]++;
            counter.cc_on_C_plastic[ECAL_COINCIDENCE]++;

            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_POSITIVE_PLASTIC, REACTION_ANY, 0, true_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_POSITIVE_PLASTIC, REACTION_ANY, 1, reco_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_PLASTIC_C, CC_ON_CARBON, 0, true_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_PLASTIC_C, CC_ON_CARBON, 1, reco_nu_E);
            
            LoopOverCells1D(cell_data, flight_length_neutron_bkg);
            LoopOverCells2D(cell_data, flight_length_neutron_bkg_vs_res_time, flight_length_neutron_bkg_vs_res_space);

            break;

        case SELECTED_FALSE_POSITIVE_PLASTIC_H:
            
            counter.ccres_on_H[FIDUCIAL_VOLUME]++;
            counter.ccres_on_H[ECAL_COINCIDENCE]++;

            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_POSITIVE_PLASTIC, REACTION_ANY, 0, true_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_POSITIVE_PLASTIC, REACTION_ANY, 1, reco_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_PLASTIC_H, CCRES_ON_H, 0, true_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_PLASTIC_H, CCRES_ON_H, 1, reco_nu_E);

            break;

        case SELECTED_FALSE_POSITIVE_OTHER:
            counter.cc_on_others[FIDUCIAL_VOLUME]++;
            break;

        // TRUE NEGATIVES
        case TRUE_NEGATIVE_GRAPHITE:
            // efficiency & matrix unfold_____________________________
            if(event_number%2){
                response_matrix[BKG_SAMPLE2]->Miss(true_nu_E);
                efficiency[BKG_SAMPLE2] -> Fill(false, true_nu_E);
                hist4unfold[BKG_SAMPLE2][0] -> Fill(true_nu_E); 
            } else {
                hist4unfold[BKG_SAMPLE3][0] -> Fill(true_nu_E); 
            }
            response_matrix[BKG_SAMPLE1]->Miss(true_nu_E);
            efficiency[BKG_SAMPLE1] -> Fill(false, true_nu_E);
            hist4unfold[BKG_SAMPLE1][0] -> Fill(true_nu_E); 
            // _________________________________________________
            
            counter.cc_on_C_graphite[FIDUCIAL_VOLUME]++;
            
            LoopOverCells1D(cell_data, flight_length_neutron_bkg);
            LoopOverCells2D(cell_data, flight_length_neutron_bkg_vs_res_time, flight_length_neutron_bkg_vs_res_space);

            break;

        case TRUE_NEGATIVE_PLASTIC_C:
            // efficiency & matrix unfold_____________________________
            if(event_number%2){
                response_matrix[SIGNAL_BKG_SAMPLE2]->Miss(true_nu_E);
                efficiency[SIGNAL_BKG_SAMPLE2] -> Fill(false, true_nu_E);
                hist4unfold[SIGNAL_BKG_SAMPLE2][0] -> Fill(true_nu_E); 
            } else {
                hist4unfold[SIGNAL_BKG_SAMPLE3][0] -> Fill(true_nu_E); 
            }
            response_matrix[SIGNAL_BKG_SAMPLE1]->Miss(true_nu_E);
            efficiency[SIGNAL_BKG_SAMPLE1] -> Fill(false, true_nu_E);
            hist4unfold[SIGNAL_BKG_SAMPLE1][0] -> Fill(true_nu_E); 
           
            // _________________________________________________
            
            counter.cc_on_C_plastic[FIDUCIAL_VOLUME]++;
            
            LoopOverCells1D(cell_data, flight_length_neutron_bkg);
            LoopOverCells2D(cell_data, flight_length_neutron_bkg_vs_res_time, flight_length_neutron_bkg_vs_res_space);

            break;

        case TRUE_NEGATIVE_PLASTIC_H:
            counter.ccres_on_H[FIDUCIAL_VOLUME]++;
            break;

        case TRUE_NEGATIVE_OTHERS:
            counter.cc_on_others[FIDUCIAL_VOLUME]++;
            break;

        case SELECTION_NONE:
            break;
    }
}

std::vector<int> are_cells_compatibles(
                                    const std::vector<double>& flight_length,
                                    const std::vector<double>& tof,
                                     const std::vector<double>& space_residuals,
                                      const std::vector<double>& time_residuals,
                                      const std::vector<double>& reco_energy,
                                      const std::vector<int>& is_earliest
                                      ){
    /*
      Input: time/space residuals between the expected and the reconstructed hit in the cell
      and look for those cell compatible in the space and in the time
    */
    std::vector<int> isCompatible(space_residuals.size());
    
    for (size_t i = 0; i < space_residuals.size(); i++){
        // if(!is_earliest[i]) continue;    
        if((space_residuals[i] == -999.)){
            isCompatible[i] = 0;
        }else{
            float expected_sigma_time = 0.05 / sqrt(reco_energy[i]*1e-3);
            float expected_sigma_space = 10. / sqrt(reco_energy[i]*1e-3);
            bool space_pass = (fabs(space_residuals[i]) <= 150);
            bool time_pass = (fabs(time_residuals[i])) < 3*expected_sigma_time;
            // bool space_pass = (fabs(space_residuals[i])) < 3*expected_sigma_space;
            if(space_pass && time_pass)
            {
                isCompatible[i] = 1;
            }else{
                isCompatible[i] = 0;
            }
        }
    }
    return isCompatible;
}

bool has_compatible_cell(double best_time_residual,
                               double best_space_residual,
                               int is_cell_complete,
                               double reco_deposited_energy){
    
    float expected_time_res = 0.06 / sqrt(reco_deposited_energy * 1e-3); // ns
    float expected_space_res = 10; // mm
    
    bool time_pass = (fabs(best_time_residual) <= 3. * expected_time_res);
    bool space_pass = (fabs(best_space_residual) <= expected_space_res);
    
    if(is_cell_complete){
        // if(time_pass | space_pass){
        if(time_pass && space_pass){
            return true;
        }else{
            return false;
        }
    }else{
        return false;
    }
}

/***
 * EFFICIENCIES:
 */
TEfficiency* eff_neutron_vs_reco_neutrino_E = new TEfficiency("eff",";reco #bar{#nu}_{#mu} E [GeV];#epsilon", 24, 0, 6);
TH1D* h_pass = new TH1D("pass", "pass", 24, 0, 6);
TH1D* h_tot = new TH1D("tot", "tot", 24, 0, 6);
TEfficiency* ECAL_neutron_efficiency_vs_neutron_K = new TEfficiency("eff1",";true neutron K [GeV];#epsilon", 20, 0, 1);
TEfficiency* neutrons_w_hit_in_ECAL = new TEfficiency("eff2",";true neutron K [GeV];#epsilon", 20, 0, 1);
TEfficiency* neutrons_w_complete_cells_in_ECAL = new TEfficiency("eff3",";true neutron K [GeV];#epsilon", 20, 0, 1);
TEfficiency* neutrons_w_compatible_cells_in_ECAL = new TEfficiency("eff4",";true neutron K [GeV];#epsilon", 20, 0, 1);

TEfficiency* neutrons_w_hit_in_ECAL_vs_neutrino_E = new TEfficiency("eff5", "; reco E_#nu [GeV]; #epsilon", 24, 0, 6.);
TEfficiency* neutrons_w_compatible_cells_in_ECAL_vs_neutrino_E = new TEfficiency("eff6", "; reco E_#nu [GeV]; #epsilon", 24, 0, 6.);
TEfficiency* neutrons_w_compatible_cells_in_ECAL_vs_neutrino_E_ = new TEfficiency("eff7", "; reco E_#nu [GeV]; #epsilon", 24, 0, 6.);

TH1D* h_nof_cells_fired = new TH1D("", "", 50, 0, 50);
TH1D* h_nof_cells_fired_neutron = new TH1D("", "", 50, 0, 50);
TH1D* h_nof_cells_fired_antimu = new TH1D("", "", 50, 0, 50);

TH1D* adc_count_neutron = new TH1D("adc_neutron", "", 300, 0, 300);
TH1D* adc_count_antimu = new TH1D("adc_count_antimu", "", 300 , 0, 300);

TH1D* cell_reco_energy_neutron = new TH1D("cell_reco_energy_neutron", ";reco energy deposit;", 60, 0, 60); 
TH1D* cell_reco_energy_antimu = new TH1D("cell_reco_energy_antimu", ";reco energy deposit;", 60, 0, 60); 

TH1D* nof_events_w_zero_complete_cells = new TH1D("nof_events_w_zero_complete_cells", "", 2, 0, 1);

TH2D* adc1_vs_adc2_incomplete = new TH2D("adc1_vs_adc2_incomplete", "", 20, 0, 200, 20, 0, 200);

TH2D* best_time_residual_vs_any_space_residual = new TH2D("best_time_vs_any_space", "", 100, 0, 450, 100, -4, 4);
TH2D* any_time_residual_vs_best_space_residual = new TH2D("any_time_vs_best_space", "", 100, 0, 450, 100, -4, 4);

TH2D* res_time_vs_res_space_3cells = new TH2D("res_time_vs_res_space_3cells", "", 100, 0, 450, 100, -4, 4);
TH2D* res_time_vs_res_space_3cells_ = new TH2D("res_time_vs_res_space_3cells_", "", 100, 0, 450, 100, -4, 4);

TH1D* best_time_residuals = new TH1D("best_time_residuals", "", 300, -4, 4);
TH1D* best_space_residuals = new TH1D("best_space_residuals", "", 300, 0, 450);

// FLIGHT LENGTH

std::vector<double> GetBestTime(std::vector<int> mod_id,
                                std::vector<int> is_complete,
                                std::vector<double> time_residuals, 
                                std::vector<double> space_residuals,
                                std::vector<double> reco_energy_deposited,
                                bool& complete,
                                float& energy_deposited){
    double min_residual = 999999;
    double space_residual = -1.;
    double time_residual = -1.;
    
    for(size_t i=0; i< time_residuals.size(); i++){
        
        // if(mod_id[i] < 30) continue;
        if(!is_complete[i]) continue;

        float this_residual = fabs(time_residuals[i]);
        
        if(this_residual < min_residual){
            min_residual = this_residual;
            space_residual = space_residuals[i];
            time_residual = time_residuals[i];
            complete = true;
            energy_deposited = reco_energy_deposited[i];
        }
    }
    
    return {time_residual, space_residual};
}

std::vector<double> GetBestSpace(std::vector<int> mod_id,
                                 std::vector<int> is_complete,
                                 std::vector<double> time_residuals, 
                                 std::vector<double> space_residuals){
    double min_residual = 999999;
    double space_residual = -1.;
    double time_residual = -1.;
    
    for(size_t i=0; i< space_residuals.size(); i++){

        if(mod_id[i] < 30) continue;
        if(!is_complete[i]) continue;
        
        float this_residual = fabs(space_residuals[i]);
        
        if(this_residual < min_residual){
            min_residual = this_residual;
            space_residual = space_residuals[i];
            time_residual = time_residuals[i];
        }
    }
    
    return {time_residual, space_residual};
}

int efficiency_plots_test(){

    auto PRODUCTION_FOLDER = TString::Format("/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc_01/production_0016/preunfold/");
    auto DATA_PREUNFOLD = TString::Format("/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc_01/data_preunfold/");
    auto PRODUCTION = TString::Format("%smerged_output_analysis.0016.to.0016.root", DATA_PREUNFOLD.Data());
    auto FLUX = TString::Format("%sfluxes_merged.root", DATA_PREUNFOLD.Data());
    
    // single file to run tests
    // auto PRODUCTION = TString::Format("%sevents-in-SANDtracker.16000.to.16009.output_analysis.root", PRODUCTION_FOLDER.Data());
    
    std::cout << "Using file : " << PRODUCTION << "\n";
    
    TFile *file = TFile::Open(PRODUCTION.Data());
    TFile *fFlux = TFile::Open(FLUX.Data());

    TTree* flux = (TTree*)fFlux->Get("flux");
    if(!flux){
        return -1;
    }
    TH1F* hFlux = new TH1F("hFlux", ";neutrino E [GeV];Counts", nof_bin_nu_E, 0., 6.);
    flux->Draw("E >> hFlux");

    TTree* tree = (TTree*)file->Get("tAnalysis");

    tree->SetBranchAddress("CCQEonHydrogen", &CCQEonHydrogen);
    tree->SetBranchAddress("NofFinalStateChargedParticles", &NofFinalStateChargedParticles);
    // tree->SetBranchAddress("pass_nof_wires_cut", &pass_nof_wires_cut);
    tree->SetBranchAddress("nof_fired_wires", &nof_fired_wires);
    // tree->SetBranchAddress("candidate_signal_event", &candidate_signal_event);
    tree->SetBranchAddress("InteractionTarget", &InteractionTarget);
    tree->SetBranchAddress("InteractionVolume_short", &InteractionVolume_short);
    tree->SetBranchAddress("IncomingNeutrinoP4", &IncomingNeutrinoP4);
    tree->SetBranchAddress("Neutrino_reconstructed_P4_GeV", &Neutrino_reconstructed_P4_GeV);
    tree->SetBranchAddress("Interaction_vtxX", &Interaction_vtxX);
    tree->SetBranchAddress("Interaction_vtxY", &Interaction_vtxY);
    tree->SetBranchAddress("Interaction_vtxZ", &Interaction_vtxZ);
    tree->SetBranchAddress("FinalStateHadronicSystemTotal4Momentum", &FinalStateHadronicSystemTotal4Momentum);
    tree->SetBranchAddress("Antimuon_p_true", &Antimuon_p_true); // MeV
    tree->SetBranchAddress("Antimuon_reconstructed_P4", &Antimuon_reconstructed_P4); // MeV

    CellData cell_data;

    tree->SetBranchAddress("Fired_Cells_mod", &cell_data.Fired_Cells_mod);
    tree->SetBranchAddress("Fired_Cells_id", &cell_data.Fired_Cells_id);
    tree->SetBranchAddress("Fired_Cells_x", &cell_data.Fired_Cells_x);
    tree->SetBranchAddress("Fired_Cells_y", &cell_data.Fired_Cells_y);
    tree->SetBranchAddress("Fired_Cells_z", &cell_data.Fired_Cells_z);
    tree->SetBranchAddress("isCellComplete", &cell_data.isCellComplete);
    tree->SetBranchAddress("Fired_Cells_adc1", &cell_data.Fired_Cells_adc1);
    tree->SetBranchAddress("Fired_Cells_adc2", &cell_data.Fired_Cells_adc2);
    tree->SetBranchAddress("Fired_Cells_tdc1", &cell_data.Fired_Cells_tdc1);
    tree->SetBranchAddress("Fired_Cells_tdc2", &cell_data.Fired_Cells_tdc2);
    tree->SetBranchAddress("Fired_Cell_true_Hit_x", &cell_data.Fired_Cell_true_Hit_x);
    tree->SetBranchAddress("Fired_Cell_true_Hit_y", &cell_data.Fired_Cell_true_Hit_y);
    tree->SetBranchAddress("Fired_Cell_true_Hit_z", &cell_data.Fired_Cell_true_Hit_z);
    tree->SetBranchAddress("Fired_Cell_true_Hit_t", &cell_data.Fired_Cell_true_Hit_t);
    tree->SetBranchAddress("Fired_Cell_true_Hit_e", &cell_data.Fired_Cell_true_Hit_e);
    tree->SetBranchAddress("Fired_by_primary_antimu", &cell_data.Fired_by_primary_antimu);
    tree->SetBranchAddress("Fired_by_primary_neutron", &cell_data.Fired_by_primary_neutron);
    tree->SetBranchAddress("IsEarliestCell_neutron", &cell_data.IsEarliestCell_neutron);
    tree->SetBranchAddress("Reconstructed_HitPosition_x", &cell_data.Reconstructed_HitPosition_x);
    tree->SetBranchAddress("Reconstructed_HitPosition_y", &cell_data.Reconstructed_HitPosition_y);
    tree->SetBranchAddress("Reconstructed_HitPosition_z", &cell_data.Reconstructed_HitPosition_z);
    tree->SetBranchAddress("Reconstructed_HitTime", &cell_data.Reconstructed_HitTime);
    tree->SetBranchAddress("Reconstructed_Energy", &cell_data.Reconstructed_Energy);
    tree->SetBranchAddress("ExpectedNeutron_HitPosition_x_", &cell_data.ExpectedNeutron_HitPosition_x_);
    tree->SetBranchAddress("ExpectedNeutron_HitPosition_y_", &cell_data.ExpectedNeutron_HitPosition_y_);
    tree->SetBranchAddress("ExpectedNeutron_HitPosition_z_", &cell_data.ExpectedNeutron_HitPosition_z_);
    tree->SetBranchAddress("Expected_HitTime_", &cell_data.Expected_HitTime_);
    tree->SetBranchAddress("True_FlightLength", &cell_data.True_FlightLength);
    tree->SetBranchAddress("Reconstructed_FlightLength", &cell_data.Reconstructed_FlightLength);
    tree->SetBranchAddress("Residuals_HitTime_", &cell_data.Residuals_HitTime_);
    tree->SetBranchAddress("Residuals_HitSpace_", &cell_data.Residuals_HitSpace_);

    auto nof_entries = tree->GetEntries();

    std::cout << "Number of events in fiducial volume " << nof_entries << "\n";

    Counter counter;

    Init_Unfold();
    
    for (size_t i = 0; i < nof_entries; i++)
    {
        tree->GetEntry(i);

        if(i%100000 == 0 && i > 999) std::cout << "looping : " << i * 100 / nof_entries << "[%]\n"; 
        /**
            TRUE:
        */
        double IncomingNeutrino_energy = IncomingNeutrinoP4->T();
        double antimuon_true_energy = sqrt(Antimuon_p_true->Mag() * Antimuon_p_true->Mag() +  105.658 * 105.658);
        TVector3 Antimuon_p_true_Gev = (*Antimuon_p_true) * 1e-3;
        double antimuon_true_angle = Antimuon_p_true->Angle({0.,0.,1.});
        double hadron_syst_kin_energy = FinalStateHadronicSystemTotal4Momentum->T() -  0.939565;
        double Q2_true = (IncomingNeutrinoP4->Vect() - Antimuon_p_true_Gev).Mag2();
        bool is_in_graphite = (*InteractionVolume_short == "C_Target");
        bool is_in_plastic = (*InteractionVolume_short == "C3H6_Target");
        bool is_on_H = (*InteractionTarget == "proton");
        bool is_on_C = (*InteractionTarget == "C12");
        bool is_signal = (CCQEonHydrogen == 1);
        /**
            RECO:
        */
        double antimuon_reco_energy = Antimuon_reconstructed_P4->T();
        double Neutrino_reconstructed_energy_GeV = Neutrino_reconstructed_P4_GeV->T();
        double neutron_kin_energy_reco = PredictedNeutron_E_GeV -  0.939565;

        /**
            SELECTION:
        */

        std::vector<int> is_cell_compatible = are_cells_compatibles(*cell_data.Reconstructed_FlightLength,
                                                                    *cell_data.Reconstructed_HitTime, 
                                                                    *cell_data.Residuals_HitSpace_, 
                                                                    *cell_data.Residuals_HitTime_, 
                                                                    *cell_data.Reconstructed_Energy,
                                                                    *cell_data.IsEarliestCell_neutron); 
        bool pass_nof_wires_cut = (nof_fired_wires > 69);
        int nof_event_compatible_cells = std::accumulate(is_cell_compatible.begin(), is_cell_compatible.end(), 0.0);
        bool event_has_compatible_cells = (nof_event_compatible_cells > 0);
        bool event_has_charge_multi_1 = (NofFinalStateChargedParticles == 1); 
        bool is_event_selected = (pass_nof_wires_cut && event_has_charge_multi_1 && event_has_compatible_cells);
        // bool is_event_selected = (pass_nof_wires_cut && event_has_charge_multi_1 && event_has_compatible_cell);

        // ______________________________________________________________________________________________________

        SELECTION selection = DetermineSelection(is_event_selected, is_signal, is_in_graphite, is_in_plastic, is_on_C, is_on_H);

        antinu_hist.Fill(FIDUCIAL_VOLUME, SELECTION_NONE, REACTION_NONE, 0, IncomingNeutrino_energy);
        positive_mu_hist.Fill(FIDUCIAL_VOLUME, SELECTION_NONE, REACTION_NONE, 0, antimuon_true_energy);

        Fill_Final_Histos(selection, 
                          IncomingNeutrino_energy, 
                          Neutrino_reconstructed_energy_GeV, 
                          cell_data, 
                          counter,
                          i);

        if(!is_signal){
            continue;
        }
        //__________________________________________________________________________________

        positive_mu_hist.Fill(FIDUCIAL_VOLUME, SELECTION_NONE, CCQE_ON_H, 0, antimuon_true_energy);
        positive_mu_hist.Fill(FIDUCIAL_VOLUME, SELECTION_NONE, CCQE_ON_H, 1, antimuon_reco_energy);
        
        antinu_hist.Fill(FIDUCIAL_VOLUME, SELECTION_NONE, CCQE_ON_H, 0, IncomingNeutrino_energy);
        antinu_hist.Fill(FIDUCIAL_VOLUME, SELECTION_NONE, CCQE_ON_H, 1, Neutrino_reconstructed_energy_GeV);
        
        neutron_hist.Fill(FIDUCIAL_VOLUME, SELECTION_NONE, CCQE_ON_H, 0, hadron_syst_kin_energy);
        neutron_hist.Fill(FIDUCIAL_VOLUME, SELECTION_NONE, CCQE_ON_H, 1, neutron_kin_energy_reco);
        
        // // // true_nu_E_vs_true_n_kin_E[FIDUCIAL_VOLUME] -> Fill(IncomingNeutrino_energy, hadron_syst_kin_energy);
        // // // true_nu_E_vs_true_mu_E[FIDUCIAL_VOLUME] -> Fill(IncomingNeutrino_energy, antimuon_true_energy);
        // // // true_mu_E_vs_true_n_kin_E[FIDUCIAL_VOLUME] -> Fill(antimuon_true_energy, hadron_syst_kin_energy);

        if(!pass_nof_wires_cut){
            continue;
        }

        Increment_counter(counter, WIRES_CUT, is_signal, is_on_C, is_in_graphite, is_in_plastic, is_on_H);

        positive_mu_hist.Fill(WIRES_CUT, SELECTION_NONE, CCQE_ON_H, 0, antimuon_true_energy);
        positive_mu_hist.Fill(WIRES_CUT, SELECTION_NONE, CCQE_ON_H, 1, antimuon_reco_energy);

        antinu_hist.Fill(WIRES_CUT, SELECTION_NONE, CCQE_ON_H, 0, IncomingNeutrino_energy);
        antinu_hist.Fill(WIRES_CUT, SELECTION_NONE, CCQE_ON_H, 1, Neutrino_reconstructed_energy_GeV);

        neutron_hist.Fill(WIRES_CUT, SELECTION_NONE, CCQE_ON_H, 0, hadron_syst_kin_energy);
        neutron_hist.Fill(WIRES_CUT, SELECTION_NONE, CCQE_ON_H, 1, neutron_kin_energy_reco);
        
        // // true_nu_E_vs_true_n_kin_E[WIRES_CUT] -> Fill(IncomingNeutrino_energy, hadron_syst_kin_energy);
        // // true_nu_E_vs_true_mu_E[WIRES_CUT] -> Fill(IncomingNeutrino_energy, antimuon_true_energy);
        // // true_mu_E_vs_true_n_kin_E[WIRES_CUT] -> Fill(antimuon_true_energy, hadron_syst_kin_energy);

        if(!event_has_charge_multi_1){
            continue;
        }

        Increment_counter(counter, CHARGE_MULTIPLICITY, is_signal, is_on_C, is_in_graphite, is_in_plastic, is_on_H);

        positive_mu_hist.Fill(CHARGE_MULTIPLICITY, SELECTION_NONE, CCQE_ON_H, 0, antimuon_true_energy);
        positive_mu_hist.Fill(CHARGE_MULTIPLICITY, SELECTION_NONE, CCQE_ON_H, 1, antimuon_reco_energy);
        
        antinu_hist.Fill(CHARGE_MULTIPLICITY, SELECTION_NONE, CCQE_ON_H, 0, IncomingNeutrino_energy);
        antinu_hist.Fill(CHARGE_MULTIPLICITY, SELECTION_NONE, CCQE_ON_H, 1, Neutrino_reconstructed_energy_GeV);

        neutron_hist.Fill(CHARGE_MULTIPLICITY, SELECTION_NONE, CCQE_ON_H, 0, hadron_syst_kin_energy);
        neutron_hist.Fill(CHARGE_MULTIPLICITY, SELECTION_NONE, CCQE_ON_H, 1, neutron_kin_energy_reco);
        
        // // true_nu_E_vs_true_n_kin_E[CHARGE_MULTIPLICITY] -> Fill(IncomingNeutrino_energy, hadron_syst_kin_energy);
        // // true_nu_E_vs_true_mu_E[CHARGE_MULTIPLICITY] -> Fill(IncomingNeutrino_energy, antimuon_true_energy);
        // // true_mu_E_vs_true_n_kin_E[CHARGE_MULTIPLICITY] -> Fill(antimuon_true_energy, hadron_syst_kin_energy);

        if(!event_has_compatible_cells){
            neutrons_w_compatible_cells_in_ECAL_vs_neutrino_E_ -> Fill(0, IncomingNeutrino_energy);
            continue;
        }
        neutrons_w_compatible_cells_in_ECAL_vs_neutrino_E_ -> Fill(1, IncomingNeutrino_energy);
        
        positive_mu_hist.Fill(ECAL_COINCIDENCE, SELECTION_NONE, CCQE_ON_H, 0, antimuon_true_energy);
        positive_mu_hist.Fill(ECAL_COINCIDENCE, SELECTION_NONE, CCQE_ON_H, 1, antimuon_reco_energy);

        antinu_hist.Fill(ECAL_COINCIDENCE, SELECTION_NONE, CCQE_ON_H, 0, IncomingNeutrino_energy);
        antinu_hist.Fill(ECAL_COINCIDENCE, SELECTION_NONE, CCQE_ON_H, 1, Neutrino_reconstructed_energy_GeV);

        neutron_hist.Fill(ECAL_COINCIDENCE, SELECTION_NONE, CCQE_ON_H, 0, hadron_syst_kin_energy);
        neutron_hist.Fill(ECAL_COINCIDENCE, SELECTION_NONE, CCQE_ON_H, 1, neutron_kin_energy_reco);

        // // true_nu_E_vs_true_n_kin_E[ECAL_COINCIDENCE] -> Fill(IncomingNeutrino_energy, hadron_syst_kin_energy);
        // // true_nu_E_vs_true_mu_E[ECAL_COINCIDENCE] -> Fill(IncomingNeutrino_energy, antimuon_true_energy);
        // // true_mu_E_vs_true_n_kin_E[ECAL_COINCIDENCE] -> Fill(antimuon_true_energy, hadron_syst_kin_energy);

    }

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    /**
     * NOTE: compare true spectra of antinu, mu+ and neutron from signal at each stage of the analysis
    */
    positive_mu_hist.CompareStages("muon_spectra");
    antinu_hist.CompareStages("neutrino_spectra");
    neutron_hist.CompareStages("neutron_spectra");

    /**
    @relations: plot true relations
     */
    // plot_histograms2(true_nu_E_vs_true_n_kin_E[FIDUCIAL_VOLUME], "c1"); // nu vs n
    // plot_histograms2(true_nu_E_vs_true_mu_E[FIDUCIAL_VOLUME], "c2"); // nu vs mu
    // plot_histograms2(true_mu_E_vs_true_n_kin_E[FIDUCIAL_VOLUME], "c3"); // mu vs n
    
    //________________________________________________________________________________________
    // rate in plastic
    TH1D* ccqe_on_H_true = antinu_hist.GetHistogram(FIDUCIAL_VOLUME, SELECTION_NONE, CCQE_ON_H, 0);
    TH1D* selected_plastic_reco = antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_POSITIVE_PLASTIC, REACTION_ANY, 1);
    TH1D* ccqe_on_H_selected_plastic_true = antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_TRUE_POSITIVE, CCQE_ON_H, 0);
    TH1D* ccqe_on_H_selected_plastic_reco = antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_TRUE_POSITIVE, CCQE_ON_H, 1);
    TH1D* ccres_in_H_selected_plastic_reco = antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_PLASTIC_H, CCRES_ON_H, 1);
    TH1D* cc_on_C_selected_plastic_reco = antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_PLASTIC_C, CC_ON_CARBON, 1);
    TH1D* selected_graphite_reco = antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_GRAPHITE, CC_ON_CARBON, 1);
    TH1D* selected_plastic_carbon_subtracted = (TH1D*)selected_plastic_reco -> Clone("selected_plastic_carbon_subtracted");

    THStack* total_rate_plastic = new THStack("stack", "Stacked Histogram; #bar{#nu}_{#mu} Reco Energy [GeV]; Counts");
    total_rate_plastic -> Add(ccqe_on_H_selected_plastic_reco);
    total_rate_plastic -> Add(cc_on_C_selected_plastic_reco);
    total_rate_plastic -> Add(ccres_in_H_selected_plastic_reco);

    // to calculated errors on TH1D ratios
    ccqe_on_H_true -> Sumw2();
    selected_plastic_reco -> Sumw2();
    ccqe_on_H_selected_plastic_true -> Sumw2();
    ccqe_on_H_selected_plastic_reco -> Sumw2();
    ccres_in_H_selected_plastic_reco -> Sumw2();
    cc_on_C_selected_plastic_reco -> Sumw2();
    selected_graphite_reco -> Sumw2();
    selected_plastic_carbon_subtracted -> Sumw2();

    // SCALE
    selected_graphite_reco -> Scale(scale_factor);
    selected_plastic_carbon_subtracted -> Add(selected_graphite_reco, -1.);
    selected_plastic_carbon_subtracted->Add(ccres_in_H_selected_plastic_reco, -1.);

    // EFFICIENCY
    TH1D* total_selection_efficiency = (TH1D*)ccqe_on_H_selected_plastic_true->Clone("ratio");
    total_selection_efficiency->Divide(ccqe_on_H_true);
    
     
    // RECONSTRUCTED RATE
    TH1D* final_rate_reco = (TH1D*)selected_plastic_carbon_subtracted->Clone("final_rate_reco");
    final_rate_reco -> Divide(total_selection_efficiency);

    // NORMALIZED
    TH1D* ccqe_on_H_true_norm = (TH1D*)ccqe_on_H_true->Clone("ccqe_on_H_true_norm");
    TH1D* ccqe_on_H_selected_plastic_true_norm = (TH1D*)ccqe_on_H_selected_plastic_true->Clone("ccqe_on_H_selected_plastic_true_norm");
    TH1D* selected_plastic_carbon_subtracted_norm = (TH1D*)selected_plastic_carbon_subtracted->Clone("selected_plastic_carbon_subtracted_norm");
    TH1D* final_rate_reco_norm = (TH1D*)final_rate_reco->Clone("final_rate_reco_norm");
    TH1D* flux_norm = (TH1D*)hFlux->Clone("flux_norm");
    //
    ccqe_on_H_true_norm -> Sumw2();
    ccqe_on_H_selected_plastic_true_norm -> Sumw2();
    selected_plastic_carbon_subtracted_norm -> Sumw2();
    final_rate_reco_norm -> Sumw2();
    flux_norm -> Sumw2();
    //
    ccqe_on_H_true_norm->Scale(1. / ccqe_on_H_true->Integral());
    ccqe_on_H_selected_plastic_true_norm->Scale(1. / ccqe_on_H_selected_plastic_true->Integral());
    selected_plastic_carbon_subtracted_norm->Scale(1. / selected_plastic_carbon_subtracted->Integral());
    final_rate_reco_norm -> Scale(1. / final_rate_reco->Integral());
    flux_norm -> Scale(1. / flux_norm->Integral());

    // STYLE
    ccqe_on_H_true->SetLineColor(kRed);
    ccqe_on_H_selected_plastic_reco->SetFillColor(kBlue);
    cc_on_C_selected_plastic_reco->SetFillColor(kRed);
    ccres_in_H_selected_plastic_reco->SetFillColor(kGreen);

    // reconstructed style
    selected_plastic_reco->SetLineColor(kBlack);
    selected_plastic_reco->SetMarkerStyle(20);
    selected_plastic_reco->SetMarkerSize(0.8);
    selected_plastic_reco->SetMaximum(0.25);
    selected_graphite_reco->SetLineColor(kBlack);
    selected_graphite_reco->SetMarkerStyle(20);
    selected_graphite_reco->SetMarkerSize(0.8);
    selected_plastic_carbon_subtracted->SetLineColor(kBlack);
    selected_plastic_carbon_subtracted->SetMarkerStyle(20);
    selected_plastic_carbon_subtracted->SetMarkerSize(0.8);
    selected_plastic_carbon_subtracted_norm->SetLineColor(kBlack);
    selected_plastic_carbon_subtracted_norm->SetMarkerStyle(20);
    selected_plastic_carbon_subtracted_norm->SetMarkerSize(0.8);
    final_rate_reco->SetLineColor(kBlack);
    final_rate_reco->SetMarkerStyle(20);
    final_rate_reco->SetMarkerSize(0.8);
    final_rate_reco_norm->SetLineColor(kBlack);
    final_rate_reco_norm->SetMarkerStyle(20);
    final_rate_reco_norm->SetMarkerSize(0.8);
    final_rate_reco_norm->GetXaxis()->SetTitle("#bar{#nu}_{#mu} Reco Energy [GeV]");

    // ratio style
    TH1D* total_selection_efficiency_relative_flux = (TH1D*)ccqe_on_H_selected_plastic_true_norm->Clone("ratio");
    total_selection_efficiency_relative_flux->Divide(ccqe_on_H_true_norm);
    total_selection_efficiency_relative_flux->SetTitle(""); // Remove title for ratio plot
    total_selection_efficiency_relative_flux->GetYaxis()->SetTitle("Bin Content Ratio");
    total_selection_efficiency_relative_flux->GetXaxis()->SetTitle("#bar{#nu}_{#mu} Reco Energy [GeV]");
    total_selection_efficiency_relative_flux->GetYaxis()->SetNdivisions(505);
    total_selection_efficiency_relative_flux->GetYaxis()->SetTitleSize(0.1);
    total_selection_efficiency_relative_flux->GetYaxis()->SetTitleOffset(0.4);
    total_selection_efficiency_relative_flux->GetYaxis()->SetLabelSize(0.08);
    total_selection_efficiency_relative_flux->GetYaxis()->SetRangeUser(0, 0.3);
    total_selection_efficiency_relative_flux->GetXaxis()->SetTitleSize(0.1);
    total_selection_efficiency_relative_flux->GetXaxis()->SetTitleOffset(0.8);
    total_selection_efficiency_relative_flux->GetXaxis()->SetLabelSize(0.08);
    total_selection_efficiency_relative_flux->SetLineColor(kBlue);

    //________________________________________________________________________________________
    /**
     * THSTACK: stack plot of components in selected CCQE-on-H like in plastic
    */

    TCanvas* canvas_stack = new TCanvas("Stack", "Stack Plot", 900, 700);
    
    total_rate_plastic->Draw("E HIST");
    selected_plastic_reco->Draw("E1 SAME");

    TLegend* legendstack = new TLegend(0.7, 0.7, 0.9, 0.9);
    legendstack->AddEntry(ccqe_on_H_selected_plastic_reco, "CCQE on H (signal)", "f");
    legendstack->AddEntry(cc_on_C_selected_plastic_reco, "CC on Carbon (bkg)", "f");
    legendstack->AddEntry(ccres_in_H_selected_plastic_reco, "CCRES on H (bkg)", "f");
    legendstack->Draw();

    canvas_stack->Update();
    canvas_stack->SaveAs("plots/measured_rate_plastic_compontents.pdf");

    //_________________________________________________________________________________________
    /**
     * BACKGROUND: compare selected background in grafite vs selected backgroung of C in plastic
    */
    TCanvas* canvas_bkg = new TCanvas("c", "", 900, 700);

    // Define pads
    TPad* pad1_bkg = new TPad("pad1", "Main Plot", 0.0, 0.3, 1.0, 1.0); // Top pad
    TPad* pad2_bkg = new TPad("pad2", "Ratio Plot", 0.0, 0.0, 1.0, 0.3); // Bottom pad
    pad1_bkg->SetBottomMargin(0.02); // Reduce bottom margin for top pad
    pad2_bkg->SetTopMargin(0.02);   // Reduce top margin for bottom pad
    pad2_bkg->SetBottomMargin(0.3); // Increase bottom margin for labels
    pad1_bkg->Draw();
    pad2_bkg->Draw();

    // Draw the main plot
    pad1_bkg->cd();
    cc_on_C_selected_plastic_reco->Draw("E HIST");
    selected_graphite_reco->Draw("E1 SAME");

    // Add legend
    TLegend* legendbkg = new TLegend(0.6, 0.6, 0.9, 0.9);
    legendbkg->AddEntry(cc_on_C_selected_plastic_reco, "CC on Carbon Plastic (signal-like)", "f");
    legendbkg->AddEntry(selected_graphite_reco, "CC on Carbon Graphite (signal-like) X 4", "l");
    legendbkg->Draw();

    pad2_bkg->cd();
    TH1D* ratio_bkg = (TH1D*)cc_on_C_selected_plastic_reco->Clone("ratio");
    ratio_bkg->Divide(selected_graphite_reco);
    ratio_bkg->SetTitle(""); // Remove title for ratio plot
    ratio_bkg->GetYaxis()->SetTitle("Bin Content Ratio");
    ratio_bkg->GetXaxis()->SetTitle("Energy (GeV)");
    ratio_bkg->GetYaxis()->SetNdivisions(505);
    ratio_bkg->GetYaxis()->SetTitleSize(0.1);
    ratio_bkg->GetYaxis()->SetTitleOffset(0.4);
    ratio_bkg->GetYaxis()->SetLabelSize(0.08);
    ratio_bkg->GetXaxis()->SetTitleSize(0.1);
    ratio_bkg->GetXaxis()->SetTitleOffset(0.8);
    ratio_bkg->GetXaxis()->SetLabelSize(0.08);
    ratio_bkg->SetLineColor(kBlack);
    ratio_bkg->Draw("E1");

    canvas_bkg->Update();
    canvas_bkg->SaveAs("plots/graphite_2_plastic_bkg_comparison.pdf");
    //_________________________________________________________________________________________

    /**
     * SUBTRACTION:
    */

    TLegend* legendsub = new TLegend(0.6, 0.6, 0.9, 0.9);
    legendsub->AddEntry(ccqe_on_H_true, "#Phi_{#bar{#nu}_{#mu}} #sigma_{S} (CCQE on H)", "l");
    legendsub->AddEntry(ccqe_on_H_selected_plastic_true, "#Phi_{#bar{#nu}_{#mu}} #sigma_{S} #epsilon_{S} (CCQE on H like)", "lp");
    legendsub->AddEntry(selected_plastic_carbon_subtracted, "#Phi_{#bar{#nu}_{#mu}} #sigma_{S} #epsilon_{S} R_{det} (w/ bkg subtraction)", "lp");

    TCanvas* canvas_sub = new TCanvas("sub", "", 900, 700);
    TPad* pad1 = new TPad("pad1", "Main Plot", 0.0, 0.3, 1.0, 1.0); // Top pad
    TPad* pad2 = new TPad("pad2", "Ratio Plot", 0.0, 0.0, 1.0, 0.3); // Bottom pad

    pad1->SetBottomMargin(0.02); // Reduce bottom margin for top pad
    pad2->SetTopMargin(0.02);   // Reduce top margin for bottom pad
    pad2->SetBottomMargin(0.3); // Increase bottom margn for labels

    canvas_sub->cd();
    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    ccqe_on_H_true_norm->Draw("E HIST");
    ccqe_on_H_selected_plastic_true_norm->Draw("E HIST SAME");
    selected_plastic_carbon_subtracted_norm->Draw("E1 SAME");
    // hFlux -> Draw("H SAME");
    legendsub->Draw("SAME");

    pad2->cd();
    TH1D* effect_efficiency_and_detector_and_subtraction = (TH1D*)selected_plastic_carbon_subtracted_norm->Clone("ratio");
    effect_efficiency_and_detector_and_subtraction->Divide(ccqe_on_H_true_norm);
    effect_efficiency_and_detector_and_subtraction ->SetLineColor(kBlack);
    effect_efficiency_and_detector_and_subtraction->SetMarkerStyle(20);
    effect_efficiency_and_detector_and_subtraction->SetMarkerSize(0.8);
    total_selection_efficiency_relative_flux->Draw("E HIST");
    effect_efficiency_and_detector_and_subtraction->Draw("E1 SAME");

    canvas_sub->Update();

    /***
    POSTSUBTRACTION:
     */
    TCanvas* canvas_post_sub = new TCanvas("Stack", "Stack Plot", 900, 700);
    
    ccqe_on_H_selected_plastic_reco->GetXaxis()->SetTitle("#bar{#nu}_{#mu} Reco Energy [GeV]");    
    ccqe_on_H_selected_plastic_reco->Draw("E HIST");
    selected_plastic_carbon_subtracted->Draw("E1 SAME");

    TLegend* legendpost = new TLegend(0.7, 0.7, 0.9, 0.9);
    legendpost->AddEntry(ccqe_on_H_selected_plastic_reco, "CCQE on H (signal)", "f");
    legendpost->AddEntry(selected_plastic_carbon_subtracted, "plastic w/ bkg subtraction", "lp");
    legendpost->Draw();

    canvas_post_sub->Update();
    canvas_post_sub->SaveAs("plots/measured_rate_plastic_post_subtraction.pdf");

    /***
    RESPONSEMATRIX: response matrix true ccqe on H selected
    the response matrix account also for efficiency
     */
    TH2* hResponse = response_matrix[SIGNAL_SAMPLE1]->Hresponse();
    hResponse->SetTitle("response matrix true CCQE on H selected; true neutrino energy");
    TCanvas* canvas_response = new TCanvas("canvas_response", "", 900, 700);
    hResponse->Draw("col z");
    canvas_response->SaveAs("plots/response_matrix_ccqe_On_h_selected.pdf");

    /***
    UNFOLDING:
    - ccqe_on_H_true: true event rate as function of true neutrino energy
    - selected_plastic_carbon_subtracted: reconstructed event rate from carbon sub as function of reconstructed neutrino energy
     */
    RooUnfoldBayes   unfold(response_matrix[SIGNAL_SAMPLE1], selected_plastic_carbon_subtracted, 30);
    auto selected_plastic_carbon_subtracted_UNFOLDED = (TH1D*) unfold.Hreco();
    TCanvas* canvas_unfold = new TCanvas("unfolded", "", 900, 700);
    // unfold.PrintTable (cout, ccqe_on_H_true);
    ccqe_on_H_true->GetXaxis()->SetTitle("#bar{#nu}_{#mu} Reco Energy [GeV]");
    ccqe_on_H_true->Draw("E HIST");
    selected_plastic_carbon_subtracted_UNFOLDED->Draw("E HIST SAME");
    selected_plastic_carbon_subtracted->Draw("E1 SAME");

    TLegend* legend_unf = new TLegend(0.58, 0.58, 0.9, 0.9);
    legend_unf->AddEntry(ccqe_on_H_true, "#m_H #Phi_{#bar{#nu}_{#mu}} #sigma_{S} (CCQE on H)", "l");
    legend_unf->AddEntry(selected_plastic_carbon_subtracted, "plastic w/ bkg subtraction", "lp");
    legend_unf->AddEntry(selected_plastic_carbon_subtracted_UNFOLDED, "unfolded", "lp");
    legend_unf->Draw("SAME");
    
    canvas_unfold->Update();
    canvas_unfold->SaveAs("plots/unfolded_rate.pdf");


    /***
    CORRECTEDRATE:
    */
    // TCanvas* canvas_final_rate = new TCanvas("canvas_final_rate", "", 900, 700);
    // final_rate_reco_norm -> Draw("E1");
    // ccqe_on_H_true_norm -> Draw("E HIST SAME");

    // TLegend* legend_rates = new TLegend(0.7, 0.7, 0.9, 0.9);
    // legend_rates->AddEntry(ccqe_on_H_true_norm, "#Phi_{#bar{#nu}_{#mu}} #sigma_{S} (CCQE on H)", "l");
    // legend_rates->AddEntry(final_rate_reco_norm, "#Phi_{#bar{#nu}_{#mu}} #sigma_{S} R_{det} (w/ bkg subtraction)", "l");
    // legend_rates->Draw();

    /***
    RELATIVEFLUX:
    */
    TCanvas* canvas_flux = new TCanvas("canvas_flux", "", 900, 700);
    TPad* pad1_flux = new TPad("pad1", "Main Plot", 0.0, 0.3, 1.0, 1.0); // Top pad
    TPad* pad2_flux = new TPad("pad2", "Ratio Plot", 0.0, 0.0, 1.0, 0.3); // Bottom pad

    pad1_flux->SetBottomMargin(0.02); // Reduce bottom margin for top pad
    pad2_flux->SetTopMargin(0.02);   // Reduce top margin for bottom pad
    pad2_flux->SetBottomMargin(0.3); // Increase bottom margn for labels

    canvas_flux->cd();
    pad1_flux->Draw();
    pad2_flux->Draw();

    pad1_flux->cd();

    hFlux -> Draw("E HIST");
    // flux_norm -> Draw("E HIST");
    ccqe_on_H_true -> Draw("E HIST SAME");
    // ccqe_on_H_true_norm -> Draw("E HIST SAME");

    TLegend* legend_flux = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_flux->AddEntry(hFlux, "#Phi_{#bar{#nu}_{#mu}} ", "l");
    legend_flux->AddEntry(ccqe_on_H_true, "#Phi_{#bar{#nu}_{#mu}} #sigma_{S} (CCQE on H)", "l");
    legend_flux->Draw();

    pad2_flux->cd();
    TH1D* xsec_shape = (TH1D*)ccqe_on_H_true->Clone("xsec_shape"); 
    // TH1D* xsec_shape = (TH1D*)ccqe_on_H_true_norm->Clone("xsec_shape"); 
    xsec_shape->Divide(hFlux);
    xsec_shape->Draw("E HIST");
    canvas_flux->Update();

    /**
    UNFOLDEDCROSSEX:
         */
    RooUnfoldResponse response(xsec_shape, hResponse);
    RooUnfoldBinByBin unfoldBinByBin(&response, xsec_shape);

    // Perform the unfolding and get the result
    auto unfoldedBinByBin = (TH1D*) unfoldBinByBin.Hreco();

    // Set up a canvas for plotting
    TCanvas* canvas_xsec = new TCanvas("canvas_xsec", "", 900, 700);
    xsec_shape->Draw("E HIST");
    unfoldedBinByBin->Draw("E HIST SAME");

    // Add a legend for clarity
    TLegend* legend = new TLegend(0.75, 0.75, 0.9, 0.9);
    legend->AddEntry(xsec_shape, "Input Spectrum", "l");
    legend->AddEntry(unfoldedBinByBin, "Unfolded Spectrum", "l");
    legend->Draw();

    // Save the canvas as a PDF
    canvas_xsec->Print("plots/unfolded_XSEC.pdf");

    /**
    RECOFLUX: 
     */
    TCanvas* canvas_flux_ = new TCanvas("canvas_flux_", "", 900, 700);
    hFlux -> Draw("E HIST");
    // flux_norm -> Draw("E HIST");
    TH1D* flux_reco = (TH1D*)selected_plastic_carbon_subtracted_UNFOLDED->Clone("flux_reco");
    // TH1D* flux_norm_reco = (TH1D*)final_rate_reco_norm->Clone("flux_norm_reco");
    flux_reco -> Divide(xsec_shape);
    flux_reco->SetLineColor(kBlack);
    flux_reco->SetMarkerStyle(20);
    flux_reco->SetMarkerSize(0.8);
    flux_reco -> Draw("E1 SAME");
    canvas_flux_->SaveAs("plots/reconstructed_flux.pdf");

    /**
    SYSTEMATICSFLUX:
    */
    TCanvas* canvas_systematics_ = new TCanvas("canvas_systematics_", "", 900, 700);

    // Creazione di un istogramma per memorizzare le differenze relative
    TH1D* flux_systematics = (TH1D*)hFlux->Clone("flux_systematics");
    flux_systematics->Add(flux_reco, -1.0); // Differenza tra hFlux e flux_reco
    flux_systematics->Divide(flux_reco);   // Calcolo della differenza relativa

    // Prendi il valore assoluto dei contenuti
    for (int i = 1; i <= flux_systematics->GetNbinsX(); ++i) {
        double content = flux_systematics->GetBinContent(i);
        flux_systematics->SetBinContent(i, std::abs(content)); // Valore assoluto
    }

    // Calcola i limiti dell'asse Y
    double min_y = 0;
    double max_y = 0;
    for (int i = 1; i <= flux_systematics->GetNbinsX(); ++i) {
        double content = flux_systematics->GetBinContent(i);
        double error = flux_systematics->GetBinError(i);
        max_y = std::max(max_y, content + error);
    }
    max_y *= 1.1; // Aggiungi un margine al massimo

    // Imposta i limiti dell'asse Y
    flux_systematics->SetMinimum(min_y);
    flux_systematics->SetMaximum(max_y);

    // Disegna l'istogramma
    flux_systematics->GetYaxis()->SetRangeUser(0., 1.);
    flux_systematics->SetTitle(";#bar{#nu}_{#mu} Reco Energy [GeV]; Flux Systematic Uncertainties");
    flux_systematics->SetLineColor(kBlack);
    flux_systematics->Draw("HIST"); // Disegna prima l'istogramma

    // Creazione e disegno della banda di errore
    for (int i = 1; i <= flux_systematics->GetNbinsX(); ++i) {
        double x_low = flux_systematics->GetBinLowEdge(i);
        double x_high = x_low + flux_systematics->GetBinWidth(i);
        double y = flux_systematics->GetBinContent(i);
        double error = flux_systematics->GetBinError(i);

        double y_min = std::max(0.0, y - error);
        double y_max = y + error;

        // Crea la banda di errore
        TBox* box = new TBox(x_low, y_min, x_high, y_max);
        box->SetFillColorAlpha(kRed, 0.35); // Rosso chiaro con trasparenza
        box->SetLineWidth(0);
        box->Draw("SAME");
    }

    // Ridisegna l'istogramma sopra le bande di errore
    flux_systematics->Draw("HIST SAME");

    // Salva il canvas
    canvas_systematics_->SaveAs("plots/systematics_on_reconstructed_flux.pdf");



    /***
    RECAP: final recap on number of events
    */

    PrintCounter(counter);

    /**
    UNFOLDMAGIC:
     */
    PrintAllUnfoldInfos();

    return 0;
}