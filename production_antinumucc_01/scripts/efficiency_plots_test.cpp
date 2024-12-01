#include <string>
#include <vector>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

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
const double scale_factor = 4.084;

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
    uint nof_bins = 12;
    double low = 0.;
    double up = 6.;
    h_bins(uint bins, double l, double u) : nof_bins(bins), low(l), up(u) {};
};

h_bins Get_h_bins(PARTICLE p)
{
    if(p == ANTIMUON)
    {
        return {24, 0., 6000.};
    }else if(p == NEUTRINO)
    {
        return {24, 0., 6.};
    }else if(p == NEUTRON){
        return {40, 0., 1.};
    }else
    {
        return {24, 0., 6.};
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
            auto h3 = energy_spectrum[CHARGE_MULTIPLICITY][SELECTION_NONE][CCQE_ON_H][0];
            auto h4 = energy_spectrum[ECAL_COINCIDENCE][SELECTION_NONE][CCQE_ON_H][0];
        
            h1 -> SetLineColor(kRed);
            h2 -> SetLineColor(kMagenta);
            h3 -> SetLineColor(kBlue);
            h4 -> SetLineColor(kBlack);
        
            h1 -> Draw("E HIST");
            h2 -> Draw("E HIST SAME");
            h3 -> Draw("E HIST SAME");
            h4 -> Draw("E HIST SAME");

            TLegend* legend = new TLegend(0.5, 0.5, 0.9, 0.9);
            std::string particle_name = GetParticleName();
            legend->AddEntry(h1, TString::Format("%s (STAGE FIDUCIAL)", particle_name.c_str()).Data(), "l");
            legend->AddEntry(h2, TString::Format("%s (STAGE WIRES CUT)", particle_name.c_str()).Data(), "l");
            legend->AddEntry(h3, TString::Format("%s (STAGE CHARGE MULTI)", particle_name.c_str()).Data(), "l");
            legend->AddEntry(h4, TString::Format("%s (STAGE ECAL COINC.)", particle_name.c_str()).Data(), "l");
            
            legend->Draw();
            
            // Calculate the ratio
            pad2->cd();
            h1->Sumw2();
            h2->Sumw2();
            h3->Sumw2();
            h4->Sumw2();
            TH1D* ratio21 = (TH1D*)h2->Clone("ratio21");
            TH1D* ratio32 = (TH1D*)h3->Clone("ratio32");
            TH1D* ratio43 = (TH1D*)h4->Clone("ratio43");
            
            ratio21->Divide(h1);
            ratio32->Divide(h2);
            ratio43->Divide(h3);
            
            ratio21->SetTitle(""); // Remove title for ratio plot
            
            ratio21->SetLineColor(kMagenta);
            ratio32->SetLineColor(kBlue);
            ratio43->SetLineColor(kBlack);
            
            ratio21->GetYaxis()->SetTitle("Cut Efficiency");
            std::string x_label = (particle == ANTIMUON) ? "Energy (MeV)" : "Energy (GeV)";
            ratio21->GetXaxis()->SetTitle(x_label.c_str());
            ratio21->GetYaxis()->SetNdivisions(505);
            ratio21->GetYaxis()->SetTitleSize(0.1);
            ratio21->GetYaxis()->SetTitleOffset(0.4);
            ratio21->GetYaxis()->SetLabelSize(0.08);
            ratio21->GetXaxis()->SetTitleSize(0.1);
            ratio21->GetXaxis()->SetTitleOffset(0.8);
            ratio21->GetXaxis()->SetLabelSize(0.08);
            
            ratio21->Draw("E1");
            ratio32->Draw("E1 SAME");
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


void Fill_Final_Histos(SELECTION selection, 
                      double true_nu_E, double reco_nu_E
                    //   double true_mu_E, double reco_mu_E,
                    //   double true_n_K, double reco_n_K,
                      ){
    switch (selection)
    {
        case SELECTED_TRUE_POSITIVE:
            
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_POSITIVE_PLASTIC, REACTION_ANY, 0, true_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_POSITIVE_PLASTIC, REACTION_ANY, 1, reco_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_TRUE_POSITIVE, CCQE_ON_H, 0, true_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_TRUE_POSITIVE, CCQE_ON_H, 1, reco_nu_E);
            
            break;

        case FALSE_NEGATIVE:
            break;

        // FALSE POSITIVES
        case SELECTED_FALSE_POSITIVE_GRAPHITE:
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_GRAPHITE, CC_ON_CARBON, 0, true_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_GRAPHITE, CC_ON_CARBON, 1, reco_nu_E);

            break;

        case SELECTED_FALSE_POSITIVE_PLASTIC_C:
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_POSITIVE_PLASTIC, REACTION_ANY, 0, true_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_POSITIVE_PLASTIC, REACTION_ANY, 1, reco_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_PLASTIC_C, CC_ON_CARBON, 0, true_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_PLASTIC_C, CC_ON_CARBON, 1, reco_nu_E);

            break;

        case SELECTED_FALSE_POSITIVE_PLASTIC_H:
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_POSITIVE_PLASTIC, REACTION_ANY, 0, true_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_POSITIVE_PLASTIC, REACTION_ANY, 1, reco_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_PLASTIC_H, CCRES_ON_H, 0, true_nu_E);
            antinu_hist.Fill(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_PLASTIC_H, CCRES_ON_H, 1, reco_nu_E);

            break;

        case SELECTED_FALSE_POSITIVE_OTHER:
            break;

        // TRUE NEGATIVES
        case TRUE_NEGATIVE_GRAPHITE:
            break;

        case TRUE_NEGATIVE_PLASTIC_C:
            break;

        case TRUE_NEGATIVE_PLASTIC_H:
            break;

        case TRUE_NEGATIVE_OTHERS:
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
                                      const std::vector<double>& reco_energy
                                    //   const std::vector<int>& is_complete
                                      ){
    /*
      Input: time/space residuals between the expected and the reconstructed hit in the cell
      and look for those cell compatible in the space and in the time
    */
    std::vector<int> isCompatible(space_residuals.size());
    
    for (size_t i = 0; i < space_residuals.size(); i++){
        if((space_residuals[i] == -999.)){
            isCompatible[i] = 0;
        }else{
            bool space_pass = (fabs(space_residuals[i]) <= 150);
            float expected_sigma = 0.05 / sqrt(reco_energy[i]*1e-3);
            bool time_pass = (fabs(time_residuals[i])) < 3*expected_sigma;
            // if(space_pass && time_pass && tof[i]>15.)
            if(space_pass && time_pass && flight_length[i]>500.)
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
    TH1F* hFlux = new TH1F("hFlux", "Distribuzione di E;E (GeV);Conteggi", 24, 0, 6);
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

    int nof_events_w_zero_complete_cells = 0;
    
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

        // bool is_complete = false;
        // float reco_energy_deposit = 0;
        
        // std::vector<double> bests_time = GetBestTime(*cell_data.Fired_Cells_mod, 
        //                                              *cell_data.isCellComplete, 
        //                                              *cell_data.Residuals_HitTime_, 
        //                                              *cell_data.Residuals_HitSpace_,
        //                                              *cell_data.Reconstructed_Energy,
        //                                              is_complete,
        //                                              reco_energy_deposit);

        // std::vector<double> bests_space = GetBestSpace(*cell_data.Fired_Cells_mod, 
        //                                                *cell_data.isCellComplete, 
        //                                                *cell_data.Residuals_HitTime_, 
        //                                                *cell_data.Residuals_HitSpace_
        //                                                );
        
        // double best_residual_time = bests_time[0];
        // double best_residual_space = bests_space[0];
        // bool event_has_compatible_cell = has_compatible_cell(best_residual_time, 
        //                                                      best_residual_space, 
        //                                                      is_complete,
        //                                                      reco_energy_deposit);

        // int nof_cells_fired_by_neutron = std::accumulate(cell_data.Fired_by_primary_neutron->begin(),
        //                                                  cell_data.Fired_by_primary_antimu->end(),
        //                                                  0.);
        std::vector<int> is_cell_compatible = are_cells_compatibles(*cell_data.Reconstructed_FlightLength,
                                                                    *cell_data.Reconstructed_HitTime, 
                                                                    *cell_data.Residuals_HitSpace_, 
                                                                    *cell_data.Residuals_HitTime_, 
                                                                    *cell_data.Reconstructed_Energy); 
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

        Fill_Final_Histos(selection, IncomingNeutrino_energy, Neutrino_reconstructed_energy_GeV);

        // if(is_on_C)
        // {
        //     for(size_t j = 0; j < cell_data.Fired_by_primary_neutron->size(); j++)
        //     {
        //         if(cell_data.isCellComplete->at(j))
        //         {
        //             res_time_vs_res_space_3cells_ -> Fill(cell_data.Residuals_HitSpace_->at(j), cell_data.Residuals_HitTime_->at(j));
        //         }
        //     }
        // }
        
        if(!is_signal){
            continue;
        }

        // std::vector<double> event_cell_best_time = GetBestTime(*cell_data.Fired_Cells_mod, *cell_data.isCellComplete, *cell_data.Residuals_HitTime_, *cell_data.Residuals_HitSpace_, *cell_data.Reconstructed_Energy,
        //                                              is_complete,reco_energy_deposit);
        // std::vector<double> event_cell_best_space = GetBestSpace(*cell_data.Fired_Cells_mod, *cell_data.isCellComplete, *cell_data.Residuals_HitTime_, *cell_data.Residuals_HitSpace_);
        
        // best_time_residual_vs_any_space_residual -> Fill(event_cell_best_time[1], event_cell_best_time[0]);
        // any_time_residual_vs_best_space_residual -> Fill(event_cell_best_space[1], event_cell_best_space[0]);

        // best_time_residuals -> Fill(event_cell_best_time[0]);
        // best_space_residuals -> Fill(event_cell_best_space[0]);

        // int nof_tot_fired_cells = cell_data.Fired_Cells_mod->size();
        // int nof_tot_fired_cells_neutron = std::accumulate(cell_data.Fired_by_primary_neutron->begin(), cell_data.Fired_by_primary_neutron->end(), 0.0);
        // int nof_tot_fired_cells_antimu = std::accumulate(cell_data.Fired_by_primary_antimu->begin(), cell_data.Fired_by_primary_antimu->end(), 0.0);

        // int counter_ev_w_zero_complete_cells = 0;

        // h_nof_cells_fired -> Fill(nof_tot_fired_cells);
        // h_nof_cells_fired_neutron -> Fill(nof_tot_fired_cells_neutron);
        // h_nof_cells_fired_antimu -> Fill(nof_tot_fired_cells_antimu);

        // std::cout << "CCQE on H event, tot nof cells : " << nof_tot_fired_cells
        //           << ", fired by neutron " << nof_tot_fired_cells_neutron
        //           << ", fired by antimu " << nof_tot_fired_cells_antimu
        //           << "\n";
        
        // bool found_neutron_complete = false;
        // for(size_t j = 0; j < cell_data.Fired_by_primary_neutron->size(); j++)
        // {
        //     if(cell_data.Fired_by_primary_neutron->at(j))
        //     {
        //         if(cell_data.isCellComplete->at(j))
        //         {
        //             found_neutron_complete = true;
        //             if(1)
        //             {
        //                 res_time_vs_res_space_3cells -> Fill(cell_data.Residuals_HitSpace_->at(j), cell_data.Residuals_HitTime_->at(j));
        //             }
        //         }
        //         adc_count_neutron -> Fill(cell_data.Fired_Cells_adc1->at(j));
        //         adc_count_neutron -> Fill(cell_data.Fired_Cells_adc2->at(j));
        //         cell_reco_energy_neutron -> Fill(cell_data.Reconstructed_Energy->at(j));
        //     }else{
        //         adc_count_antimu -> Fill(cell_data.Fired_Cells_adc1->at(j));
        //         adc_count_antimu -> Fill(cell_data.Fired_Cells_adc2->at(j));
        //         cell_reco_energy_antimu -> Fill(cell_data.Reconstructed_Energy->at(j));
        //     }
        //     if(!(cell_data.isCellComplete->at(j)))
        //     {
        //         adc1_vs_adc2_incomplete -> Fill(cell_data.Fired_Cells_adc1->at(j), cell_data.Fired_Cells_adc2->at(j));
        //     }

        // }

        // if(!found_neutron_complete) nof_events_w_zero_complete_cells++;

        // // FILL EFFICIENCIES __________________________________________________________________________________
        
        // int nof_cell_fired_by_neutron = 0;
        // int nof_complete_cell_fired_by_neutron = 0;
        // int nof_compatible_cells_fired_by_neutron = 0;
        
        // for(size_t j = 0; j < is_cell_fired_by_neutron->size(); j++)
        // {
        //     if(!(is_cell_fired_by_neutron->at(j))) continue;
            
        //     nof_cell_fired_by_neutron++;
            
        //     if(!(isCellComplete->at(j))) continue;

        //     nof_complete_cell_fired_by_neutron++;

        //     if(!(IsCompatible->at(j))) continue;

        //     nof_compatible_cells_fired_by_neutron++;
        // }

        // if(nof_cell_fired_by_neutron > 0)
        // {
        //     neutrons_w_hit_in_ECAL -> Fill(1, hadron_syst_kin_energy);
        //     neutrons_w_hit_in_ECAL_vs_neutrino_E -> Fill(1, IncomingNeutrino_energy);
        // }else{
        //     neutrons_w_hit_in_ECAL -> Fill(0, hadron_syst_kin_energy);
        //     neutrons_w_hit_in_ECAL_vs_neutrino_E -> Fill(0, IncomingNeutrino_energy);
        // }

        // if(nof_complete_cell_fired_by_neutron > 0)
        // {
        //     neutrons_w_complete_cells_in_ECAL -> Fill(1, hadron_syst_kin_energy);
        // }else{
        //     neutrons_w_complete_cells_in_ECAL -> Fill(0, hadron_syst_kin_energy);
        // }
        
        // if(nof_compatible_cells_fired_by_neutron > 0)
        // {
        //     neutrons_w_compatible_cells_in_ECAL -> Fill(1, hadron_syst_kin_energy);
        //     neutrons_w_compatible_cells_in_ECAL_vs_neutrino_E -> Fill(1, Neutrino_reconstructed_energy_GeV);
        // }else{
        //     neutrons_w_compatible_cells_in_ECAL -> Fill(0, hadron_syst_kin_energy);
        //     neutrons_w_compatible_cells_in_ECAL_vs_neutrino_E -> Fill(0, Neutrino_reconstructed_energy_GeV);
        // }
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
        positive_mu_hist.Fill(CHARGE_MULTIPLICITY, SELECTION_NONE, CCQE_ON_H, 0, antimuon_true_energy);
        positive_mu_hist.Fill(CHARGE_MULTIPLICITY, SELECTION_NONE, CCQE_ON_H, 1, antimuon_reco_energy);
        
        antinu_hist.Fill(CHARGE_MULTIPLICITY, SELECTION_NONE, CCQE_ON_H, 0, IncomingNeutrino_energy);
        antinu_hist.Fill(CHARGE_MULTIPLICITY, SELECTION_NONE, CCQE_ON_H, 1, Neutrino_reconstructed_energy_GeV);

        neutron_hist.Fill(CHARGE_MULTIPLICITY, SELECTION_NONE, CCQE_ON_H, 0, hadron_syst_kin_energy);
        neutron_hist.Fill(CHARGE_MULTIPLICITY, SELECTION_NONE, CCQE_ON_H, 1, neutron_kin_energy_reco);
        
        // // true_nu_E_vs_true_n_kin_E[CHARGE_MULTIPLICITY] -> Fill(IncomingNeutrino_energy, hadron_syst_kin_energy);
        // // true_nu_E_vs_true_mu_E[CHARGE_MULTIPLICITY] -> Fill(IncomingNeutrino_energy, antimuon_true_energy);
        // // true_mu_E_vs_true_n_kin_E[CHARGE_MULTIPLICITY] -> Fill(antimuon_true_energy, hadron_syst_kin_energy);
        for(size_t j = 0; j < cell_data.Fired_by_primary_neutron->size(); j++)
        {
            if(!(cell_data.Fired_by_primary_neutron->at(j))) continue;
            if(!(cell_data.isCellComplete->at(j))) continue;
            // double reco_cell_space_res = sqrt(pow((Fired_Cell_true_Hit_x->at(j) - Reconstructed_HitPosition_x->at(j)), 2) 
            //                                 + pow((Fired_Cell_true_Hit_y->at(j) - Reconstructed_HitPosition_y->at(j)), 2) 
            //                                 + pow((Fired_Cell_true_Hit_z->at(j) - Reconstructed_HitPosition_z->at(j)), 2));
            // double reco_cell_space_res = sqrt(pow((Reconstructed_HitPosition_x->at(j) - ExpectedNeutron_HitPosition_x_->at(j)), 2) 
            //                                 + pow((Reconstructed_HitPosition_y->at(j) - ExpectedNeutron_HitPosition_y_->at(j)), 2) 
            //                                 + pow((Reconstructed_HitPosition_z->at(j) - ExpectedNeutron_HitPosition_z_->at(j)), 2));
            neutron_space_time_residulas -> Fill(cell_data.Residuals_HitSpace_->at(j), cell_data.Residuals_HitTime_->at(j));
        }

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


    //________________________________________________________________________________________
    /**
     * THSTACK: stack plot of components in selected CCQE-on-H like in plastic
    */
    THStack* stack = new THStack("stack", "Stacked Histogram; Reco Energy (GeV); Counts");

    int colors[REACTION_NONE] = {kRed, kBlue, kGreen};

    auto h1 = antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_TRUE_POSITIVE, CCQE_ON_H, 1);
    h1->SetFillColor(kBlue);

    auto h2 = antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_PLASTIC_C, CC_ON_CARBON, 1);
    h2->SetFillColor(kRed);

    auto h3 = antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_PLASTIC_H, CCRES_ON_H, 1);
    h3->SetFillColor(kGreen);

    stack->Add(h1);
    stack->Add(h2);
    stack->Add(h3);
    
    TCanvas* canvas_stack = new TCanvas("Stack", "Stack Plot", 900, 700);
    
    stack->Draw("E HIST");
    auto hm = antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_POSITIVE_PLASTIC, REACTION_ANY, 1);
    hm->SetLineColor(kBlack);
    hm->SetMarkerStyle(20);
    hm->SetMarkerSize(0.8);
    hm->SetMaximum(0.25);
    hm->Draw("E1 SAME");

    // Add legend
    TLegend* legendstack = new TLegend(0.7, 0.7, 0.9, 0.9);
    legendstack->AddEntry(h1, "CCQE on H (signal)", "f");
    legendstack->AddEntry(h2, "CC on Carbon (bkg)", "f");
    legendstack->AddEntry(h3, "CCRES on H (bkg)", "f");
    legendstack->Draw();

    // Update canvas
    canvas_stack->Update();

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
    auto h2_ = antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_GRAPHITE, CC_ON_CARBON, 1);
    h2_->SetLineColor(kBlue);
    
    h2->Draw("E HIST");
    h2_->Draw("E1 SAME");

    // Add legend
    TLegend* legendbkg = new TLegend(0.6, 0.6, 0.9, 0.9);
    legendbkg->AddEntry(h2, "CC on Carbon Plastic (signal-like)", "f");
    legendbkg->AddEntry(h2_, "CC on Carbon Graphite (signal-like)", "l");
    legendbkg->Draw();

    // Calculate the ratio
    pad2_bkg->cd();
    h2->Sumw2();
    h2_->Sumw2();
    h2_ -> SetLineColor(kBlack);
    h2_ -> SetLineWidth(2);
    h2_ -> SetMarkerStyle(20);
    TH1D* ratio_bkg = (TH1D*)h2->Clone("ratio");
    ratio_bkg->Divide(h2_);
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

    // Update the canvas
    canvas_bkg->Update();
    //_________________________________________________________________________________________

    /**
     * SUBTRACTION:
    */
    auto h_true_signal = antinu_hist.GetHistogram(FIDUCIAL_VOLUME, SELECTION_NONE, CCQE_ON_H, 0);
    auto h_plastic =  antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_POSITIVE_PLASTIC, REACTION_ANY, 1);
    auto h_graphite =  antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_GRAPHITE, CC_ON_CARBON, 1);
    auto h_residual = antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_PLASTIC_H, CCRES_ON_H, 1);
    
    auto h_selected_true_positive = antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_TRUE_POSITIVE, CCQE_ON_H, 0);

    h_plastic->Sumw2();
    h_graphite->Sumw2();
    h_true_signal->Sumw2();
    h_selected_true_positive->Sumw2();

    h_plastic->Add(h_graphite, -scale_factor);
    h_plastic->Add(h_residual, -1.);
    
    // h_true_signal->Scale(1. / h_true_signal->Integral());
    // h_selected_true_positive->Scale(1. / h_selected_true_positive->Integral());
    // h_plastic->Scale(1. / h_plastic->Integral());

    // STYLE
    // h_true_signal->GetYaxis()->SetTitle("Hist Normalized");
    h_true_signal->SetLineColor(kRed);
    h_plastic->SetLineColor(kBlack);
    h_plastic->SetMarkerStyle(20);
    h_plastic->SetMarkerSize(0.8);
    // h_true_signal->SetMaximum(0.13);
    h_plastic->SetMaximum(0.13);
    // STYLE

    TLegend* legendsub = new TLegend(0.6, 0.6, 0.9, 0.9);
    legendsub->AddEntry(h_true_signal, "#Phi_{#bar{#nu}_{#mu}} #sigma_{S} (CCQE on H)", "l");
    legendsub->AddEntry(h_selected_true_positive, "#Phi_{#bar{#nu}_{#mu}} #sigma_{S} #epsilon_{S} (CCQE on H like)", "lp");
    legendsub->AddEntry(h_plastic, "#Phi_{#bar{#nu}_{#mu}} #sigma_{S} #epsilon_{S} R_{det} (w/ bkg subtraction)", "lp");

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
    h_true_signal->Draw("E HIST");
    h_selected_true_positive->Draw("E HIST SAME");
    h_plastic->Draw("E1 SAME");
    hFlux -> Draw("H SAME");
    legendsub->Draw("SAME");

    pad2->cd();
    TH1D* ratio = (TH1D*)h_plastic->Clone("ratio");
    TH1D* ratio2 = (TH1D*)h_selected_true_positive->Clone("ratio2");
    
    ratio->Divide(h_true_signal);
    ratio->SetTitle(""); // Remove title for ratio plot
    ratio->GetYaxis()->SetTitle("Bin Content Ratio");
    ratio->GetXaxis()->SetTitle("Energy (GeV)");
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitleSize(0.1);
    ratio->GetYaxis()->SetTitleOffset(0.4);
    ratio->GetYaxis()->SetLabelSize(0.08);
    ratio->GetYaxis()->SetRangeUser(0, 0.3);
    ratio->GetXaxis()->SetTitleSize(0.1);
    ratio->GetXaxis()->SetTitleOffset(0.8);
    ratio->GetXaxis()->SetLabelSize(0.08);
    ratio->SetLineColor(kBlack);
    ratio->Draw("E1");

    ratio2->Divide(h_true_signal);
    ratio2->Draw("E HIST SAME");

    // Aggiorna il canvas
    canvas_sub->Update();

    /***
        EVENTRATE:
    */
    TCanvas* canvas_rate = new TCanvas("canvas_rate", "", 900, 700);
    TH1D* reco_signal_rate = (TH1D*)h_plastic->Clone("reco_signal_rate");
    reco_signal_rate->Sumw2();
    reco_signal_rate->Divide(ratio2);
    
    TPad* pad1rate = new TPad("pad1", "Main Plot", 0.0, 0.3, 1.0, 1.0); // Top pad
    TPad* pad2rate = new TPad("pad2", "Ratio Plot", 0.0, 0.0, 1.0, 0.3); // Bottom pad

    TLegend* legend_rate = new TLegend(0.6, 0.6, 0.9, 0.9);
    legend_rate->AddEntry(h_true_signal, "#Phi_{#bar{#nu}_{#mu}} #sigma_{S} (CCQE on H)", "l");
    legend_rate->AddEntry(reco_signal_rate, "reconstructed rate", "lp");

    pad1rate->SetBottomMargin(0.02); // Reduce bottom margin for top pad
    pad2rate->SetTopMargin(0.02);   // Reduce top margin for bottom pad
    pad2rate->SetBottomMargin(0.3); // Increase bottom margn for labels

    canvas_rate->cd();
    pad1rate->Draw();
    pad2rate->Draw();

    pad1rate->cd();
    h_true_signal->Draw("E HIST");
    reco_signal_rate->Draw("E1 SAME");
    legend_rate->Draw("SAME");

    pad2rate->cd();
    TH1D* ratio_rate_true_reco = (TH1D*)reco_signal_rate->Clone("ratio_rate_true_reco");
    ratio_rate_true_reco->Sumw2();
    ratio_rate_true_reco->Divide(h_true_signal);

    ratio_rate_true_reco->Draw("E1");
    canvas_rate->Update();


    /***
     * EFFICIENCIES: ECAL reconstruction efficiency on neutron as function of K neutron
    */

    // TCanvas* c3 = new TCanvas("eff1","",900,700);

    // neutrons_w_hit_in_ECAL->Draw("AP");

    // neutrons_w_compatible_cells_in_ECAL->SetLineColor(kGreen);
    // neutrons_w_complete_cells_in_ECAL->Draw("AP");
    
    // neutrons_w_compatible_cells_in_ECAL->SetLineColor(kRed);
    // neutrons_w_compatible_cells_in_ECAL->Draw("EP SAME"); 

    // TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Posizione: in alto a destra
    // legend->AddEntry(neutrons_w_hit_in_ECAL, "fraction of neutrons w hits in ECAL", "l"); // "l" per linea
    // legend->AddEntry(neutrons_w_complete_cells_in_ECAL, "Complete cells", "l"); // "l" per linea
    // legend->AddEntry(neutrons_w_compatible_cells_in_ECAL, "fraction of neutrons with ECAL coincidence", "l"); // "p" per marker
    // legend->SetTextSize(0.03); // Dimensione del testo
    // legend->Draw();

    // AS FUNCTION OF NEUTRINO ENERGY _____________________________________________

    // TCanvas* c = new TCanvas("eff2","",900,700);
    // neutrons_w_hit_in_ECAL_vs_neutrino_E -> Draw("AP");
    // neutrons_w_compatible_cells_in_ECAL_vs_neutrino_E->SetLineColor(kRed);
    // neutrons_w_compatible_cells_in_ECAL_vs_neutrino_E -> Draw("EP SAME");

    // TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Posizione: in alto a destra
    // legend->AddEntry(neutrons_w_hit_in_ECAL_vs_neutrino_E, "fraction of neutrons w hits in ECAL", "l"); // "l" per linea
    // legend->AddEntry(neutrons_w_compatible_cells_in_ECAL_vs_neutrino_E, "fraction of neutrons with ECAL coincidence", "l"); // "p" per marker
    // legend->SetTextSize(0.03); // Dimensione del testo
    // legend->Draw();

    // neutrons_w_compatible_cells_in_ECAL_vs_neutrino_E_->Draw("AP");

    /***
        CUT: optimize ECAL cut
    */
    // TCanvas* c = new TCanvas("","",900,700);
    // neutron_space_time_residulas -> Draw("COL Z");

    /***
        CELLS:
    */
    // TCanvas* c = new TCanvas("","",900,700);

    // // h_nof_cells_fired -> Draw("HIST");
    // h_nof_cells_fired_neutron -> SetLineColor(kRed);
    // h_nof_cells_fired_antimu -> SetLineColor(kBlack);
    // h_nof_cells_fired_neutron -> Draw("HIST");
    // h_nof_cells_fired_antimu -> Draw("HIST SAME");
    // std::cout << "total : " << h_nof_cells_fired_neutron->Integral() << "\n";


    // // 
    // TCanvas* c2 = new TCanvas("c2","",900,700);
    // adc_count_antimu -> Draw("HIST");
    // adc_count_neutron -> SetLineColor(kRed);
    // adc_count_neutron -> Draw("HIST SAME");

    // //
    // TCanvas* c3 = new TCanvas("c3","",900,700);
    // adc1_vs_adc2_incomplete->Draw("col z");

    // // nod events with 0 complete cells
    // std::cout << "nof of events with zero complete cells fired by neutrons : " << nof_events_w_zero_complete_cells << "\n";

    // // 
    // TCanvas* c4 = new TCanvas("c4","",900,700);
    // cell_reco_energy_neutron -> Draw("HIST");
    // cell_reco_energy_antimu -> Draw("HIST SAME");

    /***
    XSEC:
     */
    // TCanvas* cx = new TCanvas("c3","",900,700);
    // h_true_signal->Draw("HIST");
    // hFlux->Draw("HIST SAME");

    // TPad* pad1x = new TPad("pad1", "Main Plot", 0.0, 0.3, 1.0, 1.0); // Top pad
    // TPad* pad2x = new TPad("pad2", "Ratio Plot", 0.0, 0.0, 1.0, 0.3); // Bottom pad
    
    // pad1x->cd();
    // h_true_signal->Draw("HIST");
    // hFlux -> Draw("HIST SAME");

    // pad2x->cd();
    // TH1D* ratiox = (TH1D*)h_true_signal->Clone("ratio");
    // hFlux->Sumw2();
    // ratiox->Divide(hFlux);
    
    // ratiox->SetTitle(""); // Remove title for ratio plot
    // ratiox->GetYaxis()->SetTitle("Bin Content Ratio");
    // ratiox->GetXaxis()->SetTitle("Energy (GeV)");
    // ratiox->GetYaxis()->SetNdivisions(505);
    // ratiox->GetYaxis()->SetTitleSize(0.1);
    // ratiox->GetYaxis()->SetTitleOffset(0.4);
    // ratiox->GetYaxis()->SetLabelSize(0.08);
    // ratiox->GetYaxis()->SetRangeUser(0, 0.3);
    // ratiox->GetXaxis()->SetTitleSize(0.1);
    // ratiox->GetXaxis()->SetTitleOffset(0.8);
    // ratiox->GetXaxis()->SetLabelSize(0.08);
    // ratiox->SetLineColor(kBlack);
    // ratiox->Draw("E1");

    // cx->Update();
    
    // TCanvas* cx2 = new TCanvas("c3","",900,700);
    // ratiox->Draw("E");


    /***
    PLOT:
     */
    // TCanvas* canvas_residulas = new TCanvas("cr1","",900,700);
    // best_time_residual_vs_any_space_residual -> Draw("col z");
    
    // TCanvas* canvas_residulas_2 = new TCanvas("cr2","",900,700);
    // any_time_residual_vs_best_space_residual -> Draw("col z");

    // TCanvas* canvas_proj1 = new TCanvas("cp1","",900,700);
    // best_time_residuals -> Draw("HIST");

    // TCanvas* canvas_proj2 = new TCanvas("cp2","",900,700);
    // best_space_residuals -> Draw("HIST");

    // TCanvas* canvas_3cells = new TCanvas("canvas_3cells","",900,700);
    // res_time_vs_res_space_3cells->Draw("col z");

    // TCanvas* canvas_3cells_c = new TCanvas("canvas_3cells_","",900,700);
    // res_time_vs_res_space_3cells_->Draw("col z");

    /***
    RECAP: final recap on number of events
    */

    int total_ccqe_on_H = antinu_hist.GetHistogram(FIDUCIAL_VOLUME, SELECTION_NONE, CCQE_ON_H, 0)->GetEntries();
    int total_ccqe_on_H_selected = antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_TRUE_POSITIVE, CCQE_ON_H, 1)->GetEntries();
    
    int total_cc_on_Carbon_plastic = 0;
    int total_cc_on_Carbon_plastic_selected = antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_PLASTIC_C, CC_ON_CARBON, 1)->GetEntries();
    int total_cc_on_Carbon_graphite = 0;
    int total_cc_on_Carbon_graphite_selected = antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_GRAPHITE, CC_ON_CARBON, 1)->GetEntries();

    int total_ccres_on_h = 0;
    int total_ccres_on_h_selected = antinu_hist.GetHistogram(ECAL_COINCIDENCE, SELECTED_FALSE_POSITIVE_PLASTIC_H, CCRES_ON_H, 1)->GetEntries();

    std::cout << 
            "count CCQE on H : " << total_ccqe_on_H << 
                   ", selected " << total_ccqe_on_H_selected << "\n"
            
            "count CC on Carbon plastic : " << total_cc_on_Carbon_plastic <<
                              ", selected " << total_cc_on_Carbon_plastic_selected << "\n"
            
            "count CCRES on H plastic: " << total_ccres_on_h << 
                           ", selected " << total_ccres_on_h_selected << "\n"
            
            "count CC on Carbon graphite : " << total_cc_on_Carbon_graphite <<  
                               ", selected " << total_cc_on_Carbon_graphite_selected << "\n";
            

    return 0;
}