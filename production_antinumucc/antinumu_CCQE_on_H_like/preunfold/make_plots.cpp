#include <string>
#include <vector>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

// #include "RooUnfoldResponse.h"
// #include "RooUnfoldBayes.h"

struct PlotConfig {
    std::vector<TH1D*> histograms;
    std::string canvasTitle;
    std::vector<std::string> legendLabels;
    std::string pdf_name;
    bool normalize;
    bool draw_ratio;
};

double z_c_target[8] = {
    25.4422, // X0
    25.0776, // X1
    24.7130, // C0
    24.3484, // B0
    23.9838, // A0
    23.6192, // A1
    23.2546, // B1
    22.8900, // C1
};

bool is_close_2_c_target(const double z, double epsilon = 0.04){ // 4 cm
    for (double value : z_c_target) {
        if (std::fabs(value - z) < epsilon) {
            return true; // Trovato un valore vicino
        }
    }
    return false; // Nessun valore vicino trovato
}

enum Stage{
    FIDUCIAL_VOLUME,
    WIRES_CUT,
    CHARGE_MULTIPLICITY,
    ECAL_COINCIDENCE
};

enum EventType{
    SELECTED_TRUE_POSITIVE,
    FALSE_NEGATIVE,
    // FALSE POSITIVES
    SELECTED_FALSE_POSITIVE_GRAPHITE,
    SELECTED_FALSE_POSITIVE_PLASTIC_C,
    SELECTED_FALSE_POSITIVE_PLASTIC_H,
    SELECTED_FALSE_POSITIVE_OTHER,
    // TRUE NEGATIVES
    TRUE_NEGATIVE_GRAPHITE,
    TRUE_NEGATIVE_PLASTIC_C,
    TRUE_NEGATIVE_PLASTIC_H,
    TRUE_NEGATIVE_OTHERS,
    DEFAULT
};

void CondifgurePad(TVirtualPad& c, std::string title){
    c.DrawFrame(0,0,8,0.3,title.c_str());
}

EventType determineEventType(bool is_event_selected, bool is_signal, bool is_in_graphite, bool is_in_plastic, bool is_on_C, bool is_on_H){
    if(is_event_selected)
    {
        if(is_signal)
        {
            return SELECTED_TRUE_POSITIVE;
        }else
        {
            if(is_in_graphite)
            {
                return SELECTED_FALSE_POSITIVE_GRAPHITE;
            }else if (is_in_plastic)
            {
                if(is_on_C)
                {
                    return SELECTED_FALSE_POSITIVE_PLASTIC_C;
                }else if(is_on_H)
                {
                    return SELECTED_FALSE_POSITIVE_PLASTIC_H;
                }
            }else
            {
                return SELECTED_FALSE_POSITIVE_OTHER;
            }
        }
    }else
    {
        if(is_signal)
        {
            return FALSE_NEGATIVE;
        }else
        {
            if(is_in_graphite)
            {
                return TRUE_NEGATIVE_GRAPHITE;
            }else if(is_in_plastic)
            {
                if(is_on_C)
                {
                    return TRUE_NEGATIVE_PLASTIC_C;
                }else if(is_on_H)
                {
                    return TRUE_NEGATIVE_PLASTIC_H;
                }
            }else
            {
                return TRUE_NEGATIVE_OTHERS;
            }
        }
    }
    return DEFAULT;
}
void plotHistograms(const std::vector<TH1D*>& histograms, 
                    const std::string& canvasTitle, 
                    const std::vector<std::string>& legendLabels, 
                    bool normalize = false, 
                    bool draw_ratio = false,
                    std::string pdf_name = "") {
    
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    std::vector<TH1D*> clones(histograms.size());

    if (histograms.empty()) {
        std::cerr << "Errore: nessun istogramma fornito." << std::endl;
        return;
    }

    double maxY = 0;
    
    for (size_t i = 0; i < histograms.size(); ++i) 
    {

        TH1D* hist = static_cast<TH1D*>(histograms[i]->Clone());
        
        if (normalize && hist->Integral() > 0) hist->Scale(1.0 / hist->Integral());

        double maxVal = (hist->Integral() > 0) ? hist->GetMaximum() : 0.;

        if (maxVal > maxY) maxY = maxVal;

        clones[i] = hist;
        
    }

    maxY *= 1.1; 
    
    TCanvas* canvas = new TCanvas("canvas", "", 900, 700);

    TLegend* legend = new TLegend(0.6, 0.6, 0.9, 0.9);

    std::vector<int> colors = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kYellow+1};
    
    if (draw_ratio && clones.size() == 2) {
        auto rp = new TRatioPlot(clones[0], clones[1]);
        rp->Draw();
        rp->GetLowerRefYaxis()->SetTitle("hists ratio");
        rp->GetLowerRefGraph()->SetMinimum(0.0);
        rp->GetLowerRefGraph()->SetMaximum(1.5);
        rp->GetUpperPad()->cd();
    }
    canvas->Update();
    
    for (size_t i = 0; i < clones.size(); ++i)
    // for (size_t i = 0; i < 2; ++i)
    {
        
        if (i == 0){
            clones[i]->SetTitle(canvasTitle.c_str());
        }else{
            clones[i]->SetTitle("");
            clones[i]->SetName("");
        }
        legend->AddEntry(clones[i], legendLabels[i].c_str(), "l");
        int color = colors[i % colors.size()];
        clones[i]->SetMaximum(maxY);
        clones[i]->SetLineColor(color);
        clones[i]->SetLineWidth(2);
        clones[i]->Draw((i == 0) ? "E1 HIST" : "E1 HIST SAME");
        // clones[i]->Draw((i == 0) ? "E1" : "E1 SAME");
    }


    legend->Draw();
    canvas->Update();

    if (!pdf_name.empty()) {
        canvas->SaveAs(pdf_name.c_str());
    }
}

void plotEfficiencies(const std::vector<TEfficiency*>& efficiencies, const std::vector<std::string>& labels) {

    TCanvas* canvas = new TCanvas("canvas", "Efficienza", 800, 600);

    TLegend* legend = new TLegend(0.8, 0.8, 0.9, 0.9);
    legend->SetTextSize(0.03);  // Dimensione del testo nella legenda

    std::vector<int> colors = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kYellow+1};
    
    for (size_t i = 0; i < efficiencies.size(); ++i) {
        // Assegna un colore alla curva
        int color = colors[i % colors.size()];
        efficiencies[i]->SetLineColor(color);
        efficiencies[i]->SetMarkerColor(color);

        // Disegna l'oggetto TEfficiency nella canvas
        if (i == 0) {
            efficiencies[i]->Draw("AP");  // Disegna il primo oggetto con assi
        } else {
            efficiencies[i]->Draw("P SAME");  // Disegna gli altri sopra
        }

        legend->AddEntry(efficiencies[i], labels[i].c_str(), "lep");
    }
    
    legend->Draw();
    canvas->Update();
}

int make_plots(){
  
  TFile *file = TFile::Open("/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/antinumu_CCQE_on_H_like/preunfold/merged_productions.0000.to.0059.root");

  int CCQEonHydrogen;
  int NofFinalStateChargedParticles;
  bool pass_nof_wires_cut;
  int candidate_signal_event; // 
  std::string* InteractionTarget = new std::string();
  std::string* InteractionVolume_short = new std::string();
  double IncomingNeutrino_energy;
  double Neutrino_reconstructed_energy_GeV;
  double Interaction_vtxX;
  double Interaction_vtxY;
  double Interaction_vtxZ;
  TVector3* Antimuon_p_true = nullptr;; 
  TLorentzVector* Antimuon_reconstructed_P4 = nullptr;; 
  TLorentzVector* FinalStateHadronicSystemTotal4Momentum = nullptr; 
  TVector3* PredictedNeutron_P3_GeV = nullptr;; 
  double PredictedNeutron_E_GeV; 

  TTree* tree = (TTree*)file->Get("preunfold");

  tree->SetBranchAddress("CCQEonHydrogen", &CCQEonHydrogen);
  tree->SetBranchAddress("NofFinalStateChargedParticles", &NofFinalStateChargedParticles);
  tree->SetBranchAddress("pass_nof_wires_cut", &pass_nof_wires_cut);
  tree->SetBranchAddress("candidate_signal_event", &candidate_signal_event);
  tree->SetBranchAddress("InteractionTarget", &InteractionTarget);
  tree->SetBranchAddress("InteractionVolume_short", &InteractionVolume_short);
  tree->SetBranchAddress("IncomingNeutrino_energy", &IncomingNeutrino_energy);
  tree->SetBranchAddress("Neutrino_reconstructed_energy_GeV", &Neutrino_reconstructed_energy_GeV);
  tree->SetBranchAddress("Interaction_vtxX", &Interaction_vtxX);
  tree->SetBranchAddress("Interaction_vtxY", &Interaction_vtxY);
  tree->SetBranchAddress("Interaction_vtxZ", &Interaction_vtxZ);
  tree->SetBranchAddress("Interaction_vtxZ", &Interaction_vtxZ);
  tree->SetBranchAddress("FinalStateHadronicSystemTotal4Momentum", &FinalStateHadronicSystemTotal4Momentum);
  tree->SetBranchAddress("Antimuon_p_true", &Antimuon_p_true); // MeV
  tree->SetBranchAddress("Antimuon_reconstructed_P4", &Antimuon_reconstructed_P4); // MeV

  TH1D* h_true = new TH1D ("true", "Neutrino energy; [GeV]", 16, 0., 8.);

  auto nof_entries = tree->GetEntries();

  std::cout << "Number of events in fiducial volume " << nof_entries << "\n";

  /*** ALL HIST HAVE FIDUCIAL VOLUME CUT
  TRUE COMPONENTS_________________________________________________________________________________________
  @h_true_signal: true antineutrino spectrum for signal rate
  @h_true_background: true antineutrino spectrum for background rate
  @h_true_c_background_plastic: true antineutrino spectrum for C background rate in graphite
  @h_true_c_background_graphite: true antineutrino spectrum for C background rate in plastic
  @h_true_plastic: true antineutrino spectrum for all interactions in plastic
  @h_true_graphite: true antineutrino spectrum for  all interactions in graphite

  CCQE-LIKE SELECTED______________________________________________________________________________________
  @h_reco_signal_like: reconstructed antineutrino spectrum for selected events
  @h_reco_signal_like_true_positive: reconstructed antineutrino spectrum for selected true positive events
  @h_reco_signal_like_plastic: reconstructed antineutrino spectrum for selected events in plastic
  @h_reco_signal_like_graphite: reconstructed antineutrino spectrum for selected events in graphite
  @h_reco_signal_like_plastic_graphite: reconstructed antineutrino spectrum for selected events in plastic with statistical subtraction of graphite background
  @h_reco_signal_like_on_c_plastic: reconstructed antineutrino spectrum for selected events in plastic interaction on Carbon nuclei

  ANTIMUON SPECTRA__________________________________________________________________________________________
  @h_true_antimu_c_background_plastic: true antimuon spectrum for C background rate in graphite
  @h_true_antimu_c_background_graphite: true antimuon spectrum for C background rate in plastic
  @h_reco_antimu_c_background_plastic: reco antimuon spectrum for C background rate in graphite
  @h_reco_antimu_c_background_graphite: reco antimuon spectrum for C background rate in plastic
   */
  uint nof_bins = 20;
  double xlow = 0.;
  double xup = 8.;
    // TRUE COMPONENTS
  TH1D* h_true_signal = new TH1D ("h_true_signal", "", nof_bins, xlow, xup);
  TH1D* h_true_background = new TH1D ("h_true_background", "", nof_bins, xlow, xup);
  TH1D* h_true_C_background_plastic = new TH1D ("h_true_C_background_plastic", "", nof_bins, xlow, xup);
  TH1D* h_true_C_background_graphite = new TH1D ("h_true_C_background_graphite", "", nof_bins, xlow, xup);
  TH1D* h_true_plastic = new TH1D ("h_true_plastic", "", nof_bins, xlow, xup);
  TH1D* h_true_graphite = new TH1D ("h_true_graphite", "", nof_bins, xlow, xup);
   // CCQE - LIKE SELECTED
  TH1D* h_reco_signal_like = new TH1D ("h_reco_signal_like", "", nof_bins, xlow, xup);
  TH1D* h_reco_signal_like_true_positive = new TH1D ("h_reco_signal_like_true_positive", "", nof_bins, xlow, xup);
  TH1D* h_reco_signal_like_plastic = new TH1D ("h_reco_signal_like_plastic", "", nof_bins, xlow, xup);
  TH1D* h_reco_signal_like_graphite = new TH1D ("h_reco_signal_like_graphite", "", nof_bins, xlow, xup);
  TH1D* h_reco_signal_like_plastic_graphite = new TH1D ("h_reco_signal_like_plastic_graphite", "", nof_bins, xlow, xup);
  TH1D* h_reco_signal_like_on_C_plastic = new TH1D ("h_reco_signal_like_on_C_plastic", "", nof_bins, xlow, xup);
  TH1D* h_reco_signal_like_on_C_plastic_1st_after_gra = new TH1D ("h_reco_signal_like_on_C_plastic_1st_after_gra", "", nof_bins, xlow, xup);
  // ANTIMUON SPECTRA
  TH1D* h_true_antimu_C_background_plastic = new TH1D ("h_true_antimu_C_background_plastic", "", nof_bins, xlow, xup);
  TH1D* h_true_antimu_C_background_graphite = new TH1D ("h_true_antimu_C_background_graphite", "", nof_bins, xlow*1e3, xup*1e3);
  TH1D* h_reco_antimu_C_background_plastic = new TH1D ("h_reco_antimu_C_background_plastic", "", nof_bins, xlow*1e3, xup*1e3);
  TH1D* h_reco_antimu_C_background_graphite = new TH1D ("h_reco_antimu_C_background_graphite", "", nof_bins, xlow*1e3, xup*1e3);
  /***
   EFFICIENCIES:
   @selection_eff_on_graphite: selection efficiency of carbon in graphite
   @selection_eff_on_plastic: selection efficiency of carbon in plastic
   */
   TEfficiency* sel_eff_graphite = new TEfficiency("sel_eff_graphite","CCQE-like selection efficiency;E [Gev];#epsilon",10,0,5);
   TEfficiency* sel_eff_carbon_in_plastic = new TEfficiency("eff","CCQE-like selection efficiency;E [Gev];#epsilon",10,0,5);
   TEfficiency* sel_eff_H_in_plastic = new TEfficiency("eff","Selection efficiency events in plastic; E_{#bar_{#nu}_{#mu}} [Gev];#epsilon",5,0,5);

  /***
   STATISTICS:
   */
   int total_in_fv = 0;
   int total_signal = 0;

   int total_in_graphite = 0;
   int total_in_plastic_on_C = 0;
   int total_in_plastic_on_H = 0;

   int total_selected_in_plastic_C = 0;
   int total_selected_in_plastic_H = 0;
   int total_selected_in_graphite = 0;

   int this_cout =0;

   /***
   UNFOLDING:
    */
    // RooUnfoldResponse response(16, 0., 8.);
    uint train_sample_size = 5000000;
    TH1D* h_reco_test_sample_signal_like = new TH1D ("h_reco_test_sample_signal_like", "", nof_bins, xlow, xup);
    TH1D* h_true_test_sample_signal_like = new TH1D ("h_true_test_sample_signal_like", "", nof_bins, xlow, xup);

  for (size_t i = 0; i < nof_entries; i++)
  {
    
    tree->GetEntry(i);

    bool is_signal = (CCQEonHydrogen==1);
    bool is_event_selected = (NofFinalStateChargedParticles==1 && pass_nof_wires_cut==1 && candidate_signal_event==1);
    bool is_in_graphite = (*InteractionVolume_short == "C_Target");
    bool is_in_plastic = (*InteractionVolume_short == "C3H6_Target");
    bool is_on_H = (*InteractionTarget == "proton");
    bool is_on_C = (*InteractionTarget == "C12");

    EventType eventType = determineEventType(is_event_selected, is_signal, is_in_graphite, is_in_plastic, is_on_C, is_on_H);
    
    total_in_fv++;

    double antimuon_true_energy = sqrt(Antimuon_p_true->Mag() * Antimuon_p_true->Mag() +  105.658 * 105.658);

    
    switch (eventType)
    {
        case SELECTED_TRUE_POSITIVE:
            total_signal++;
            total_in_plastic_on_H++;
            h_true_signal -> Fill(IncomingNeutrino_energy);
            h_reco_signal_like_true_positive -> Fill(Neutrino_reconstructed_energy_GeV);
            h_reco_signal_like -> Fill(Neutrino_reconstructed_energy_GeV);
            h_reco_signal_like_plastic -> Fill(Neutrino_reconstructed_energy_GeV);
            sel_eff_H_in_plastic -> Fill(1, IncomingNeutrino_energy);
            if(i < train_sample_size)
            { // construct ufolding matrix
                // response.Fill(Neutrino_reconstructed_energy_GeV, IncomingNeutrino_energy);
            }else
            { // fill test sample
                h_true_test_sample_signal_like->Fill(IncomingNeutrino_energy);
                h_reco_test_sample_signal_like->Fill(Neutrino_reconstructed_energy_GeV);
            }
            break;
        
        case FALSE_NEGATIVE:
            total_signal++;
            total_in_plastic_on_H++;
            h_true_signal -> Fill(IncomingNeutrino_energy);
            sel_eff_H_in_plastic -> Fill(0, IncomingNeutrino_energy);
            if(i < train_sample_size)
            {// construct ufolding matrix
                // response.Miss(IncomingNeutrino_energy);
            } else
            { // fill test sample
                h_true_test_sample_signal_like->Fill(IncomingNeutrino_energy);
            }
            break;

        case SELECTED_FALSE_POSITIVE_GRAPHITE:
            total_in_graphite++;
            total_selected_in_graphite++;
            h_true_background -> Fill(IncomingNeutrino_energy);
            h_true_C_background_graphite -> Fill(IncomingNeutrino_energy);
            h_true_antimu_C_background_graphite -> Fill(antimuon_true_energy);
            h_reco_antimu_C_background_graphite -> Fill(Antimuon_reconstructed_P4->T());
            h_true_graphite -> Fill(IncomingNeutrino_energy);
            h_reco_signal_like -> Fill(Neutrino_reconstructed_energy_GeV);
            h_reco_signal_like_graphite -> Fill(Neutrino_reconstructed_energy_GeV);
            sel_eff_graphite -> Fill(1, IncomingNeutrino_energy);
            break;
        
        case SELECTED_FALSE_POSITIVE_PLASTIC_C:
            total_in_plastic_on_C++;
            total_selected_in_plastic_C++;
            h_true_background -> Fill(IncomingNeutrino_energy);
            h_true_C_background_plastic -> Fill(IncomingNeutrino_energy);
            h_true_antimu_C_background_plastic -> Fill(antimuon_true_energy);
            h_reco_antimu_C_background_plastic -> Fill(Antimuon_reconstructed_P4->T());
            h_true_plastic -> Fill(IncomingNeutrino_energy);
            h_reco_signal_like -> Fill(Neutrino_reconstructed_energy_GeV);
            h_reco_signal_like_plastic -> Fill(Neutrino_reconstructed_energy_GeV);
            h_reco_signal_like_on_C_plastic -> Fill(Neutrino_reconstructed_energy_GeV);
            // std::cout << "Interaction_vtxZ " << Interaction_vtxZ <<"\n";
            // if(is_close_2_c_target(Interaction_vtxZ)){
            //     h_reco_signal_like_on_C_plastic_1st_after_gra -> Fill(Neutrino_reconstructed_energy_GeV);
            //     this_cout++;
            // }
            sel_eff_carbon_in_plastic -> Fill(1, IncomingNeutrino_energy);
            break;
        
        case SELECTED_FALSE_POSITIVE_PLASTIC_H:
            h_true_background -> Fill(IncomingNeutrino_energy);
            h_true_plastic -> Fill(IncomingNeutrino_energy);
            h_reco_signal_like -> Fill(Neutrino_reconstructed_energy_GeV);
            h_reco_signal_like_plastic -> Fill(Neutrino_reconstructed_energy_GeV);
            break;
        
        case SELECTED_FALSE_POSITIVE_OTHER:
            h_true_background -> Fill(IncomingNeutrino_energy);
            h_reco_signal_like -> Fill(Neutrino_reconstructed_energy_GeV);
            break;
        
        case TRUE_NEGATIVE_GRAPHITE:
            total_in_graphite++;
            h_true_background -> Fill(IncomingNeutrino_energy);
            h_true_C_background_graphite -> Fill(IncomingNeutrino_energy);
            h_true_antimu_C_background_graphite -> Fill(antimuon_true_energy);
            h_true_graphite -> Fill(IncomingNeutrino_energy);
            sel_eff_graphite -> Fill(0, IncomingNeutrino_energy);
            break;
        
        case TRUE_NEGATIVE_PLASTIC_C:
            total_in_plastic_on_C++;
            h_true_background -> Fill(IncomingNeutrino_energy);
            h_true_plastic -> Fill(IncomingNeutrino_energy);
            h_true_C_background_plastic -> Fill(IncomingNeutrino_energy);
            h_true_antimu_C_background_plastic -> Fill(antimuon_true_energy);
            sel_eff_carbon_in_plastic -> Fill(0, IncomingNeutrino_energy);

            break;
        
        case TRUE_NEGATIVE_PLASTIC_H:
            total_in_plastic_on_H++;
            h_true_background -> Fill(IncomingNeutrino_energy);
            h_true_plastic -> Fill(IncomingNeutrino_energy);
            break;
        
        case TRUE_NEGATIVE_OTHERS:
            h_true_background -> Fill(IncomingNeutrino_energy);
            break;
        
        default :
            break;
    }
    
  }
  std::vector<PlotConfig> plots = {
    /**
    PLOT:
        @variable: antineutrino Energy spectrum (shape)
        @stage_a: MC truth
        @hist_a: intearaction on Carbon nucleus of Plastic
        @stage_b: MC truth
        @hist_b: intearaction on Carbon nucleus of Graphite
     */
    {
        {h_true_C_background_plastic, h_true_C_background_graphite},
        " ; E_{#bar{#nu}_{#mu}}^{True} [GeV]; density",
        {"#bar{#nu}_{#mu} on Carbon (Plastic)", "#bar{#nu}_{#mu} on Carbon (Graphite)"},
        "plots/TRUE_POSIVITE_C_GRAPHITE_Etrue.TRUE_POSIVITE_C_PLASTIC_Etrue.pdf",
        true, false
    },
    
    /**
    PLOT:
        @variable: antineutrino Energy spectrum (shape)
        @stage_a: Selected CCQE-like
        @hist_a: intearaction in Plastic
        @stage_b: Selected CCQE-like
        @hist_b: intearaction on Carbon nucleus of Plastic
     */
    {
        {h_reco_signal_like_on_C_plastic, h_reco_signal_like_graphite},
        " ; E_{#bar{#nu}_{#mu}}^{Reco} [GeV]; density",
        {"#bar{#nu}_{#mu} on Carbon CCQE-like (Plastic)", "#bar{#nu}_{#mu} on Carbon CCQE-like (Graphite)"},
        "plots/CCQELIKE_C_GRAPHITE_Ereco.CCQELIKE_C_PLASTIC_Ereco.pdf",
        true, false
    },
    
    /**
    PLOT:
        @variable: antineutrino Energy spectrum (shape)
        @stage_a: MC truth
        @hist_a: interaction on any nuclus and any target
        @stage_b: Selected CCQE-like
        @hist_b: interaction on any nuclus of plastic
     */
    {
        {h_true_signal, h_reco_signal_like_plastic},
        " ; E_{#bar{#nu}_{#mu}}^{} [GeV]; density",
        {"True #bar{#nu}_{#mu} CCQE on H", "#bar{#nu}_{#mu} CCQE-like (Plastic)"},
        "plots/TRUE_POSIVITE_Etrue.SEL_CCQELIKE_PLASTIC_Ereco.pdf",
        true, false
    },

    /**
    PLOT:
        @variable: antimuon Energy spectrum (shape)
        @stage_a: MC truth
        @hist_a: intearaction on Carbon nucleus of Plastic
        @stage_b: MC truth
        @hist_b: intearaction on Carbon nucleus of Graphite
     */

    {
        {h_true_antimu_C_background_plastic, h_true_antimu_C_background_graphite},
        "; E_{#mu^+}^{True} [GeV]; density",
        {"#bar{#nu}_{#mu} on Carbon (Plastic)", "#bar{#nu}_{#mu} on Carbon (Graphite)"},
        "plots/antimu.CCQELIKE_C_GRAPHITE_Etrue.CCQELIKE_C_PLASTIC_Etrue.pdf",
        true, false
    },

    /**
    PLOT:
        @variable: antimuon Energy spectrum (shape)
        @stage_a: Selected CCQE-like
        @hist_a: intearaction on Carbon nucleus of Plastic
        @stage_b: Selected CCQE-like
        @hist_b: intearaction on Carbon nucleus of Graphite
     */

    {
        {h_reco_antimu_C_background_plastic, h_reco_antimu_C_background_graphite},
        "; E_{#mu^+}^{True} [GeV]; density",
        {"#bar{#nu}_{#mu} on Carbon CCQE-like (Plastic)", "#bar{#nu}_{#mu} on Carbon CCQE-like (Graphite)"},
        "plots/antimu.CCQELIKE_C_GRAPHITE_Ereco.CCQELIKE_C_PLASTIC_Ereco.pdf",
        true, false
    },
    
    };
   for (const auto& plot : plots) {
       plotHistograms(plot.histograms, plot.canvasTitle, plot.legendLabels, 
                      plot.normalize, plot.draw_ratio, plot.pdf_name);
   }

   return 0;
}







    /**
    HYDROGEN: statistical subtratcion of interactions in graphite from interactions in palstic
     */

    // const double scale_factor = 4.084;

    // h_reco_signal_like_plastic_graphite = h_reco_signal_like_plastic;
    // h_reco_signal_like_plastic_graphite->Sumw2();
    // h_reco_signal_like_plastic_graphite->Add(h_reco_signal_like_graphite, - scale_factor);

    // legends = {"E_{ #bar{#nu}_{#mu}}^{True} True CCQE on H", "E_{ #bar{#nu}_{#mu}}^{Reco} Selected CCQE-like on H (subraction)"};

    // plotHistograms({h_true_signal, h_reco_signal_like_plastic_graphite}, // histots
    //                 "#bar{#nu}_{#mu} Energy; E [GeV]; density", // canvas name 
    //                 legends, // hist labels
    //                 true, // normalize
    //                 true, // show_errors
    //                 true, // ratio_plot
    //                 "plots/TRUE_POSIVITE_Etrue.SEL_CCQELIKE_H_SUBTRACTION_Ereco.pdf" // pdf_name
    //                 );



  /***
  STATISTICS:
   */
    //    std::cout << "number of interactions in graphite : " << h_true_graphite -> GetEntries() << "\n";
    //    std::cout << "number of interactions in plastic : " << h_true_plastic -> GetEntries() << "\n";

    //   /***
    //   EFFICIENCIES:
    //    */
    //   std::vector<TEfficiency*> effs = {};
    //   std::vector<std::string> labels = {"Carbon in Plastic", "Carbon in Graphite",};
    //   plotEfficiencies(effs, labels);
    
    //   std::vector<std::string> labels = {"Plastic"};
    //   std::vector<TEfficiency*> effs = {sel_eff_H_in_plastic};
    //   plotEfficiencies(effs, labels);

    /***
    UNFOLDING:
     */
    //   auto* R = response.HresponseNoOverflow();
    //   auto* c1 = new TCanvas();
    //   R->SetStats(0);
    //   R->Draw("colz");
    //   c1->Draw();
    //   c1->SaveAs("response.png"); 






// TGeoManager* geo = TGeoManager::Import("/storage/gpfs_data/neutrino/users/gi/dunendggd/SAND_opt3_DRIFT1.gdml");
// TString path2tracker("/volWorld_1/rockBox_lv_0/volDetEnclosure_0/volSAND_0/MagIntVol_volume_0/sand_inner_volume_0/SANDtracker_0");

// void FindCarbonTargetZ(){
//     auto path2CTarget_X0 = TString::Format("%s/SuperMod_X0_0/CMod_X0_0/CTarget_X0_0", path2tracker.Data());
//     auto path2CTarget_X1 = TString::Format("%s/SuperMod_X1_0/CMod_X1_0/CTarget_X1_0", path2tracker.Data());
//     auto path2CTarget_A0 = TString::Format("%s/SuperMod_A_0/CMod_A_0/CTarget_A_0", path2tracker.Data());
//     auto path2CTarget_B0 = TString::Format("%s/SuperMod_B_0/CMod_B_0/CTarget_B_0", path2tracker.Data());
//     auto path2CTarget_C0 = TString::Format("%s/SuperMod_C_0/CMod_C_0/CTarget_C_0", path2tracker.Data());

//     std::vector<TString> paths = {
//         path2CTarget_X0, 
//         path2CTarget_X1, 
//         path2CTarget_A0, 
//         path2CTarget_B0,
//         path2CTarget_C0};

//     for(auto& path : paths){

//         if (gGeoManager->cd(path.Data())) {
//                 // Ottieni il nodo attuale
//                 TGeoNode *node = gGeoManager->GetCurrentNode();

//                 if (node) {
//                     // Ottieni la matrice globale per ottenere le coordinate globali
//                     const TGeoMatrix *matrix = gGeoManager->GetCurrentMatrix();

//                     // Estrai le coordinate di traslazione globale del volume
//                     Double_t globalCoords[3] = {0, 0, 0};
//                     matrix->LocalToMaster(globalCoords, globalCoords);

//                     std::cout << "Coordinate globali del volume " << ": ("
//                               << globalCoords[0] << ", " << globalCoords[1] << ", " << globalCoords[2] << ")" << std::endl;
//                 } else {
//                     std::cerr << "Nodo non trovato!" << std::endl;
//                 }
//             } else {
//                 std::cerr << "Percorso non valido!" << std::endl;
//             }

//     }
// }

