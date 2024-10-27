#include <string>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <vector>

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
    TRUE_NEGATIVE_OTHERS
};

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
}

void plotHistograms(const std::vector<TH1D*>& histograms, 
                    const std::string& canvasTitle, 
                    const std::vector<std::string>& legendLabels, 
                    bool normalize = false, 
                    bool showErrorBars = false,
                    bool draw_ratio = true) {
    gStyle->SetOptStat(0);

    // Controllo per il numero di istogrammi
    if (histograms.empty()) {
        std::cerr << "Errore: nessun istogramma fornito." << std::endl;
        return;
    }

    // Normalizzazione opzionale degli istogrammi
    if (normalize) {
        for (auto* hist : histograms) {
            if (hist->Integral() > 0) hist->Scale(1.0 / hist->Integral());
        }
    }

    // Se si desidera mostrare le barre di errore, abilita la sommatoria dei quadrati per i pesi
    if (showErrorBars) {
        for (auto* hist : histograms) {
            hist->Sumw2();
        }
    }

    // Calcolo del massimo tra tutti gli istogrammi
    double maxY = 0;
    for (auto* hist : histograms) {
        double maxVal = hist->GetMaximum();
        if (maxVal > maxY) {
            maxY = maxVal;
        }
    }
    maxY *= 1.1;  // Aggiungi un margine del 10%

    // Creazione della canvas
    TCanvas* canvas = new TCanvas("canvas", canvasTitle.c_str(), 900, 700);
    canvas->SetTitle(canvasTitle.c_str());

    // Impostazioni di stile e disegno degli istogrammi
    std::vector<int> colors = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kYellow+1};
    for (size_t i = 0; i < histograms.size(); ++i) {
        TH1D* hist = histograms[i];
        int color = colors[i % colors.size()];
        hist->SetMaximum(maxY);  // Imposta lo stesso massimo su tutti gli istogrammi
        hist->SetLineColor(color);  // Assegna un colore diverso per ogni istogramma
        hist->SetLineWidth(2);

        if (showErrorBars) {
            hist->Draw((i == 0) ? "E HIST" : "E HIST SAME");
        } else {
            hist->Draw((i == 0) ? "HIST" : "HIST SAME");
        }
    }

    // Creazione della legenda
    // TLegend* legend = new TLegend(0.5, 0.5, 0.9, 0.9);
    TLegend* legend = new TLegend();
    for (size_t i = 0; i < histograms.size(); ++i) {
        legend->AddEntry(histograms[i], legendLabels[i].c_str(), "l");
    }

    if(draw_ratio && histograms.size()==2){
        auto rp = new TRatioPlot(histograms[0], histograms[1]);
        rp->Draw();
        rp->GetLowerRefGraph()->SetMinimum(0.);
        rp->GetLowerRefGraph()->SetMaximum(1.5);
    }

    legend->Draw();
    // canvas->Update();
    canvas->Draw();
}

void plotEfficiencies(const std::vector<TEfficiency*>& efficiencies, const std::vector<std::string>& labels) {

    TCanvas* canvas = new TCanvas("canvas", "Efficienza", 800, 600);

    TLegend* legend = new TLegend(0.6, 0.6, 0.9, 0.9);
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
  
  TFile *file = TFile::Open("/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/antinumu_CCQE_on_H_like/preunfold/merged_productions.0000.to.0050.root");

  int CCQEonHydrogen;
  int NofFinalStateChargedParticles;
  bool pass_nof_wires_cut;
  int candidate_signal_event; // 
  std::string* InteractionTarget = new std::string();
  std::string* InteractionVolume_short = new std::string();
  double IncomingNeutrino_energy;
  double Neutrino_reconstructed_energy_GeV;

  TTree* tree = (TTree*)file->Get("preunfold");

  tree->SetBranchAddress("CCQEonHydrogen", &CCQEonHydrogen);
  tree->SetBranchAddress("NofFinalStateChargedParticles", &NofFinalStateChargedParticles);
  tree->SetBranchAddress("pass_nof_wires_cut", &pass_nof_wires_cut);
  tree->SetBranchAddress("candidate_signal_event", &candidate_signal_event);
  tree->SetBranchAddress("InteractionTarget", &InteractionTarget);
  tree->SetBranchAddress("InteractionVolume_short", &InteractionVolume_short);
  tree->SetBranchAddress("IncomingNeutrino_energy", &IncomingNeutrino_energy);
  tree->SetBranchAddress("Neutrino_reconstructed_energy_GeV", &Neutrino_reconstructed_energy_GeV);

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
  @h_reco_singal_like_plastic_graphite: reconstructed antineutrino spectrum for selected events in plastic with statistical subtraction of graphite background
  @h_reco_signal_like_on_c_plastic: reconstructed antineutrino spectrum for selected events in plastic interaction on Carbon nuclei
   */
  uint nof_bins = 20;
  double xlow = 0.;
  double xup = 8.;
    // TRUE COMPONENTS
  TH1D* h_true_signal = new TH1D ("h_true_signal", "; [GeV]", nof_bins, xlow, xup);
  TH1D* h_true_background = new TH1D ("h_true_background", "; [GeV]", nof_bins, xlow, xup);
  TH1D* h_true_C_background_plastic = new TH1D ("h_true_C_background_plastic", "; [GeV]", nof_bins, xlow, xup);
  TH1D* h_true_C_background_graphite = new TH1D ("h_true_C_background_graphite", "; [GeV]", nof_bins, xlow, xup);
  TH1D* h_true_plastic = new TH1D ("h_true_plastic", "; [GeV]", nof_bins, xlow, xup);
  TH1D* h_true_graphite = new TH1D ("h_true_graphite", "; [GeV]", nof_bins, xlow, xup);
   // CCQE - LIKE SELECTED
  TH1D* h_reco_signal_like = new TH1D ("h_reco_signal_like", "; [GeV]", nof_bins, xlow, xup);
  TH1D* h_reco_signal_like_true_positive = new TH1D ("h_reco_signal_like_true_positive", "; [GeV]", nof_bins, xlow, xup);
  TH1D* h_reco_signal_like_plastic = new TH1D ("h_reco_signal_like_plastic", "; [GeV]; [GeV]", nof_bins, xlow, xup);
  TH1D* h_reco_signal_like_graphite = new TH1D ("h_reco_signal_like_graphite", "; [GeV]", nof_bins, xlow, xup);
  TH1D* h_reco_singal_like_plastic_graphite = new TH1D ("h_reco_singal_like_plastic_graphite", "; [GeV]", nof_bins, xlow, xup);
  TH1D* h_reco_signal_like_on_C_plastic = new TH1D ("h_reco_signal_like_on_C_plastic", "; [GeV]", nof_bins, xlow, xup);

  /***
   EFFICIENCIES:
   @selection_eff_on_graphite: selection efficiency of carbon in graphite
   @selection_eff_on_plastic: selection efficiency of carbon in plastic
   */
   TEfficiency* sel_eff_graphite = new TEfficiency("sel_eff_graphite","CCQE-like selection efficiency;E [Gev];#epsilon",10,0,5);
   TEfficiency* sel_eff_carbon_in_plastic = new TEfficiency("eff","CCQE-like selection efficiency;E [Gev];#epsilon",10,0,5);

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

  for (size_t i = 0; i < nof_entries; i++){
  // for (size_t i = 0; i < 10; i++){
    
    tree->GetEntry(i);
    bool is_signal = (CCQEonHydrogen==1);
    bool is_event_selected = (NofFinalStateChargedParticles==1 && pass_nof_wires_cut==1 && candidate_signal_event==1);
    bool is_in_graphite = (*InteractionVolume_short == "C_Target");
    bool is_in_plastic = (*InteractionVolume_short == "C3H6_Target");
    bool is_on_H = (*InteractionTarget == "proton");
    bool is_on_C = (*InteractionTarget == "C12");

    EventType eventType = determineEventType(is_event_selected, is_signal, is_in_graphite, is_in_plastic, is_on_C, is_on_H);
    
    total_in_fv++;
    
    switch (eventType)
    {
        case SELECTED_TRUE_POSITIVE:
            total_signal++;
            total_in_plastic_on_H++;
            h_true_signal -> Fill(IncomingNeutrino_energy);
            h_reco_signal_like_true_positive -> Fill(Neutrino_reconstructed_energy_GeV);
            h_reco_signal_like -> Fill(Neutrino_reconstructed_energy_GeV);
            h_reco_signal_like_plastic -> Fill(Neutrino_reconstructed_energy_GeV);
            break;
        
        case FALSE_NEGATIVE:
            total_signal++;
            total_in_plastic_on_H++;
            h_true_signal -> Fill(IncomingNeutrino_energy);
            break;

        case SELECTED_FALSE_POSITIVE_GRAPHITE:
            total_in_graphite++;
            total_selected_in_graphite++;
            h_true_background -> Fill(IncomingNeutrino_energy);
            h_true_C_background_graphite -> Fill(IncomingNeutrino_energy);
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
            h_true_plastic -> Fill(IncomingNeutrino_energy);
            h_reco_signal_like -> Fill(Neutrino_reconstructed_energy_GeV);
            h_reco_signal_like_plastic -> Fill(Neutrino_reconstructed_energy_GeV);
            h_reco_signal_like_on_C_plastic -> Fill(Neutrino_reconstructed_energy_GeV);
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
            h_true_graphite -> Fill(IncomingNeutrino_energy);
            sel_eff_graphite -> Fill(0, IncomingNeutrino_energy);
            break;
        
        case TRUE_NEGATIVE_PLASTIC_C:
            total_in_plastic_on_C++;
            h_true_background -> Fill(IncomingNeutrino_energy);
            h_true_plastic -> Fill(IncomingNeutrino_energy);
            h_true_C_background_plastic -> Fill(IncomingNeutrino_energy);
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

   // TRUE NEGATIVE, PLASTIC VS GRAPHITE, CARBON TARGET

  std::vector<std::string> legends = {"interactions on C of Graphite", "interactions on C of Plastic"};
  
//   plotHistograms({h_true_C_background_plastic, 
//                  h_true_C_background_graphite},
//                  "True antineutrino energy", 
//                  legends, 
//                  true,
//                  true);
    
    // CCQE-like FALSE POSITIVE, PLASTIC VS GRAPHITE, CARBON TARGET

//     legends = {"CCQE-like Carbon of Graphite", "CCQE-like Carbon of Plastic"};

//     plotHistograms({h_reco_signal_like_graphite, 
//                h_reco_signal_like_on_C_plastic},
//                "Reco antineutrino energy", 
//                legends, 
//                true,
//                true);
    
    // TRUE SIGNAL VS SELECTED CCQE-LIKE ON PLASTIC

    legends = {"True CCQE on H", "Selected CCQE-like Plastic", "Selected CCQE-like Plastic [Carbon nuclei]"};

    plotHistograms({h_true_signal, 
                    h_reco_signal_like_plastic,
                    // h_reco_signal_like_on_C_plastic
                    },
                    "Reco antineutrino energy", 
                    legends, 
                    true,
                    true,
                    true);

    // TRUE SIGNAL VS SELECTED TRUE POSITIVE

    // legends = {"True CCQE on H", "Selected CCQE-like true positive"};

    // plotHistograms({h_true_signal, 
    //                 h_reco_signal_like_true_positive,
    //                 // h_reco_signal_like_on_C_plastic
    //                 },
    //                 "Reco antineutrino energy", 
    //                 legends, 
    //                 true,
    //                 true);
  /***
  STATISTICS:
   */
//    std::cout << "number of interactions in graphite : " << h_true_graphite -> GetEntries() << "\n";
//    std::cout << "number of interactions in plastic : " << h_true_plastic -> GetEntries() << "\n";

//   /***
//   EFFICIENCIES:
//    */
//   std::vector<TEfficiency*> effs = {sel_eff_graphite, sel_eff_carbon_in_plastic};
//   std::vector<std::string> labels = {"Carbon in Plastic", "Carbon in Graphite",};
//   plotEfficiencies(effs, labels);


  return 0;
}
