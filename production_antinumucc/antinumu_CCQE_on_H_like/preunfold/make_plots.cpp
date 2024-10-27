#include <string>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

void plotHistograms(TH1D* hist1, TH1D* hist2, const std::string& canvasTitle, const std::string legendLabels[2], bool normalize = false,  bool showErrorBars = false) {

    if (hist1->GetEntries() == 0 || hist2->GetEntries() == 0) {
        std::cerr << "Errore: uno degli istogrammi non contiene voci valide." << std::endl;
        return;
    }

    // Normalizzazione opzionale degli istogrammi
    if (normalize) {
        if (hist1->Integral() > 0) hist1->Scale(1.0 / hist1->Integral());
        if (hist2->Integral() > 0) hist2->Scale(1.0 / hist2->Integral());
    }

    if (showErrorBars) {
        hist1->Sumw2();
        hist2->Sumw2();
    }

    // Calcolo del massimo tra i due istogrammi
    double max1 = hist1->GetMaximum();
    double max2 = hist2->GetMaximum();
    double maxY = std::max(max1, max2) * 1.1;  // Aggiungi un margine del 10% per lasciare spazio

    hist1->SetMaximum(maxY);  // Imposta lo stesso massimo su entrambi gli istogrammi
    hist2->SetMaximum(maxY);

    TCanvas* canvas = new TCanvas("canvas", canvasTitle.c_str(), 800, 600);
    canvas->SetTitle(canvasTitle.c_str());

    hist1->SetTitle(canvasTitle.c_str());  // Imposta il titolo per l'istogramma 1
    hist2->SetTitle(canvasTitle.c_str());  // Imposta il titolo per l'istogramma 2

    hist1->SetLineColor(kRed);
    hist2->SetLineColor(kBlue);
    hist1->SetLineWidth(2);
    hist2->SetLineWidth(2);

    // Disegna gli istogrammi con o senza barre di errore
    if (showErrorBars) {
        hist1->Draw("E HIST");  // "E" per mostrare le barre di errore
        hist2->Draw("E HIST SAME");
    } else {
        hist1->Draw("HIST");
        hist2->Draw("HIST SAME");
    }

    // Creazione della legenda
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(hist1, legendLabels[0].c_str(), "l");
    legend->AddEntry(hist2, legendLabels[1].c_str(), "l");
    legend->Draw();

    canvas->Update();
    canvas->Draw();
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
  @h_true_signal: true antineutrino spectrum for signal rate
  @h_true_background: true antineutrino spectrum for background rate
  @h_true_c_background_plastic: true antineutrino spectrum for C background rate in graphite
  @h_true_c_background_graphite: true antineutrino spectrum for C background rate in plastic

  @h_reco_singal_like: reconstructed antineutrino spectrum for selected events
  @h_reco_singal_like_plastic: reconstructed antineutrino spectrum for selected events in plastic
  @h_reco_signal_like_graphite: reconstructed antineutrino spectrum for selected events in graphite
  @h_reco_singal_like_plastic_graphite: reconstructed antineutrino spectrum for selected events in plastic with statistical subtraction of graphite background
  
  @h_reco_signal_like_on_c_plastic: reconstructed antineutrino spectrum for selected events in plastic interaction on Carbon nuclei
  @h_reco_signal_like_on_h_plastic: reconstructed antineutrino spectrum for selected events in plastic interaction on H nuclei
   */
  TH1D* h_true_signal = new TH1D ("h_true_signal", "; [GeV]", 16, 0., 8.);
  TH1D* h_true_background = new TH1D ("h_true_background", "; [GeV]", 16, 0., 8.);
  TH1D* h_true_C_background_plastic = new TH1D ("h_true_C_background_plastic", "; [GeV]", 16, 0., 8.);
  TH1D* h_true_C_background_graphite = new TH1D ("h_true_C_background_graphite", "; [GeV]", 16, 0., 8.);
  
  TH1D* h_reco_singal_like = new TH1D ("h_reco_singal_like", "; [GeV]", 16, 0., 8.);
  TH1D* h_reco_singal_like_plastic = new TH1D ("h_reco_singal_like_plastic", "; [GeV]; [GeV]", 16, 0., 8.);
  TH1D* h_reco_signal_like_graphite = new TH1D ("h_reco_signal_like_graphite", "; [GeV]", 16, 0., 8.);
  TH1D* h_reco_singal_like_plastic_graphite = new TH1D ("h_reco_singal_like_plastic_graphite", "; [GeV]", 16, 0., 8.);
  
  TH1D* h_reco_signal_like_on_C_plastic = new TH1D ("h_reco_signal_like_on_C_plastic", "; [GeV]", 16, 0., 8.);
  TH1D* h_reco_signal_like_on_H_plastic = new TH1D ("h_reco_signal_like_on_H_plastic", "; [GeV]", 16, 0., 8.);

  for (size_t i = 0; i < nof_entries; i++){
  // for (size_t i = 0; i < 10; i++){
    
    tree->GetEntry(i);
    bool is_signal = (CCQEonHydrogen==1);
    bool is_event_selected = (NofFinalStateChargedParticles==1 && pass_nof_wires_cut==1 && candidate_signal_event==1);
    bool is_in_graphite = (*InteractionVolume_short == "C_Target");
    bool is_in_plastic = (*InteractionVolume_short == "C3H6_Target");
    bool is_on_H = (*InteractionTarget == "proton");
    bool is_on_C = (*InteractionTarget == "C12");
    // std::cout << "*InteractionVolume_short " << *InteractionVolume_short 
    //           << ", is_in_graphite : "         << is_in_graphite
    //           << ", is_in_plastic : "          << is_in_plastic
    //           <<"\n";
    if(is_event_selected)
    { // selected
      h_reco_singal_like -> Fill(Neutrino_reconstructed_energy_GeV);
      if(is_signal) // selected true positive
      {
        h_true_signal ->Fill(IncomingNeutrino_energy);
      }else
      { // false positive
        h_true_background -> Fill(IncomingNeutrino_energy); // true
        if(is_in_graphite)
        { // false positive in graphite
          h_true_C_background_graphite -> Fill(IncomingNeutrino_energy); // true
          h_reco_signal_like_graphite -> Fill(Neutrino_reconstructed_energy_GeV); // reco
        }else if(is_in_plastic)
        { // false positvie in plastic
           h_reco_singal_like_plastic -> Fill(Neutrino_reconstructed_energy_GeV); // reco
           if(is_on_C)
           { // false positive in carbon of plastic
              h_true_C_background_plastic -> Fill(IncomingNeutrino_energy); // true
              h_reco_signal_like_on_C_plastic -> Fill(Neutrino_reconstructed_energy_GeV); // reco
           }else if(is_on_H)
           { // false positive in H of plastic
           }
        }
      }
    
    }else
    { // not selected
      if(is_signal)
      { // false negative
        h_true_signal ->Fill(IncomingNeutrino_energy); // true
      }else
      { // true negative
        h_true_background->Fill(IncomingNeutrino_energy); // true
        if(is_in_graphite)
        { // true negative in graphite 
          h_true_C_background_graphite -> Fill(IncomingNeutrino_energy);
        }else if(is_in_plastic)
        { //  true negative in plastic
          if(is_on_C)
          { // true negative on C of plastic
            h_true_C_background_plastic -> Fill(IncomingNeutrino_energy);
          }
        }
      }
    }
  }

  std::string legendLabels[2] = {"plastic", "graphite"};

  std::cout << "";
  std::cout << "h_true_C_background_plastic entries " << h_true_C_background_plastic->GetEntries() << "\n";
  std::cout << "h_true_C_background_graphite entries " << h_true_C_background_graphite->GetEntries() << "\n";
  plotHistograms(h_true_C_background_plastic, 
                 h_true_C_background_graphite, 
                 "True antineutrino energy interaction on C nuclei", 
                 legendLabels, 
                 true,
                 true);
  
  // std::cout << "";
  // std::cout << "h_true_C_background_plastic entries " << h_true_C_background_plastic->GetEntries() << "\n";
  // std::cout << "h_true_C_background_graphite entries " << h_reco_signal_like_graphite->GetEntries() << "\n";
  
  // plotHistograms(h_reco_signal_like_on_C_plastic, 
  //                h_reco_signal_like_graphite, 
  //                "Reconstructed antineutrino energy interaction on C nuclei", 
  //                legendLabels, 
  //                true,
  //                true);

  return 0.;
}