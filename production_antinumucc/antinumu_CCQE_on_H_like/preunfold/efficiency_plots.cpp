#include <string>
#include <vector>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

enum STAGE{
    FIDUCIAL_VOLUME = 0,
    WIRES_CUT = 1,
    CHARGE_MULTIPLICITY = 2,
    ECAL_COINCIDENCE = 3,
    STAGE_NONE = 4
};

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
    SELECTION_NONE
};

enum PARTICLE{
    ANTIMUON = 0,
    NEUTRON = 1,
    NEUTRINO = 2,
    PARTICLE_NONE = 3
};

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

std::vector<int> colors = {kRed, kBlue, kGreen, kBlack, kMagenta, kCyan, kOrange, kYellow+1};

TH1D* calculate_bin_ratio(TH1D* hist1, TH1D* hist2) {
    // Controlla che entrambi gli istogrammi abbiano lo stesso numero di bin
    if (hist1->GetNbinsX() != hist2->GetNbinsX()) {
        std::cerr << "Error: Histograms have different binning!" << std::endl;
        return nullptr;
    }

    // Clona il primo istogramma per creare l'istogramma del rapporto
    TH1D* ratio_hist = (TH1D*)hist1->Clone("ratio_hist");
    ratio_hist->SetTitle("Bin-by-bin ratio (hist2 / hist1)"); // Imposta il titolo
    ratio_hist->SetLineColor(kBlack); // Imposta il colore della linea del rapporto

    // Cicla su ogni bin e calcola il rapporto
    for (int bin = 1; bin <= hist1->GetNbinsX(); bin++) {
        double val1 = hist1->GetBinContent(bin); // Valore del primo istogramma
        double val2 = hist2->GetBinContent(bin); // Valore del secondo istogramma
        double err1 = hist1->GetBinError(bin); // Errore sul primo istogramma
        double err2 = hist2->GetBinError(bin); // Errore sul secondo istogramma

        if (val1 != 0) { // Calcola il rapporto solo se il valore nel primo istogramma non è zero
            double ratio = val2 / val1;
            // Calcola l'errore del rapporto usando la propagazione degli errori
            double ratio_error = ratio * sqrt((err1 / val1) * (err1 / val1) + (err2 / val2) * (err2 / val2));
            ratio_hist->SetBinContent(bin, ratio);
            ratio_hist->SetBinError(bin, ratio_error);
        } else {
            ratio_hist->SetBinContent(bin, 0); // Se val1 è zero, imposta il rapporto a 0 o un altro valore speciale
            ratio_hist->SetBinError(bin, 0); // Errore a 0 in questo caso
        }
    }

    return ratio_hist;
}


void plot_histograms(TH1D* histos[], int num_stages, double max_y_norm = 0.2) {

    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);  // Posizione della legenda (x1, y1, x2, y2)
    legend->AddEntry(histos[FIDUCIAL_VOLUME], "FIDUCIAL_VOLUME", "l");
    legend->AddEntry(histos[WIRES_CUT], "WIRES_CUT", "l");
    legend->AddEntry(histos[CHARGE_MULTIPLICITY], "CHARGE_MULTIPLICITY", "l");
    legend->AddEntry(histos[ECAL_COINCIDENCE], "ECAL_COINCIDENCE", "l");

    TCanvas* canvas = new TCanvas("canvas", "", 1300, 700);
    canvas->Divide(2);
    canvas->cd(1);
    
    TH1D* ratio_hist[3];

    for (int i = 0; i < num_stages; i++) {
        histos[i]->SetLineColor(colors[i % colors.size()]); 
        histos[i]->SetLineWidth(2);
        if (i == 0) {
            histos[i]->Draw("E HIST");
            histos[i]->SetTitle("test");
        } else {
            histos[i]->Draw("E HIST SAME");
            ratio_hist[i-1] = calculate_bin_ratio(histos[i-1], histos[i]);
            ratio_hist[i-1] -> SetLineColor(colors[i % colors.size()]);
            ratio_hist[i-1] -> SetMarkerStyle(kFullDiamond);
        }
    }

    legend->Draw();

    canvas->cd(2);
    TH1D* clone = nullptr;
    for (int i = 0; i < num_stages; i++) {
        clone = reinterpret_cast<TH1D*>(histos[i]->Clone(TString::Format("clone_%d",i).Data()));
        clone->Scale(1.0 / histos[i]->Integral());
        clone->SetMaximum(max_y_norm);
        clone->SetLineColor(colors[i % colors.size()]);  // Ciclo sui colori se ci sono più stadi
        clone->SetLineWidth(2);
        if (i == 0) {
            clone->Draw("E HIST");
        } else {
            clone->Draw("E HIST SAME");
        }
    }

    legend->Draw();

    // Aggiorna il canvas
    canvas->Update();

    TCanvas* canvas_eff = new TCanvas("canvas_eff", "", 900, 700);
    for (size_t i = 0; i < 3; i++)
    {
        if(i==0){
            ratio_hist[i] -> Draw("E");
        }else{
            ratio_hist[i] -> Draw("E SAME");
        }
    }
    
}

void plot_histograms2(TH2D* h2, std::string canvas_name){
    TCanvas* c_h2 = new TCanvas(canvas_name.c_str(), "", 900, 700);
    h2->Draw("COLZ");
    c_h2->Draw();
}

void plot_reco(TH1D* h_true, TH1D* h_reco, std::string canvas_name = "", std::string pdf_name = ""){
    
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(h_true, "true", "l");
    legend->AddEntry(h_reco, "reco", "l");

    TCanvas* c_tr = new TCanvas(canvas_name.c_str(), "", 900, 700);
    

    h_true -> SetLineWidth(2);
    h_true -> SetLineColor(kRed);

    h_reco -> SetLineColor(kBlack);
    h_reco -> SetLineWidth(2);
    h_reco -> SetMarkerStyle(20);
    
    h_true -> Draw("E HIST");
    h_reco -> Draw("E SAME");

    if (!canvas_name.empty()) {
        TPaveText* pave = new TPaveText(0.1, 0.93, 0.9, 0.98, "NDC");
        pave->AddText(canvas_name.c_str());
        pave->SetTextAlign(22);  // Centro all'interno della TPaveText
        pave->SetFillColor(0);   // Colore di riempimento trasparente
        pave->SetBorderSize(0);  // Nessun bordo
        pave->Draw();
    }

    c_tr -> Draw();
    legend -> Draw();

    if (!pdf_name.empty()) {
        c_tr->SaveAs(pdf_name.c_str());
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

/***
COMPARE: STAGES, array of TH1D to compare true spectra at any stage -> estimate efficiencies
 */

TH1D* antineutrino_true_energy[STAGE_NONE] = {
    new TH1D("antineutrino_true_energy_FIDUCIAL_VOLUME", "", 12u, 0., 6.),
    new TH1D("antineutrino_true_energy_WIRES_CUT", "", 12u, 0., 6.),
    new TH1D("antineutrino_true_energy_CHARGE_MULTIPLICITY", "", 12u, 0., 6.),
    new TH1D("antineutrino_true_energy_ECAL_COINCIDENCE", "", 12u, 0., 6.)
};

TH1D* antimu_true_energy[STAGE_NONE] = {
    new TH1D("antimuon_true_energy_FIDUCIAL_VOLUME", "", 20u, 0., 6000.),
    new TH1D("antimuon_true_energy_WIRES_CUT", "", 20u, 0., 6000.),
    new TH1D("antimuon_true_energy_CHARGE_MULTIPLICITY", "", 20u, 0., 6000.),
    new TH1D("antimuon_true_energy_ECAL_COINCIDENCE", "", 20u, 0., 6000.)
};

TH1D* antimu_reco_energy[STAGE_NONE] = {
    new TH1D("antimuon_reco_energy_FIDUCIAL_VOLUME", "", 20u, 0., 6000.),
    new TH1D("antimuon_reco_energy_WIRES_CUT", "", 20u, 0., 6000.),
    new TH1D("antimuon_reco_energy_CHARGE_MULTIPLICITY", "", 20u, 0., 6000.),
    new TH1D("antimuon_reco_energy_ECAL_COINCIDENCE", "", 20u, 0., 6000.)
};

TH1D* neutron_true_kin_energy[STAGE_NONE] = {
    new TH1D("neutron_true_kin_energy_FIDUCIAL_VOLUME", "", 20u, 0., 1.),
    new TH1D("neutron_true_kin_energy_WIRES_CUT", "", 20u, 0., 1.),
    new TH1D("neutron_true_kin_energy_CHARGE_MULTIPLICITY", "", 20u, 0., 1.),
    new TH1D("neutron_true_kin_energy_ECAL_COINCIDENCE", "", 20u, 0., 1.)
};

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

/***
RECO: Reconstructed quantities after selection
 */
TH1D* antineutrino_reco_energy_CCQEonHLIKE[TARGET_NONE] = {
    new TH1D("antineutrino_reco_energy_SELECTED_GRAPHITE", "", 12u, 0., 6.),
    new TH1D("antineutrino_reco_energy_SELECTED_PLASTIC", "", 12u, 0., 6.),
    new TH1D("antineutrino_reco_energy_SELECTED_OTHER", "", 12u, 0., 6.),
    new TH1D("antineutrino_reco_energy_SELECTED_ANY_TARGET", "", 12u, 0., 6.),
};


void EnableSumw2(TH1D* histograms[], size_t size) {
    for (size_t i = 0; i < size; ++i) {
        if (histograms[i]) { // Controlla che l'istogramma esista
            histograms[i]->Sumw2(); // Abilita Sumw2
        }
    }
}

void EnableSumw2_all_hists(){
    EnableSumw2(antineutrino_true_energy, STAGE_NONE);
    EnableSumw2(antimu_true_energy, STAGE_NONE);
    EnableSumw2(antimu_reco_energy, STAGE_NONE);
    EnableSumw2(neutron_true_kin_energy, STAGE_NONE);
    EnableSumw2(antineutrino_reco_energy_CCQEonHLIKE, TARGET_NONE);
}


bool IsInsideSpot(const double vtxX, const double vtxY, const double vtxZ){
    TVector3 spot_center = {0., -2.5, 23.5};
    TVector3 vtx = {vtxX, vtxY, vtxZ};
    return ((spot_center - vtx).Mag() < 1.);
}

void Fill_Final_Histos(SELECTION selection, double reco_nu_E){
    switch (selection)
    {
        case SELECTED_TRUE_POSITIVE:
            antineutrino_reco_energy_CCQEonHLIKE[PLASTIC] -> Fill(reco_nu_E);
            break;

        case FALSE_NEGATIVE:
            break;

        // FALSE POSITIVES
        case SELECTED_FALSE_POSITIVE_GRAPHITE:
            break;

        case SELECTED_FALSE_POSITIVE_PLASTIC_C:
            antineutrino_reco_energy_CCQEonHLIKE[PLASTIC] -> Fill(reco_nu_E);
            break;

        case SELECTED_FALSE_POSITIVE_PLASTIC_H:
            antineutrino_reco_energy_CCQEonHLIKE[PLASTIC] -> Fill(reco_nu_E);
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

int efficiency_plots(){

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
    TVector3* Antimuon_p_true = nullptr;
    TLorentzVector* Antimuon_reconstructed_P4 = nullptr;
    TLorentzVector* FinalStateHadronicSystemTotal4Momentum = nullptr; 
    TVector3* PredictedNeutron_P3_GeV = nullptr; 
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
    tree->SetBranchAddress("FinalStateHadronicSystemTotal4Momentum", &FinalStateHadronicSystemTotal4Momentum);
    tree->SetBranchAddress("Antimuon_p_true", &Antimuon_p_true); // MeV
    tree->SetBranchAddress("Antimuon_reconstructed_P4", &Antimuon_reconstructed_P4); // MeV

    auto nof_entries = tree->GetEntries();

    std::cout << "Number of events in fiducial volume " << nof_entries << "\n";
    
    EnableSumw2_all_hists();
    
    for (size_t i = 0; i < nof_entries; i++)
    // for (size_t i = 0; i < 1000; i++)
    {
        tree->GetEntry(i);

        bool is_signal = (CCQEonHydrogen==1);
        bool is_event_selected = (NofFinalStateChargedParticles==1 && pass_nof_wires_cut==1 && candidate_signal_event==1);
        bool is_in_graphite = (*InteractionVolume_short == "C_Target");
        bool is_in_plastic = (*InteractionVolume_short == "C3H6_Target");
        bool is_on_H = (*InteractionTarget == "proton");
        bool is_on_C = (*InteractionTarget == "C12");

        double antimuon_true_energy = sqrt(Antimuon_p_true->Mag() * Antimuon_p_true->Mag() +  105.658 * 105.658);
        double antimuon_reco_energy = Antimuon_reconstructed_P4->T();
        double antimuon_true_angle = Antimuon_p_true->Angle({0.,0.,1.});
        double hadron_syst_kin_energy = FinalStateHadronicSystemTotal4Momentum->T() -  0.939565;

        SELECTION selection = DetermineSelection(is_event_selected, is_signal, is_in_graphite, is_in_plastic, is_on_C, is_on_H);

        Fill_Final_Histos(selection, Neutrino_reconstructed_energy_GeV);

        if(!is_signal){
            continue;
        }
        
        // if(!IsInsideSpot(Interaction_vtxX, Interaction_vtxY, Interaction_vtxZ)){
        //     continue;
        // }

        antineutrino_true_energy[FIDUCIAL_VOLUME] -> Fill(IncomingNeutrino_energy);
        antimu_true_energy[FIDUCIAL_VOLUME] -> Fill(antimuon_true_energy);
        neutron_true_kin_energy[FIDUCIAL_VOLUME] -> Fill(hadron_syst_kin_energy);
        
        antimu_reco_energy[FIDUCIAL_VOLUME] -> Fill(antimuon_reco_energy);

        true_nu_E_vs_true_n_kin_E[FIDUCIAL_VOLUME] -> Fill(IncomingNeutrino_energy, hadron_syst_kin_energy);
        true_nu_E_vs_true_mu_E[FIDUCIAL_VOLUME] -> Fill(IncomingNeutrino_energy, antimuon_true_energy);
        true_mu_E_vs_true_n_kin_E[FIDUCIAL_VOLUME] -> Fill(antimuon_true_energy, hadron_syst_kin_energy);

        if(!pass_nof_wires_cut){
            continue;
        }

        antineutrino_true_energy[WIRES_CUT] -> Fill(IncomingNeutrino_energy);
        antimu_true_energy[WIRES_CUT] -> Fill(antimuon_true_energy);
        neutron_true_kin_energy[WIRES_CUT] -> Fill(hadron_syst_kin_energy);
        
        antimu_reco_energy[WIRES_CUT] -> Fill(antimuon_reco_energy);

        true_nu_E_vs_true_n_kin_E[WIRES_CUT] -> Fill(IncomingNeutrino_energy, hadron_syst_kin_energy);
        true_nu_E_vs_true_mu_E[WIRES_CUT] -> Fill(IncomingNeutrino_energy, antimuon_true_energy);
        true_mu_E_vs_true_n_kin_E[WIRES_CUT] -> Fill(antimuon_true_energy, hadron_syst_kin_energy);

        if(NofFinalStateChargedParticles!=1){
            continue;
        }
        antineutrino_true_energy[CHARGE_MULTIPLICITY] -> Fill(IncomingNeutrino_energy);
        antimu_true_energy[CHARGE_MULTIPLICITY] -> Fill(antimuon_true_energy);
        neutron_true_kin_energy[CHARGE_MULTIPLICITY] -> Fill(hadron_syst_kin_energy);
        
        antimu_reco_energy[CHARGE_MULTIPLICITY] -> Fill(antimuon_reco_energy);

        true_nu_E_vs_true_n_kin_E[CHARGE_MULTIPLICITY] -> Fill(IncomingNeutrino_energy, hadron_syst_kin_energy);
        true_nu_E_vs_true_mu_E[CHARGE_MULTIPLICITY] -> Fill(IncomingNeutrino_energy, antimuon_true_energy);
        true_mu_E_vs_true_n_kin_E[CHARGE_MULTIPLICITY] -> Fill(antimuon_true_energy, hadron_syst_kin_energy);

        if(!candidate_signal_event){
            continue;
        }

        antineutrino_true_energy[ECAL_COINCIDENCE] -> Fill(IncomingNeutrino_energy);
        antimu_true_energy[ECAL_COINCIDENCE] -> Fill(antimuon_true_energy);
        neutron_true_kin_energy[ECAL_COINCIDENCE] -> Fill(hadron_syst_kin_energy);
        
        antimu_reco_energy[ECAL_COINCIDENCE] -> Fill(antimuon_reco_energy);

        true_nu_E_vs_true_n_kin_E[ECAL_COINCIDENCE] -> Fill(IncomingNeutrino_energy, hadron_syst_kin_energy);
        true_nu_E_vs_true_mu_E[ECAL_COINCIDENCE] -> Fill(IncomingNeutrino_energy, antimuon_true_energy);
        true_mu_E_vs_true_n_kin_E[ECAL_COINCIDENCE] -> Fill(antimuon_true_energy, hadron_syst_kin_energy);
    }

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    /**
    @true: compare true quantitie for different
     */

    // plot_histograms(antineutrino_true_energy, STAGE_NONE, 0.16);
    // plot_histograms(antimu_true_energy, STAGE_NONE, 0.16);
    // plot_histograms(neutron_true_kin_energy, STAGE_NONE, 0.35);
    
    /**
    @reco: compare reconstructed quantities for each stage
     */
    // plot_reco(antimu_true_energy[WIRES_CUT], antimu_reco_energy[WIRES_CUT], 
    //         "antimuon energy CCQE on H [WIRES CUT]", "antimuon_energy_CCQEonH_wires_cut.png");
    // plot_reco(antimu_true_energy[ECAL_COINCIDENCE], antimu_reco_energy[ECAL_COINCIDENCE], 
    //         "antimuon energy CCQE on H [ECAL COINCIDENCE]", "antimuon_energy_CCQEonH_ecal_coincidence.png");

    /**
    @relations: plot true relations
     */
    // plot_histograms2(true_nu_E_vs_true_n_kin_E[FIDUCIAL_VOLUME], "c1"); // nu vs n
    // plot_histograms2(true_nu_E_vs_true_mu_E[FIDUCIAL_VOLUME], "c2"); // nu vs mu
    // plot_histograms2(true_mu_E_vs_true_n_kin_E[FIDUCIAL_VOLUME], "c3"); // mu vs n
    
    TCanvas *c1 = new TCanvas("c1", "TRatioPlot Normalized", 900, 700);

    // Clona gli istogrammi originali per non modificarli
    TH1D* hTrue = (TH1D*)antineutrino_true_energy[FIDUCIAL_VOLUME]->Clone("hTrue");
    TH1D* hReco = (TH1D*)antineutrino_reco_energy_CCQEonHLIKE[PLASTIC]->Clone("hReco");

    // Normalizza gli istogrammi
    if (hTrue->Integral() != 0) hTrue->Scale(1. / hTrue->Integral());
    if (hReco->Integral() != 0) hReco->Scale(1. / hReco->Integral());

    // Configura gli istogrammi normalizzati
    hTrue->SetTitle(";E_{#nu} [GeV]; counts / tot counts");
    hTrue->SetLineColor(kRed);
    hTrue->SetLineStyle(1);
    hTrue->SetLineWidth(2);
    hTrue->SetMaximum(0.25);

    hReco->SetLineColor(kBlack);
    hReco->SetMarkerStyle(20);
    hReco->SetMarkerSize(0.8);
    hReco->SetMaximum(0.25);

    // Crea il TRatioPlot con gli istogrammi normalizzati
    TRatioPlot *rp = new TRatioPlot(hTrue, hReco);

    // Configura le opzioni di disegno per il TRatioPlot
    rp->SetH1DrawOpt("E HIST"); // Disegna il primo istogramma con linea continua
    rp->SetH2DrawOpt("E1");   // Disegna il secondo istogramma con barre di errore

    // Disegna il TRatioPlot
    rp->Draw();
    rp->GetLowerRefYaxis()->SetTitle("True / Reco");
    rp->GetLowerRefGraph()->SetMinimum(0.0);
    rp->GetLowerRefGraph()->SetMaximum(1.5);

    rp->GetUpperPad()->cd();
    TLegend *legend = new TLegend(0.6, 0.6, 0.9, 0.9);
    legend->AddEntry("hTrue", "E_{ #bar{#nu}_{#mu}}^{True} True CCQE on H", "l");
    legend->AddEntry("hReco", "E_{ #bar{#nu}_{#mu}}^{Reco} Selected CCQE-like Plastic", "le");
    legend->Draw();

    // Aggiorna il canvas
    c1->Update();

    return 0;
}
