#include <string>
#include <vector>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

enum STAGE{
    FIDUCIAL_VOLUME,
    WIRES_CUT,
    CHARGE_MULTIPLICITY,
    ECAL_COINCIDENCE,
    STAGE_NONE
};

enum SELECTION{
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
    SELECTION_NONE
};

enum PARTICLE{
    ANTIMUON,
    NEUTRON,
    NEUTRINO,
    PARTICLE_NONE
};

enum KINETIC_VARIABLE{
    ENERGY,
    KINETIC_ENERGY,
    MOMENTUM,
    TRANSVERSE_MOMENTUM,
    LONGITUDINAL_MOMENTUM,
    ANGLE_OF_EMISSION,
    KINETIC_VARIABLE_NONE
};

enum TARGET {
    GRAPHITE,
    PLASTIC,
    OTHER,
    TARGET_NONE
};

struct histo {
    STAGE stage = STAGE_NONE;
    SELECTION selection = SELECTION_NONE;
    PARTICLE particle = PARTICLE_NONE;
    KINETIC_VARIABLE variable = KINETIC_VARIABLE_NONE;
    TARGET target = TARGET_NONE; 
    bool mc_truth = false;
    TH1D* hist = nullptr;
};

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
    tree->SetBranchAddress("Interaction_vtxZ", &Interaction_vtxZ);
    tree->SetBranchAddress("FinalStateHadronicSystemTotal4Momentum", &FinalStateHadronicSystemTotal4Momentum);
    tree->SetBranchAddress("Antimuon_p_true", &Antimuon_p_true); // MeV
    tree->SetBranchAddress("Antimuon_reconstructed_P4", &Antimuon_reconstructed_P4); // MeV

    auto nof_entries = tree->GetEntries();

    std::cout << "Number of events in fiducial volume " << nof_entries << "\n";
    
    /*** ALL HIST HAVE FIDUCIAL VOLUME CUT

        ANTIMUON SPECTRA__________________________________________________________________________________________
        @h_true_antimu: true antimuon spectrum from SIGNAL
        @h_reco_antimu: reco antimuon spectrum [STAGE: FIDUCIAL VOLUNE]
        @h_reco_antimu_stage_wires: reco antimuon spectrum [STAGE: WIRES CUT]
        @h_reco_antimu_stage_multi: reco antimuon spectrum [STAGE: CHARGE MULTIPLICITY]
        @h_reco_antimu_stage_final: reco antimuon spectrum [STAGE: ECAL COINCIDENCE]
    */
    
    uint nof_bins = 20;
    double xlow = 0.;
    double xup = 6000.;

    // NEUTRINO SPECTRA
    TH1D* h_true_antinumu = new TH1D ("h_true_antinumu", "", nof_bins, xlow, 6.);
    TH1D* h_reco_antinumu = new TH1D ("h_reco_antinumu", "", nof_bins, xlow, 6.);
    TH1D* h_true_antinumu_stage_wires = new TH1D ("h_true_antinumu_stage_wires", "", nof_bins, xlow, 6.);
    TH1D* h_reco_antinumu_stage_wires = new TH1D ("h_reco_antinumu_stage_wires", "", nof_bins, xlow, 6.);
    TH1D* h_true_antinumu_stage_multi = new TH1D ("h_true_antinumu_stage_wires", "", nof_bins, xlow, 6.);
    TH1D* h_reco_antinumu_stage_multi = new TH1D ("h_reco_antinumu_stage_wires", "", nof_bins, xlow, 6.);
    TH1D* h_true_antinumu_stage_final = new TH1D ("h_true_antinumu_stage_final", "", nof_bins, xlow, 6.);
    TH1D* h_reco_antinumu_stage_final = new TH1D ("h_reco_antinumu_stage_final", "", nof_bins, xlow, 6.);

    
    // ANTIMUON SPECTRA
    // energy
    TH1D* h_true_antimu = new TH1D ("h_true_antimu", "", nof_bins, xlow, xup);
    TH1D* h_reco_antimu = new TH1D ("h_reco_antimu", "", nof_bins, xlow, xup);
    TH1D* h_true_antimu_stage_wires = new TH1D ("h_true_antimu_stage_wires", "", nof_bins, xlow, xup);
    TH1D* h_reco_antimu_stage_wires = new TH1D ("h_reco_antimu_stage_wires", "", nof_bins, xlow, xup);
    TH1D* h_true_antimu_stage_multi = new TH1D ("h_true_antimu_stage_wires", "", nof_bins, xlow, xup);
    TH1D* h_reco_antimu_stage_multi = new TH1D ("h_reco_antimu_stage_wires", "", nof_bins, xlow, xup);
    TH1D* h_true_antimu_stage_final = new TH1D ("h_true_antimu_stage_final", "", nof_bins, xlow, xup);
    TH1D* h_reco_antimu_stage_final = new TH1D ("h_reco_antimu_stage_final", "", nof_bins, xlow, xup);
    // angle of emission
    TH1D* h_true_antimu_angle = new TH1D("h_true_antimu_angle", "", nof_bins, -3.14/2., 3.14/2.);
    
    // histos
    histo myHisto {
    FIDUCIAL_VOLUME,
    SELECTED_TRUE_POSITIVE,
    ANTIMUON,
    KINETIC_ENERGY,
    TARGET_NONE,
    true,
    new TH1D("hist_name", "Histogram Title", 100, 0, 10)
    };


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
        // guard clauses tecnique

        if(!is_signal){
            continue;
        }

        h_true_antinumu -> Fill(IncomingNeutrino_energy);
        h_true_antimu -> Fill(antimuon_true_energy);
        h_reco_antimu -> Fill(antimuon_reco_energy);

        h_true_antimu_angle -> Fill(antimuon_true_angle);

        if(!pass_nof_wires_cut){
            continue;
        }

        h_true_antinumu_stage_wires -> Fill(IncomingNeutrino_energy);
        h_true_antimu_stage_wires -> Fill(antimuon_true_energy);
        h_reco_antimu_stage_wires -> Fill(antimuon_reco_energy);

        if(NofFinalStateChargedParticles!=1){
            continue;
        }

        h_true_antinumu_stage_multi -> Fill(IncomingNeutrino_energy);
        h_true_antimu_stage_multi -> Fill(antimuon_true_energy);
        h_reco_antimu_stage_multi -> Fill(antimuon_reco_energy);

        if(!candidate_signal_event){
            continue;
        }

        h_true_antinumu_stage_final -> Fill(IncomingNeutrino_energy);
        h_true_antimu_stage_final -> Fill(antimuon_true_energy);
        h_reco_antimu_stage_final -> Fill(antimuon_reco_energy);
    }

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    h_true_antinumu->SetLineColor(kRed);
    h_true_antinumu->SetLineWidth(2);
    h_true_antinumu->Scale(1./h_true_antinumu->Integral());
    h_true_antinumu_stage_wires->SetLineColor(kCyan);
    h_true_antinumu_stage_wires->SetLineWidth(2);
    h_true_antinumu_stage_wires->Scale(1./h_true_antinumu_stage_wires->Integral());
    h_true_antinumu_stage_multi->SetLineColor(kMagenta);
    h_true_antinumu_stage_multi->SetLineWidth(2);
    h_true_antinumu_stage_multi->Scale(1./h_true_antinumu_stage_multi->Integral());
    h_true_antinumu_stage_final->SetLineColor(kBlue);
    h_true_antinumu_stage_final->SetLineWidth(2);
    h_true_antinumu_stage_final->Scale(1./h_true_antinumu_stage_final->Integral());
    
    h_true_antimu->SetLineColor(kRed);
    h_true_antimu->SetLineWidth(2);
    h_true_antimu->Scale(1./h_true_antimu->Integral());
    
    h_reco_antimu->SetLineColor(kBlack);
    h_reco_antimu->SetLineWidth(2);
    h_reco_antimu->SetMarkerStyle(20);
    h_reco_antimu->Scale(1./h_reco_antimu->Integral());
    
    h_reco_antimu_stage_wires->SetLineColor(kBlack);
    h_reco_antimu_stage_wires->SetLineWidth(2);
    h_reco_antimu_stage_wires->SetMarkerStyle(20);
    h_reco_antimu_stage_wires->Scale(1./h_reco_antimu_stage_wires->Integral());
    
    h_reco_antimu_stage_final->SetLineColor(kBlack);
    h_reco_antimu_stage_final->SetLineWidth(2);
    h_reco_antimu_stage_final->SetMarkerStyle(20);
    h_reco_antimu_stage_final->Scale(1./h_reco_antimu_stage_final->Integral());
    h_reco_antimu_stage_final->SetMaximum(0.14);
    
    std::vector<int> colors = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kYellow+1};
    
    //________________________________________
    // TCanvas* canvas = new TCanvas("canvas", "", 1500, 700);
    // canvas->Divide(3);
    
    // canvas->cd(1);

    // auto rp = new TRatioPlot(h_true_antimu, h_reco_antimu);
    // rp->Draw();
    // rp->GetUpperPad()->cd();
    // TLegend* legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    // legend1->AddEntry("h_true_antimu","True #mu+ Signal","l");
    // legend1->AddEntry("h_reco_antimu","Reco #mu+","le");
    // legend1->Draw();

    // canvas -> cd(2);
    
    // auto rp2 = new TRatioPlot(h_true_antimu, h_reco_antimu_stage_wires);
    // rp2->Draw();
    // rp2->GetUpperPad()->cd();
    // TLegend* legend2 = new TLegend(0.6, 0.6, 0.9, 0.9);
    // legend2->AddEntry("h_true_antimu","True #mu+ Signal","l");
    // legend2->AddEntry("h_reco_antimu_stage_wires","Reco #mu+ [nof wires > 70]","le");
    // legend2->Draw();

    // canvas -> cd(3);
    
    // auto rp3 = new TRatioPlot(h_true_antimu, h_reco_antimu_stage_final);
    // h_reco_antimu_stage_final->Draw("E1");
    // rp3->Draw("SAME");
    // rp3->GetUpperPad()->cd();
    // TLegend* legend3 = new TLegend(0.6, 0.6, 0.9, 0.9);
    // // h_true_antimu->Draw("HIST SAME"),
    // legend3->AddEntry("h_true_antimu","True #mu+ Signal","l");
    // legend3->AddEntry("h_reco_antimu_stage_wires","Reco #mu+ [CCQE-like]","le");
    // legend3->Draw();
    //________________________________________
    TCanvas* canvas_antimu_stages = new TCanvas("canvas_antimu_stages", "", 900, 700);

    h_true_antinumu->SetMaximum(0.15);
    h_true_antinumu->Draw("E HIST");
    h_true_antinumu_stage_wires->SetMaximum(0.15);
    h_true_antinumu_stage_wires->Draw("E HIST SAME");
    h_true_antinumu_stage_multi->SetMaximum(0.15);
    h_true_antinumu_stage_multi->Draw("E HIST SAME");
    h_true_antinumu_stage_final->SetMaximum(0.15);
    h_true_antinumu_stage_final->Draw("E HIST SAME");

    TLegend* l1 = new TLegend(0.6, 0.6, 0.9, 0.9);
    l1 -> AddEntry(h_true_antinumu, "STAGE 1: FIDUCIAL");
    l1 -> AddEntry(h_true_antinumu_stage_wires, "STAGE 2: WIRES");
    l1 -> AddEntry(h_true_antinumu_stage_multi, "STAGE 3: MULTIPL.");
    l1 -> AddEntry(h_true_antinumu_stage_final, "STAGE 4: ECAL");

    l1->Draw();

    //________________________________________
    // TCanvas* canvas_antimu_stages = new TCanvas("canvas_antinumu_stages", "", 900, 700);

    // h_true_antimu->Scale(1.0 / h_true_antimu->Integral());
    // h_true_antimu->SetMaximum(0.15);
    // h_true_antimu->Draw("E HIST");

    // h_true_antimu_stage_wires->SetLineColor(kCyan);
    // h_true_antimu_stage_wires->Scale(1.0 / h_true_antimu_stage_wires->Integral());
    // h_true_antimu_stage_wires->SetMaximum(0.15);
    // h_true_antimu_stage_wires->Draw("E HIST SAME");

    // h_true_antimu_stage_multi->SetLineColor(kMagenta);
    // h_true_antimu_stage_multi->Scale(1.0 / h_true_antimu_stage_multi->Integral());
    // h_true_antimu_stage_multi->SetMaximum(0.15);
    // h_true_antimu_stage_multi->Draw("E HIST SAME");

    // h_true_antimu_stage_final->SetLineColor(kBlue);
    // h_true_antimu_stage_final->Scale(1.0 / h_true_antimu_stage_final->Integral());
    // h_true_antimu_stage_final->SetMaximum(0.15);
    // h_true_antimu_stage_final->Draw("E HIST SAME");


    // TLegend* l1 = new TLegend(0.6, 0.6, 0.9, 0.9);
    // l1 -> AddEntry(h_true_antimu, "STAGE 1: FIDUCIAL");
    // l1 -> AddEntry(h_true_antimu_stage_wires, "STAGE 2: WIRES");
    // l1 -> AddEntry(h_true_antimu_stage_multi, "STAGE 3: MULTIPL.");
    // l1 -> AddEntry(h_true_antimu_stage_final, "STAGE 4: ECAL");

    // l1->Draw();

    //________________________________________


    return 0;
}