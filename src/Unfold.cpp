// gInterpreter->AddIncludePath("/opt/exp_software/neutrino/ROOUNFOLD/RooUnfold/build/");

#include <string>
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

int Unfold(){
    
    uint run_start = 101; // production_0001, run_1000 -> run_1999
    uint production = run_start / 100;
    auto FOLDER_PRODUCTION = TString::Format("/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc/production_%04d/", production);
    auto folder_reports = TString::Format("%spreunfold/", FOLDER_PRODUCTION);
    
    // auto index = 0u;
    // unsigned int files_per_jobs = 100u;
    
    int CCQEonHydrogen;
    int NofFinalStateChargedParticles;
    bool pass_nof_wires_cut;
    int candidate_signal_event; // 
    std::string* InteractionTarget = new std::string();
    std::string* InteractionVolume_short = new std::string();
    double IncomingNeutrino_energy;
    double Neutrino_reconstructed_energy_GeV;

    RooUnfoldResponse response(16, 0., 8.);
    TH1D* h_true = new TH1D ("true", "Neutrino energy; [GeV]", 16, 0., 8.);
    TH1D* h_reco = new TH1D ("reco", "Neutrino energy; [GeV]",  16, 0., 8.);
    
    int total=0, total_true=0, nof_true_selected=0, nof_true_nof_excluded=0;

    for (size_t index = 0; index < 10; index++) // 10 is the nof files to be processed (1 file <-> 100 file edep <-> 1mln events)
    {
        unsigned int file_start = index * files_per_jobs;
        unsigned int file_stop = index * files_per_jobs + files_per_jobs;
        
        auto file_report = TString::Format("%sevents-in-SANDtracker.%d.to.%d.preunfold.root", folder_reports.Data(), file_start, file_stop);
        TFile TFile_report(file_report.Data());
        TTree* tree = (TTree*)TFile_report.Get("preunfold");
        auto nof_entries = tree->GetEntries();
        
        tree->SetBranchAddress("CCQEonHydrogen", &CCQEonHydrogen);
        tree->SetBranchAddress("NofFinalStateChargedParticles", &NofFinalStateChargedParticles);
        tree->SetBranchAddress("pass_nof_wires_cut", &pass_nof_wires_cut);
        tree->SetBranchAddress("candidate_signal_event", &candidate_signal_event);
        tree->SetBranchAddress("InteractionTarget", &InteractionTarget);
        tree->SetBranchAddress("InteractionVolume_short", &InteractionVolume_short);
        tree->SetBranchAddress("IncomingNeutrino_energy", &IncomingNeutrino_energy);
        tree->SetBranchAddress("Neutrino_reconstructed_energy_GeV", &Neutrino_reconstructed_energy_GeV);

        for (size_t i = 0; i < nof_entries; i++)
        {
            tree->GetEntry(i); 
            if(CCQEonHydrogen==1){
                    total_true++;
                    h_true->Fill(IncomingNeutrino_energy);
                if (NofFinalStateChargedParticles==1 && pass_nof_wires_cut==1 && candidate_signal_event==1)
                {
                    // CCQE ON H [PLASTIC]
                    response.Fill(Neutrino_reconstructed_energy_GeV, IncomingNeutrino_energy);
                    h_reco->Fill(Neutrino_reconstructed_energy_GeV);
                    nof_true_selected++;
                }else{
                    response.Miss(IncomingNeutrino_energy);
                    nof_true_nof_excluded++;
                }
            }
            total++;
        }
        
        std::cout << "Result analysis from production genie+edep from " << file_start << " to " << file_stop << "\n";
        std::cout << "total " << total << "\n";
        std::cout << "total_true " << total_true << "\n";
        std::cout << "nof_true_selected " << nof_true_selected << "\n";
        std::cout << "nof_true_nof_excluded " << nof_true_nof_excluded << "\n";
        
        tree->Delete();
    }
    
    RooUnfoldBayes   unfold(&response, h_reco, 4);

    TH1D* hUnfold= (TH1D*)unfold.Hunfold();

    TCanvas* c1= new TCanvas("canvas","canvas");

    unfold.PrintTable (cout, h_true);
    hUnfold->Draw();
    h_reco->Draw("SAME");
    h_true->SetLineColor(8);
    h_true->Draw("SAME");

    c1->SaveAs("RooUnfoldExample.pdf");

    return 0.;
}
// root [0] gInterpreter->AddIncludePath("/opt/exp_software/neutrino/ROOUNFOLD/RooUnfold/build/");
// root [1] .L /opt/exp_software/neutrino/ROOUNFOLD/RooUnfold/build/libRooUnfold.so
