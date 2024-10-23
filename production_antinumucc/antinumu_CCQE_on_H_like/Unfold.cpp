// gInterpreter->AddIncludePath("/opt/exp_software/neutrino/ROOUNFOLD/RooUnfold/build/");

#include <string>
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

int Unfold(){
    TFile file("/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/antinumu_CCQE_on_H_like/reports/events-in-SANDtracker.0.to.100.preunfold.root","READ");

    TTree* tree = (TTree*)file.Get("preunfold");
    int CCQEonHydrogen;
    int NofFinalStateChargedParticles;
    bool pass_nof_wires_cut;
    int candidate_signal_event; // 
    std::string* InteractionTarget = new std::string();
    std::string* InteractionVolume_short = new std::string();
    double IncomingNeutrino_energy;
    double Neutrino_reconstructed_energy_GeV;

    RooUnfoldResponse response (16, 0., 8.);

    tree->SetBranchAddress("CCQEonHydrogen", &CCQEonHydrogen);
    tree->SetBranchAddress("NofFinalStateChargedParticles", &NofFinalStateChargedParticles);
    tree->SetBranchAddress("pass_nof_wires_cut", &pass_nof_wires_cut);
    tree->SetBranchAddress("candidate_signal_event", &candidate_signal_event);
    tree->SetBranchAddress("InteractionTarget", &InteractionTarget);
    tree->SetBranchAddress("InteractionVolume_short", &InteractionVolume_short);
    tree->SetBranchAddress("IncomingNeutrino_energy", &IncomingNeutrino_energy);
    tree->SetBranchAddress("Neutrino_reconstructed_energy_GeV", &Neutrino_reconstructed_energy_GeV);

    auto nof_entries = tree->GetEntries();
    int total=0, total_true=0, nof_true_selected=0, nof_true_nof_excluded=0;
    for (size_t i = 0; i < nof_entries; i++)
    {
        tree->GetEntry(i); 
        if(CCQEonHydrogen==1){
                total_true++;
            if (NofFinalStateChargedParticles==1 && pass_nof_wires_cut==1 && candidate_signal_event==1)
            {
                response.Fill(Neutrino_reconstructed_energy_GeV, IncomingNeutrino_energy);
                nof_true_selected++;
            }else{
                response.Miss(IncomingNeutrino_energy);
                nof_true_nof_excluded++;

            }
        }
        total++;
    }
    std::cout << "total " << total << "\n";
    std::cout << "total_true " << total_true << "\n";
    std::cout << "nof_true_selected " << nof_true_selected << "\n";
    std::cout << "nof_true_nof_excluded " << nof_true_nof_excluded << "\n";
    return 0.;
}