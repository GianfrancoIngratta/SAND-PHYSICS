#include <string>
#include <vector>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

// source /opt/exp_software/neutrino/ROOUNFOLD/RooUnfold/setup.sh

#include <TEfficiency.h>
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

void Unfold(RooUnfoldResponse* response, TH1D* hTrue, TH1D* hMeas, std::string pdf_name){
  RooUnfoldBayes   unfold (response, hMeas, 4);
  TH1D* hUnfold= (TH1D*) unfold.Hunfold();

  TCanvas* c1= new TCanvas("canvas","canvas");

  unfold.PrintTable (cout, hTrue);
  hUnfold->Draw();
  hMeas->Draw("SAME");
  hTrue->SetLineColor(8);
  hTrue->Draw("SAME");

//   c1->SaveAs("RooUnfoldExample.pdf");
  c1->SaveAs(pdf_name.c_str());
}


int unfold_decri(){
    TFile file_unfold("unfold_file.root","READ");
    auto response_sinal_sample1 = (RooUnfoldResponse*)file_unfold.Get("SIGNAL_SAMPLE1");
    auto htrue_sinal_sample1 = (TH1D*)file_unfold.Get("h_SIGNAL_SAMPLE1_true");
    auto hmeas_sinal_sample1 = (TH1D*)file_unfold.Get("h_SIGNAL_SAMPLE1_reco");
    Unfold(response_sinal_sample1, htrue_sinal_sample1, hmeas_sinal_sample1, "RooUnfoldExample.pdf");
    return 0;
}