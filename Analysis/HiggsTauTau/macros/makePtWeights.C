#include "TFile.h"
#include "TH2D.h"
#include "TTree.h"
#include "../Core/interface/Plotting.h"

void makePtWeights() {

  TFile f_cms("GluGluToHToTauTau_M-125_mt_2012_June11.root");
  TTree *t_cms = (TTree*)f_cms.Get("pt_tree");

  float y[4] = {-0.5, 0.5, 1.5, 2.5};
  std::vector<float> x;
  float s = 0.;
  while (s < 100.) {
    std::cout << s << "\n";
    x.push_back(s);
    s += 5.;
  }
  while (s <= 301.) {
    std::cout << s << "\n";
    x.push_back(s);
    s += 10.;
  }


  TH2D h_cms("h_cms", "h_cms", x.size() - 1, &(x[0]), 3, y);


  {
    double cms_pt;
    unsigned cms_n_jets;
    t_cms->SetBranchAddress("pt", &cms_pt);
    t_cms->SetBranchAddress("n_jets", &cms_n_jets);

    for (int i = 0; i < t_cms->GetEntries(); ++i) {
      t_cms->GetEntry(i);
      h_cms.Fill(TMath::Min(cms_pt, 300.*(1.-1E-5)), TMath::Min(cms_n_jets, unsigned(2)));
    }
  }
  h_cms.Scale(1. / h_cms.Integral());


  TFile f_atlas("HiggsPtTree.root");
  TTree *t_atlas = (TTree*)f_atlas.Get("HiggsPtTree");

  TH2D h_atlas("h_atlas", "h_atlas", x.size() - 1, &(x[0]), 3, y);
  {
    float atlas_pt;
    int atlas_n_jets;
    float atlas_wt;
    t_atlas->SetBranchAddress("H_pT", &atlas_pt);
    t_atlas->SetBranchAddress("Njets25", &atlas_n_jets);
    t_atlas->SetBranchAddress("weight_pT_dyn", &atlas_wt);

    for (int i = 0; i < t_atlas->GetEntries(); ++i) {
      t_atlas->GetEntry(i);
      h_atlas.Fill(TMath::Min(atlas_pt, float(300.*(1.-1E-5))), TMath::Min(atlas_n_jets, int(2)), atlas_wt);
    }
  }
  h_atlas.Scale(1. / h_atlas.Integral());

  TH2D h_weight = *(TH2D*)h_atlas.Clone("h_weight");
  h_weight.Divide(&h_cms);

  TFile out("weights_June11.root", "RECREATE");
  h_cms.Write();
  h_atlas.Write();
  h_weight.Write();
  out.Close();


}

