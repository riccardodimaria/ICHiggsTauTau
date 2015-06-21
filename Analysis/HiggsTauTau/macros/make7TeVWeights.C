#include "TFile.h"
#include "TH2D.h"
#include "TTree.h"
#include "../Core/interface/Plotting.h"

void make7TeVWeights() {
  TH1::SetDefaultSumw2(true);
  TFile f_cms("output/lhchcg/GluGluToHToTauTau_M-125_mt_2011.root");
  TTree *t_cms = (TTree *)f_cms.Get("pt_tree");

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
      h_cms.Fill(TMath::Min(cms_pt, 300. * (1. - 1E-5)),
                 TMath::Min(cms_n_jets, unsigned(2)));
    }
  }
  h_cms.Scale(1. / h_cms.Integral());

  // TFile f_atlas("output/lhchcg/mc11_HSG4_ggH.root");
  // TTree *t_atlas = (TTree *)f_atlas.Get("tau");
  TFile f_atlas("ATLAS_ggF_HpT_trees.root");
  TTree *t_atlas = (TTree *)f_atlas.Get("ggF_PowPy6_7TeV_nohfact");

  TH2D h_atlas("h_atlas", "h_atlas", x.size() - 1, &(x[0]), 3, y);
  TH1F h_at_inc("atlas_inclusive", "atlas_inclusive", 60, 0, 300);
  {
    float atlas_pt;
    int atlas_n_jets;
    float atlas_wt;

    t_atlas->SetBranchAddress("H_pT", &atlas_pt);
    t_atlas->SetBranchAddress("Njets25", &atlas_n_jets);
    t_atlas->SetBranchAddress("weight", &atlas_wt);

  
    for (int i = 0; i < t_atlas->GetEntries(); ++i) {
      t_atlas->GetEntry(i);
      h_atlas.Fill(TMath::Min(atlas_pt, float(300. * (1. - 1E-5))),
                   TMath::Min(atlas_n_jets, int(2)), atlas_wt);
      h_at_inc.Fill(atlas_pt, atlas_wt);
    }
  }
  h_atlas.Scale(1. / h_atlas.Integral());

  TH2D h_weight = *(TH2D *)h_atlas.Clone("h_weight");
  h_weight.Divide(&h_cms);

  TFile out("weights_June18_7TeV.root", "RECREATE");
  h_cms.Write();
  h_atlas.Write();
  h_weight.Write();
  h_at_inc.Write();
  out.Close();
}
