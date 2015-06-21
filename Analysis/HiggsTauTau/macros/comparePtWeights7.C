#include "TFile.h"
#include "TH2D.h"
#include "TTree.h"
#include "../Core/interface/Plotting.h"
#include "../Core/interface/Plotting_Style.h"

void comparePtWeights7() {
  ModTDRStyle();
  TH1::SetDefaultSumw2(true);

  TFile f_cms("output/lhchcg/GluGluToHToTauTau_M-125_mt_2011_result.root");
  TTree *t_cms = (TTree*)f_cms.Get("pt_tree");

  TString draw = "pt>>cms_pt(60,0,300)";
  // TString draw = "TMath::Min(n_jets,2)>>cms_pt(3,-0.5,2.5)";
  t_cms->Draw(draw, "TMath::Min(n_jets,2)==2", "goff");
  TH1F *cms_pt = (TH1F*)gDirectory->Get("cms_pt");
  cms_pt->SetLineWidth(2);
  cms_pt->SetLineColor(kBlue+2);
  cms_pt->SetMarkerColor(kBlue+2);
  cms_pt->SetMarkerStyle(22);
  cms_pt->Scale(1. / cms_pt->Integral());


  TString draw_cms_old = "pt>>cms_pt_oldweighted(60,0,300)";
  // TString draw_cms_old = "TMath::Min(n_jets,2)>>cms_pt_oldweighted(3,-0.5,2.5)";
  t_cms->Draw(draw_cms_old, "cms_oldweight*(TMath::Min(n_jets,2)==2)", "goff");
  TH1F *cms_pt_oldweighted = (TH1F*)gDirectory->Get("cms_pt_oldweighted");
  cms_pt_oldweighted->SetLineWidth(2);
  cms_pt_oldweighted->SetLineColor(1);
  cms_pt_oldweighted->SetMarkerColor(1);
  cms_pt_oldweighted->SetMarkerStyle(26);
  cms_pt_oldweighted->Scale(1. / cms_pt_oldweighted->Integral());

  TString draw_cms = "pt>>cms_pt_weighted(60,0,300)";
  // TString draw_cms = "TMath::Min(n_jets,2)>>cms_pt_weighted(3,-0.5,2.5)";
  t_cms->Draw(draw_cms, "cms_weight*(TMath::Min(n_jets,2)==2)", "goff");
  TH1F *cms_pt_weighted = (TH1F*)gDirectory->Get("cms_pt_weighted");
  cms_pt_weighted->SetLineWidth(2);
  cms_pt_weighted->SetLineColor(kBlue);
  cms_pt_weighted->SetMarkerColor(kBlue);
  cms_pt_weighted->Scale(1. / cms_pt_weighted->Integral());

  TFile f_atlas("ATLAS_ggF_HpT_trees.root");
  TTree *t_atlas = (TTree*)f_atlas.Get("ggF_PowPy6_7TeV_nohfact");
  TString draw_atlas = "H_pT>>atlas_pt_weighted(60,0,300)";
  // TString draw_atlas = "TMath::Min(Njets25,2)>>atlas_pt_weighted(3,-0.5,2.5)";
  t_atlas->Draw(draw_atlas, "weight*(TMath::Min(Njets25,2)==2)", "goff");
  TH1F *atlas_pt_weighted = (TH1F*)gDirectory->Get("atlas_pt_weighted");
  atlas_pt_weighted->SetLineWidth(2);
  atlas_pt_weighted->SetLineColor(kRed);
  atlas_pt_weighted->SetMarkerColor(kRed);
  atlas_pt_weighted->SetMarkerStyle(24);
  atlas_pt_weighted->Scale(1. / atlas_pt_weighted->Integral());
  
  TCanvas * canv = new TCanvas("pt_ratio", "pt_ratio");
  canv->cd();
  std::vector<TPad*> pads = TwoPadSplit(0.29, 0.005, 0.005);

  pads[0]->SetLogy(true);

  std::vector<TH1*> h = CreateAxisHists(2, cms_pt_weighted);
  SetupTwoPadSplitAsRatio(pads, h[0], h[1], "CMS/ATLAS", true, 0.67, 1.33);

  h[0]->GetYaxis()->SetTitle("a.u.");
  h[1]->GetXaxis()->SetTitle("Higgs p_{T} (GeV)");
  // h[1]->GetXaxis()->SetTitle("Njets");
  h[0]->Draw("axis");
  h[0]->SetMinimum(1E-4);
  cms_pt->Draw("SAME");
  cms_pt_oldweighted->Draw("SAME");
  cms_pt_weighted->Draw("SAME");
  atlas_pt_weighted->Draw("SAME");
  FixTopRange(pads[0], GetPadYMax(pads[0]), 0.25);

  TLegend *legend = PositionedLegend(0.35, 0.25, 3, 0.03);
  legend->AddEntry(cms_pt, "CMS Unweighted");
  legend->AddEntry(cms_pt_oldweighted, "CMS Published");
  legend->AddEntry(cms_pt_weighted, "CMS Weighted");
  legend->AddEntry(atlas_pt_weighted, "ATLAS Weighted");
  FixBoxPadding(pads[0], legend, 0.05);
  legend->Draw();

  TH1F *ratio = reinterpret_cast<TH1F *>(
      MakeRatioHist(cms_pt_weighted, atlas_pt_weighted, true, true));
  // ratio->SetLineColor(kBlack);
  // ratio->SetMarkerColor(kBlack);

  TH1F *ratio_un = reinterpret_cast<TH1F *>(
      MakeRatioHist(cms_pt, atlas_pt_weighted, true, true));

  TH1F *ratio_old = reinterpret_cast<TH1F *>(
      MakeRatioHist(cms_pt_oldweighted, atlas_pt_weighted, true, true));

  pads[1]->cd();
  h[1]->Draw("axis");
  TLine line;
  line.SetLineStyle(2);
  DrawHorizontalLine(pads[1], &line, 1.0);
  ratio->Draw("esame");
  ratio_old->Draw("esame");
  ratio_un->Draw("esame");

  DrawCMSLogo(pads[0], "CMS", "Internal", 11, 0.045, 0.035, 1.2);
  DrawTitle(pads[0], "H#rightarrow#tau#tau: 2-jet ggH 7TeV", 1);


  pads[0]->RedrawAxis();
  pads[1]->RedrawAxis();


  canv->Print(".pdf");
}

