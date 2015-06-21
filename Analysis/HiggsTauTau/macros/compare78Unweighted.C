#include "TFile.h"
#include "TH2D.h"
#include "TTree.h"
#include "../Core/interface/Plotting.h"
#include "../Core/interface/Plotting_Style.h"

void comparePtWeights() {
  ModTDRStyle();
  TH1::SetDefaultSumw2(true);

  TFile f_cms7("output/lhchcg/GluGluToHToTauTau_M-125_mt_2011.root");
  TH1F *cms7_pt = (TH1F*)gDirectory->Get("higgs_pt/nominal");
  cms7_pt->SetLineWidth(2);
  cms7_pt->SetLineColor(kBlue+2);
  cms7_pt->SetMarkerColor(kBlue+2);
  cms7_pt->SetMarkerStyle(24);
  cms7_pt->SetMarkerSize(0.5);
  cms7_pt->Scale(1. / cms7_pt->Integral());

  TFile f_cms8("output/lhchcg/GluGluToHToTauTau_M-125_mt_2012.root");
  TH1F *cms8_pt = (TH1F*)gDirectory->Get("higgs_pt/nominal");
  cms8_pt->SetLineWidth(2);
  cms8_pt->SetLineColor(kBlue+2);
  cms8_pt->SetMarkerColor(kBlue+2);
  cms8_pt->SetMarkerStyle(20);
  cms8_pt->SetMarkerSize(0.5);
  cms8_pt->Scale(1. / cms8_pt->Integral());


  TFile f_atlas7("output/lhchcg/HpT_7TeV.root");
  TH1F *atlas7_pt_weighted = (TH1F*)gDirectory->Get("PowhegATLAS_7TeV");
  atlas7_pt_weighted->SetLineWidth(2);
  atlas7_pt_weighted->SetLineColor(kRed);
  atlas7_pt_weighted->SetMarkerColor(kRed);
  atlas7_pt_weighted->SetMarkerStyle(24);
  atlas7_pt_weighted->SetMarkerSize(0.5);
  atlas7_pt_weighted->Scale(1. / atlas7_pt_weighted->Integral());
  
  TFile f_atlas8("output/lhchcg/HiggsPtTree.root");
  TTree *t_atlas8 = (TTree*)f_atlas8.Get("HiggsPtTree");
  TString draw_atlas8 = "H_pT>>atlas8_pt_weighted(60,0,300)";
  t_atlas8->Draw(draw_atlas8, "", "goff");
  TH1F *atlas8_pt_weighted = (TH1F*)gDirectory->Get("atlas8_pt_weighted");
  atlas8_pt_weighted->SetLineWidth(2);
  atlas8_pt_weighted->SetLineColor(kRed);
  atlas8_pt_weighted->SetMarkerColor(kRed);
  atlas8_pt_weighted->SetMarkerStyle(20);
  atlas8_pt_weighted->SetMarkerSize(0.5);
  atlas8_pt_weighted->Scale(1. / atlas8_pt_weighted->Integral());
  
  TCanvas * canv = new TCanvas("pt_ratio", "pt_ratio");
  canv->cd();
  std::vector<TPad*> pads = TwoPadSplit(0.29, 0.005, 0.005);

  pads[0]->SetLogy(true);

  std::vector<TH1*> h = CreateAxisHists(2, cms7_pt);
  SetupTwoPadSplitAsRatio(pads, h[0], h[1], "CMS/ATLAS", true, 0.67, 1.33);

  h[0]->GetYaxis()->SetTitle("a.u.");
  h[1]->GetXaxis()->SetTitle("Higgs p_{T} (GeV)");
  //h[1]->GetXaxis()->SetTitle("Njets");
  h[0]->Draw("axis");
  h[0]->SetMinimum(1E-4);
  cms7_pt->Draw("SAME");
  cms8_pt->Draw("SAME");
  atlas7_pt_weighted->Draw("SAME");
  atlas8_pt_weighted->Draw("SAME");
  FixTopRange(pads[0], GetPadYMax(pads[0]), 0.25);

  TLegend *legend = PositionedLegend(0.35, 0.25, 3, 0.03);
  legend->AddEntry(cms7_pt, "CMS Unweighted 7TeV");
  legend->AddEntry(cms8_pt, "CMS Unweighted 8TeV");
  legend->AddEntry(atlas7_pt_weighted, "ATLAS Unweighted 7TeV");
  legend->AddEntry(atlas8_pt_weighted, "ATLAS Unweighted 8TeV");
  FixBoxPadding(pads[0], legend, 0.05);
  legend->Draw();

  TH1F *ratio7 = reinterpret_cast<TH1F *>(
      MakeRatioHist(cms7_pt, atlas7_pt_weighted, true, true));
  TH1F *ratio8 = reinterpret_cast<TH1F *>(
      MakeRatioHist(cms8_pt, atlas8_pt_weighted, true, true));
  pads[1]->cd();
  h[1]->Draw("axis");
  TLine line;
  line.SetLineStyle(2);
  DrawHorizontalLine(pads[1], &line, 1.0);
  ratio7->Draw("esame");
  ratio8->Draw("esame");

  DrawCMSLogo(pads[0], "CMS", "Internal", 11, 0.045, 0.035, 1.2);
  DrawTitle(pads[0], "H#rightarrow#tau#tau: Inclusive ggH", 1);


  pads[0]->RedrawAxis();
  pads[1]->RedrawAxis();


  canv->Print(".pdf");
}

