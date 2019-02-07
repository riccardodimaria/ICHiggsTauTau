#include <iostream>

#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTree.h"
#include "TStyle.h"

int plotKFactors(){

  bool doZ = true;

  TFile *oldF = TFile::Open("kfactors.root");
  TFile *newF = TFile::Open(doZ?"kfactor_VBF_zjets_v2.root":"kfactor_VBF_wjets_v2.root");

  TCanvas *myc = new TCanvas("myc","myc",1);
  TCanvas *myc2 = new TCanvas("myc2","myc2",1);

  TH1F *inclEWKQCD = 0;
  TH1F *vbfQCD = 0;
  TH1F *inclQCD = 0;
  TH1F *inclQCDold = 0;
  TH1F *vbfEWKQCD = 0;

  oldF->cd("EWKcorr");
  inclEWKQCD = (TH1F*)gDirectory->Get(doZ?"Z":"W");
  oldF->cd(doZ?"ZJets_012j_NLO":"WJets_012j_NLO");
  inclQCDold = (TH1F*)gDirectory->Get("nominal");
  oldF->cd(doZ?"ZJets_LO":"WJets_LO");
  inclEWKQCD->Divide(inclEWKQCD,(TH1F*)gDirectory->Get("inv_pt"),1,1);
  inclQCDold->Divide(inclQCDold,(TH1F*)gDirectory->Get("inv_pt"),1,1);

  newF->cd();
  vbfQCD = (TH1F*)gDirectory->Get("bosonPt_NLO_vbf");
  inclQCD = (TH1F*)gDirectory->Get("bosonPt_NLO_monojet");
  vbfQCD->Divide(vbfQCD,(TH1F*)gDirectory->Get("bosonPt_LO_vbf"),1,1);
  inclQCD->Divide(inclQCD,(TH1F*)gDirectory->Get("bosonPt_LO_monojet"),1,1);

  vbfEWKQCD = (TH1F*)vbfQCD->Clone("vbfEWKQCD");
  vbfEWKQCD->Divide(vbfQCD,inclQCD,1,1);

  for (int ibin(0); ibin<vbfEWKQCD->GetNbinsX();++ibin){
    double ptmin = vbfEWKQCD->GetXaxis()->GetBinLowEdge(ibin+1);
    double ptmax = vbfEWKQCD->GetXaxis()->GetBinLowEdge(ibin+2);
    double cor = 0;
    unsigned nbins = 0;

    std::cout << "Bin " << ibin << " " << ptmin << "-" << ptmax << std::endl;
    for (int jbin(0); jbin<inclEWKQCD->GetNbinsX();++jbin){

      double valmin = inclEWKQCD->GetXaxis()->GetBinLowEdge(jbin+1);
      double valmax = inclEWKQCD->GetXaxis()->GetBinLowEdge(jbin+2);
      if (valmin>=ptmin && valmax <= ptmax) {
	std::cout << " --- Adding " << valmin << "-" << valmax << " " << inclEWKQCD->GetBinContent(jbin+1) << std::endl;
	cor += inclEWKQCD->GetBinContent(jbin+1);
	nbins++;
      }
    }
    std::cout << " -> cor = " << cor << " nbins = " << nbins << std::endl;
    vbfEWKQCD->SetBinContent(ibin+1,vbfEWKQCD->GetBinContent(ibin+1)*cor/nbins);
  }
  //vbfEWKQCD->Multiply(vbfEWKQCD,vbfQCD,1,1);
  //vbfEWKQCD->Divide(vbfEWKQCD,inclQCD,1,1);

  myc->cd();
  inclQCD->SetLineColor(3);
  inclQCD->SetMarkerColor(3);
  inclQCD->SetMarkerStyle(20);
  inclQCD->GetYaxis()->SetRangeUser(0.5,2.5);
  inclQCD->Draw("PE");

  inclQCDold->SetLineColor(4);
  inclQCDold->SetMarkerColor(4);
  inclQCDold->SetMarkerStyle(21);
  inclQCDold->Draw("PEsame");

  vbfQCD->SetLineColor(2);
  vbfQCD->SetMarkerColor(2);
  vbfQCD->SetMarkerStyle(22);
  vbfQCD->Draw("PEsame");

  inclEWKQCD->SetLineColor(6);
  inclEWKQCD->SetMarkerColor(6);
  inclEWKQCD->SetMarkerStyle(23);
  inclEWKQCD->Draw("PEsame");

  vbfEWKQCD->SetLineColor(7);
  vbfEWKQCD->SetMarkerColor(7);
  vbfEWKQCD->SetMarkerStyle(24);
  vbfEWKQCD->Draw("PEsame");

  TLegend *leg = new TLegend(0.5,0.5,0.9,0.9);
  leg->AddEntry(inclQCDold,"QCD monojet old","P");
  leg->AddEntry(inclQCD,"QCD monojet new","P");
  leg->AddEntry(vbfQCD,"QCD vbf new","P");
  leg->AddEntry(inclEWKQCD,"EWKQCD monojet old","P");
  leg->AddEntry(vbfEWKQCD,"EWKQCD vbf new","P");
  leg->Draw("same");


  myc2->cd();

  TFile *test = TFile::Open("/vols/cms/rd1715/HiggsToInv/output_lighttree_170201_ICHEP/DYht.root");
  test->cd();

  TTree *tree = (TTree*)gDirectory->Get("LightTree");

  TH1F *bosonpt = (TH1F*)vbfEWKQCD->Clone("bosonpt");
  bosonpt->Reset();

  gStyle->SetOptStat(1111111);
  tree->Draw("boson_pt>>bosonpt","weight_trig_0*weight_leptight*weight_nolepnotrig*(jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && dijet_M>1100 && jet1_pt>80 && dijet_deta>3.6 && jet2_pt>40 && metnomuons>200 && nloosephotons==0 && alljetsmetnomu_mindphi>2.3 && nselmuons>=1&&nvetomuons==2 && nvetoelectrons==0&&m_mumu>60&&m_mumu<120&&oppsign_mumu)");
  bosonpt->Multiply(bosonpt,vbfEWKQCD,1,1);
  bosonpt->Draw("PE");

  return 0;
}
