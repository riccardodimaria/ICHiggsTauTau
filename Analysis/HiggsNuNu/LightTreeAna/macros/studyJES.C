#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TGraph.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>


int studyJES(){//main

  const unsigned nP = 2;
  const unsigned nF = 3;
  TFile *f[nP][nF];
  f[0][0] = TFile::Open("/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170327_loosencuts/ZJetsNuNu.root");
  f[0][1] = TFile::Open("/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170327_loosencuts/ZJetsNuNuJESUP.root");
  f[0][2] = TFile::Open("/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170327_loosencuts/ZJetsNuNuJESDOWN.root");
  f[1][0] = TFile::Open("/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170327_loosencuts/TT.root");
  f[1][1] = TFile::Open("/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170327_loosencuts/TTJESUP.root");
  f[1][2] = TFile::Open("/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170327_loosencuts/TTJESDOWN.root");

  std::string filename[nP] = {"ZNuNu","TT"};

  TH1F *hpt[nP][nF];
  TH1F *hjes[nP][2];

  TCanvas *myc = new TCanvas("myc","myc",1);
  myc->cd();
  gStyle->SetOptStat(0);

  std::string lcutbase = "weight_nolepnotrig*(jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && dijet_M>600&&jet1_pt>80&&dijet_deta>2.0&& jet2_pt>0 && nvetomuons==0 && nvetoelectrons==0 && nvetotaus==0)";

  std::ostringstream label;
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
  for (unsigned iP(0); iP<nP; ++iP){//loop on processes
    for (unsigned iF(0); iF<nF; ++iF){//loop on files
      if (!f[iP][iF]) return 1;
      f[iP][iF]->cd();
      TTree *tree = (TTree*)gDirectory->Get("LightTree");
      if (!tree) return 1;
      label.str("");
      label << "hpt_" << iP << "_" << iF;
      hpt[iP][iF] = new TH1F(label.str().c_str(),";p_{T}^{j2} (GeV);events",100,0,500);
      tree->Draw(("jet2_pt>>"+label.str()).c_str(),lcutbase.c_str(),(iP==0&&iF==0)?"PE":"PEsame");
      hpt[iP][iF]->SetLineColor(iF+1);
      hpt[iP][iF]->SetMarkerColor(iF+1);      
      hpt[iP][iF]->SetMarkerStyle(iP+20);      
    }//loop on files
    leg->AddEntry(hpt[iP][0],filename[iP].c_str(),"P");
    
    
  }//loop on processes
  leg->Draw("same");
  myc->Update();
  myc->Print("JESPLOTS/jet2pt_all.pdf");

  for (unsigned iP(0); iP<nP; ++iP){//loop on processes
    for (unsigned iS(0); iS<2; ++iS){//loop on updown
      
      label.str("");
      label << "hjes_" << iP << "_" << iS;
      hjes[iP][iS] = new TH1F(label.str().c_str(),";p_{T}^{j2} (GeV);JESvar/central",100,0,500);
      hjes[iP][iS]->Sumw2();
      hjes[iP][iS]->Divide(hpt[iP][iS+1],hpt[iP][0],1,1);
      hjes[iP][iS]->SetLineColor(iS+2);
      hjes[iP][iS]->SetMarkerColor(iS+2);      
      hjes[iP][iS]->SetMarkerStyle(iP+20);      
      hjes[iP][iS]->Draw(iP==0&&iS==0?"PE":"PEsame");
    }
  }
  leg->Draw("same");
  myc->Update();
  myc->Print("JESPLOTS/jet2pt_jesunc.pdf");


  return 0;

}//main
