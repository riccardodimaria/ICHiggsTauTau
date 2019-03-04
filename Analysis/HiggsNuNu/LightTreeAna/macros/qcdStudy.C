#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include <cmath>

#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TGaxis.h"

double mydeltaPhi(double phi1,double phi2)
{
  double result = phi1 - phi2;
  while (result > TMath::Pi()) result -= 2*TMath::Pi();
  while (result <= -TMath::Pi()) result += 2*TMath::Pi();
  return fabs(result);
};

double mydeltaR(double eta1,double phi1,double eta2,double phi2)
{
  double deta = eta1 - eta2;
  double dphi = mydeltaPhi(phi1, phi2);
  return std::sqrt(deta*deta + dphi*dphi);
};

int qcdStudy(){//main

  const unsigned nF = 3;
  TFile *f[nF];
  f[0] = TFile::Open("/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170111/MET.root");
  f[1] = TFile::Open("/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170111/ZJetsNuNu.root");
  f[2] = TFile::Open("/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170111/QCD.root");

  std::string fileLabel[nF] = {"MET","ZJetsNuNu","inclQCD"};

  TLatex lat;

  const unsigned nV = 24;

  std::string var[nV] = {
    //"(TMath::Sqrt(TMath::Power(jet1_pt*TMath::Cos(jet1_phi)+met_x,2)+TMath::Power(jet1_pt*TMath::Sin(jet1_phi)+met_y,2))-jet2_pt)/jet2_pt",
    //"(1+2*(met-TMath::Sqrt(TMath::Power(jet1_pt*TMath::Cos(jet1_phi)+jet2_pt*TMath::Cos(jet2_phi),2)+TMath::Power(jet1_pt*TMath::Sin(jet1_phi)+jet2_pt*TMath::Sin(jet2_phi),2)))/(met+TMath::Sqrt(TMath::Power(jet1_pt*TMath::Cos(jet1_phi)+jet2_pt*TMath::Cos(jet2_phi),2)+TMath::Power(jet1_pt*TMath::Sin(jet1_phi)+jet2_pt*TMath::Sin(jet2_phi),2))))/(1-2*(met-TMath::Sqrt(TMath::Power(jet1_pt*TMath::Cos(jet1_phi)+jet2_pt*TMath::Cos(jet2_phi),2)+TMath::Power(jet1_pt*TMath::Sin(jet1_phi)+jet2_pt*TMath::Sin(jet2_phi),2)))/(met+TMath::Sqrt(TMath::Power(jet1_pt*TMath::Cos(jet1_phi)+jet2_pt*TMath::Cos(jet2_phi),2)+TMath::Power(jet1_pt*TMath::Sin(jet1_phi)+jet2_p\t*TMath::Sin(jet2_phi),2))))",
    "unclustered_et*1./jet1_pt",
    "unclustered_et*1./met",
    "jetunclet_mindphi",
    "metnomuunclet_dphi",
    //"dijet_dphi",
    "(1+(jet1_pt-jet2_pt)/(jet1_pt+jet2_pt))/(1-(jet1_pt-jet2_pt)/(jet1_pt+jet2_pt))",
    "(1+(met-jet1_pt)/(met+jet1_pt))/(1-(met-jet1_pt)/(met+jet1_pt))",
    "(1+(unclustered_et+TMath::Sqrt(TMath::Power(jet1_pt*TMath::Cos(jet1_phi)+met_x,2)+TMath::Power(jet1_pt*TMath::Sin(jet1_phi)+met_y,2))-jet2_pt)/(jet2_pt+unclustered_et+TMath::Sqrt(TMath::Power(jet1_pt*TMath::Cos(jet1_phi)+met_x,2)+TMath::Power(jet1_pt*TMath::Sin(jet1_phi)+met_y,2))))/(1-(unclustered_et+TMath::Sqrt(TMath::Power(jet1_pt*TMath::Cos(jet1_phi)+met_x,2)+TMath::Power(jet1_pt*TMath::Sin(jet1_phi)+met_y,2))-jet2_pt)/(jet2_pt+unclustered_et+TMath::Sqrt(TMath::Power(jet1_pt*TMath::Cos(jet1_phi)+met_x,2)+TMath::Power(jet1_pt*TMath::Sin(jet1_phi)+met_y,2))))",
    "jet2_eta",
    "TMath::Sqrt(TMath::Power(jet1_pt*TMath::Cos(jet1_phi)+met_x,2)+TMath::Power(jet1_pt*TMath::Sin(jet1_phi)+met_y,2))",
    "TMath::Sqrt(TMath::Power(jet2_pt*TMath::Cos(jet2_phi)+met_x,2)+TMath::Power(jet2_pt*TMath::Sin(jet2_phi)+met_y,2))",
    "mydeltaPhi(jet1_phi,TMath::ACos(met_x/met))",
    "mydeltaPhi(jet2_phi,TMath::ACos(met_x/met))",
    "met",
    "calomet",
    "TMath::Abs(met-calomet)/met",
    "jet1_phi",
    "jet2_pt",
    "jet1_pt*1./jet1_uncorpt",
    "jet1_neutralhadfrac",
    "jet1_neutralemfrac",
    "jet1_chargedemfrac",
    "jet1_chargedhadfrac",
    "jet1_chargedmult",
    "jet1_neutralmult"
  };

  unsigned nbins[nV] = {100,100,100,100,200,200,200,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,20,30};
  double min[nV] = {0,0,0,0,-0.5,-0.5,-0.5,-5,0,0,0,0,0,0,0,-3.1416,0,0,0,0,0,0,0,0};
  double max[nV] = {2,2,3.1416,3.1416,5,3,3,5,500,500,3.1416,3.1416,500,500,2,3.1416,500,2,1,1,1,1,20,30};

  TCanvas *myc2D[nV];
  for (unsigned iV(0); iV<nV;++iV){
    std::ostringstream label;
    label << "myc" << iV;
    myc2D[iV] = new TCanvas(label.str().c_str(),label.str().c_str(),1500,500);
    myc2D[iV]->Divide(3,1);
  }
  TH2F *h2D[nF][nV];

  gStyle->SetOptStat(0);

  for (unsigned iV(0); iV<nV;++iV){
    for (unsigned iF(0); iF<nF;++iF){//loop on files
      f[iF]->cd();
      TTree *tree = (TTree*)f[iF]->Get("LightTree");
      
      if (!tree) return 1;
      
      std::ostringstream label;
      label << "h2D_" << iF << "_" << iV;
      h2D[iF][iV] = new TH2F(label.str().c_str(),(";#eta^{j1};"+var[iV]+";events").c_str(),100,-5,5,nbins[iV],min[iV],max[iV]);
      myc2D[iV]->cd(iF+1);
      tree->Draw((var[iV]+":jet1_eta>>"+label.str()).c_str(),"jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && dijet_M>800&&jet1_pt>40&&dijet_deta>3.2&& jet2_pt>40&&metnomuons>200 && alljetsmetnomu_mindphi<0.5 && nvetomuons==0 && nvetoelectrons==0 && n_jets_30==2","colz");
      
      lat.DrawLatexNDC(0.1,0.92,fileLabel[iF].c_str());
      myc2D[iV]->Update();
    }//loop on files
    if (var[iV].find("/")==var[iV].npos) myc2D[iV]->Print(("QCDPLOTS2j/"+var[iV]+".pdf").c_str());
    else myc2D[iV]->Print(("QCDPLOTS2j/"+var[iV].substr(0,var[iV].find("/"))+".pdf").c_str());
    for (unsigned iF(0); iF<nF;++iF){//loop on files
      myc2D[iV]->cd(iF+1);
      gPad->SetLogz(1);
      gPad->SetGridx(1);
    }
    myc2D[iV]->Update();
    if (var[iV].find("/")==var[iV].npos) myc2D[iV]->Print(("QCDPLOTS2j/"+var[iV]+"_log.pdf").c_str());
    else myc2D[iV]->Print(("QCDPLOTS2j/"+var[iV].substr(0,var[iV].find("/"))+"_log.pdf").c_str());


  }//loop on vars

  return 0;
}//main

