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

int htStitch(){//main

  std::string plotdir="QCDPLOTSJetHT";
  const unsigned nF = 1;
  TFile *f[nF];
  f[0] = TFile::Open("/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170405/JetHT.root");

  std::string fileLabel[nF] = {"JetHT"};
  TLatex lat;

  const unsigned nV = 2;

  std::string var[nV] = {
    "jet1_pt",
    "ht30"
  };

  unsigned nbins[nV] = {100,100};
  double min[nV] = {0,0};
  double max[nV] = {700,1000};

  TCanvas *myc1D[nV];
  for (unsigned iV(0); iV<nV;++iV){
    std::ostringstream label;
    label << "myc" << iV;
    myc1D[iV] = new TCanvas(label.str().c_str(),label.str().c_str(),1500,500);
  }
  TH1F *h1Dqcd[nF][nV];
  gStyle->SetOptStat(0);

  unsigned Eht[11] = {125,200,250,300,350,400,475,600,650,800,900};
  //unsigned Ejet[10] = {40,60,80,140,200,260,320,400,450,500};

  for (unsigned iV(0); iV<nV;++iV){
    for (unsigned iF(0); iF<nF;++iF){//loop on files
      f[iF]->cd();
      TTree *tree = (TTree*)f[iF]->Get("LightTree");
      
      if (!tree) return 1;
      myc1D[iV]->cd();
      gPad->SetGridx(1);      
      gPad->SetGridy(1);      
      std::ostringstream label;
      label << "h1Dqcd_" << iF << "_" << iV;
      h1Dqcd[iF][iV] = new TH1F(label.str().c_str(),(";"+var[iV]+";events").c_str(),nbins[iV],min[iV],max[iV]);
      for (unsigned iB(0); iB<11;++iB){
	std::ostringstream lcutbase;
	lcutbase << "pass_httrigger_" << Eht[iB] << "*(jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && dijet_M>1100&&jet1_pt>80&&dijet_deta>3.6&& jet2_pt>40&& dijet_dphi<1.5 && nvetomuons==0 && nvetoelectrons==0 && nvetotaus==0 && n_jets_csv2medium==0 && pass_httrigger_" << Eht[iB] << ">0)";

	tree->Draw((var[iV]+">>+"+label.str()).c_str(),(lcutbase.str()).c_str(),"");
      }//loop on triggers

    h1Dqcd[iF][iV]->Draw("PE");
    myc1D[iV]->Update();
    myc1D[iV]->Print((plotdir+"/"+var[iV]+".pdf").c_str());
    }//loop on files
    
  }//loop on vars

  return 0;
}//main

