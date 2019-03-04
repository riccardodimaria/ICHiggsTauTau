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

int qcdStitch(){//main

  std::string outlabel = "STITCHPLOTS/FINAL";

  const unsigned nF = 2;
  TFile *f[nF];
  //f[0] = TFile::Open("/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170111/MET.root");
  f[0] = TFile::Open("/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170111/VBFQCD.root");
  f[1] = TFile::Open("/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170111/QCD.root");

  //std::string fileLabel[nF] = {"MET","VBFQCD","inclQCD"};
  std::string fileLabel[nF] = {"VBFQCD","inclQCD"};

  TLatex lat;

  const unsigned nV = 6;

  std::string var[nV] = {
    "lheHT",
    "jet1_eta",
    "jet2_eta",
    "jet1_pt",
    "jet2_pt",
    "ht30"
  };

  unsigned nbins[nV] = {120,100,100,100,100,100};
  double min[nV] = {300,-5,-5,0,0,0};
  double max[nV] = {1500,5,5,500,500,1500};

  TCanvas *myc[nV];
  for (unsigned iV(0); iV<nV;++iV){
    std::ostringstream label;
    label << "myc" << iV;
    myc[iV] = new TCanvas(label.str().c_str(),label.str().c_str(),500,1000);
    myc[iV]->Divide(1,3);
  }
  TH1F *h[nF][nV];
  TH1F *hSF[nV];
  TH1F *hSum[nV];

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  for (unsigned iV(0); iV<1;++iV){
    double maxy = 0;
    for (unsigned iF(0); iF<nF;++iF){//loop on files
      f[iF]->cd();
      TTree *tree = (TTree*)f[iF]->Get("LightTree");
      
      if (!tree) return 1;
      
      std::ostringstream label;
      label << "h_" << iF << "_" << iV;
      std::cout << label.str() << std::endl;
      h[iF][iV] = new TH1F(label.str().c_str(),(";"+var[iV]+";events").c_str(),nbins[iV],min[iV],max[iV]);
      h[iF][iV]->Sumw2();
      myc[iV]->cd(1);
      if (iF==0) {
	//recipe for VBFQCD
      tree->Draw((var[iV]+">>"+label.str()).c_str(),"(lheHT>300 && lheHT<=500 && dijet_M>1100 && jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 &&jet1_pt>50&& jet2_pt>50 && nvetomuons==0 && nvetoelectrons==0 && dijet_deta>3.2)*weight_nolepnotrig*1.0674","");
      tree->Draw((var[iV]+">>+"+label.str()).c_str(),"(lheHT>500 && lheHT<=550 && dijet_M>1100 && jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 &&jet1_pt>50&& jet2_pt>50 && nvetomuons==0 && nvetoelectrons==0 && dijet_deta>3.2)*weight_nolepnotrig*1.0674*63337753/(63337753+5027540)","");
      tree->Draw((var[iV]+">>+"+label.str()).c_str(),"(lheHT>550 && lheHT<=700 && dijet_M>1100 && jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 &&jet1_pt>50&& jet2_pt>50 && nvetomuons==0 && nvetoelectrons==0 && dijet_deta>3.2)*weight_nolepnotrig*(0.00074*lheHT+0.679)*(63337753/(63337753+5027540))","");
      }
      else {
	//recipe for inclusive
	//tree->Draw((var[iV]+">>"+label.str()).c_str(),"(lheHT>300 && lheHT<=500 && dijet_M>1100 && jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 &&jet1_pt>50&& jet2_pt>50 && nvetomuons==0 && nvetoelectrons==0 && dijet_deta>3.2)*weight_nolepnotrig","");
	tree->Draw((var[iV]+">>"+label.str()).c_str(),"(lheHT>500 && lheHT<=700 && dijet_M>1100 && jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 &&jet1_pt>50&& jet2_pt>50 && nvetomuons==0 && nvetoelectrons==0 && dijet_deta>3.2)*weight_nolepnotrig*(63337753/(63337753+5027540))","");
	tree->Draw((var[iV]+">>+"+label.str()).c_str(),"(lheHT>700 && dijet_M>1100 && jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 &&jet1_pt>50&& jet2_pt>50 && nvetomuons==0 && nvetoelectrons==0 && dijet_deta>3.2)*weight_nolepnotrig","");
      }
      if (h[iF][iV]->GetMaximum()>maxy)maxy=h[iF][iV]->GetMaximum();

      label.str("");
      label << "hSF_" << iV;
      hSF[iV] = new TH1F(label.str().c_str(),(";"+var[iV]+";events").c_str(),nbins[iV],min[iV],max[iV]);
      hSF[iV]->Sumw2();
      label.str("");
      label << "hSum_" << iV;
      hSum[iV] = new TH1F(label.str().c_str(),(";"+var[iV]+";events").c_str(),nbins[iV],min[iV],max[iV]);
      hSum[iV]->Sumw2();
      //lat.DrawLatexNDC(0.1,0.92,fileLabel[iF].c_str());
      //myc[iV]->Update();
    }//loop on files
    for (unsigned iF(0); iF<nF;++iF){//loop on files
      myc[iV]->cd(1);
      h[iF][iV]->SetLineColor(iF+1);
      h[iF][iV]->SetMarkerColor(iF+1);
      h[iF][iV]->SetMarkerStyle(iF+21);
      h[iF][iV]->SetFillColor(iF+1);
      h[iF][iV]->SetMaximum(1.1*maxy);
      h[iF][iV]->Draw(iF==0?"PE":"PEsame");
    }
    myc[iV]->cd(2);
    hSum[iV]->Add(h[0][iV],h[1][iV],1.,1.);
    hSum[iV]->Draw("PE");

    myc[iV]->cd(3);
    std::cout << " Integrals: incl " << h[1][iV]->Integral(0,h[1][iV]->GetNbinsX()+1) << " vbf " << h[0][iV]->Integral(0,h[0][iV]->GetNbinsX()+1) << std::endl;
    //hSF[iV]->Divide(h[1][iV],h[0][iV],1./h[1][iV]->Integral(0,h[1][iV]->GetNbinsX()+1),1./h[0][iV]->Integral(0,h[0][iV]->GetNbinsX()+1));
    hSF[iV]->Divide(h[1][iV],h[0][iV],1.,1.);
    hSF[iV]->Draw("PE");
    //hSF[iV]->Fit("pol0","","",300,550);
    //hSF[iV]->Fit("pol1","","",550,1500);
 
    myc[iV]->Update();
    if (var[iV].find("/")==var[iV].npos) myc[iV]->Print((outlabel+var[iV]+".pdf").c_str());
    else myc[iV]->Print((outlabel+var[iV].substr(0,var[iV].find("/"))+".pdf").c_str());
    //for (unsigned iF(0); iF<nF;++iF){//loop on files
    myc[iV]->cd(1);
    gPad->SetLogy(1);
    gPad->SetGridx(1);
    myc[iV]->cd(2);
    gPad->SetLogy(1);
    gPad->SetGridx(1);
    //}
    myc[iV]->Update();
    if (var[iV].find("/")==var[iV].npos) myc[iV]->Print((outlabel+var[iV]+"_log.pdf").c_str());
    else myc[iV]->Print((outlabel+var[iV].substr(0,var[iV].find("/"))+"_log.pdf").c_str());


  }//loop on vars

  return 0;
}//main

