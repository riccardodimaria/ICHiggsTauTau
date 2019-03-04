#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include <cmath>

#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFileCollection.h"
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

int qcdExtrap(){//main

  std::string plotdir="QCDPLOTSJetHTshape";
  //const bool plotRatio = true;

  const unsigned nF = 4;

  TChain *chain[nF];
  for (unsigned iF(0); iF<nF;++iF){//loop on files
    chain[iF] = new TChain("LightTree","");
  }

  TFileCollection fc0("dum","","JetHT.txt");
  chain[0]->AddFileInfoList((TCollection*)fc0.GetList());
  TFileCollection fc1("dum","","ZJetsToNuNu.txt");
  chain[1]->AddFileInfoList((TCollection*)fc1.GetList());
  TFileCollection fc2("dum","","WJets.txt");
  chain[2]->AddFileInfoList((TCollection*)fc2.GetList());
  TFileCollection fc3("dum","","QCD.txt");
  chain[3]->AddFileInfoList((TCollection*)fc3.GetList());

  std::string fileLabel[nF] = {"JetHT","ZJetsNuNu","WJets","inclQCD"};
  //std::string fileLabel[nF] = {"VBFQCD"};

  //jetmetdphi signal cut
  const unsigned nP = 6;
  double cutval[nP] = {0.5,0.75,1,1.25,1.5,2};

  TLatex lat;

  const unsigned nV = 3;

  std::string var[nV] = {
    "metnomuons",
    "dijet_M",
    "ht30"
  };

  unsigned nbins[nV] = {10,40,40};
  double min[nV] = {100,300,0};
  double max[nV] = {300,4300,2000};

  TCanvas *myc1D[nF];
  TCanvas *mycR[nF];
  for (unsigned iV(0); iV<nF;++iV){
    std::ostringstream label;
    label << "mycR" << iV;
    mycR[iV] = new TCanvas(label.str().c_str(),label.str().c_str(),1200,800);
    label.str("");
    label << "myc" << iV;
    myc1D[iV] = new TCanvas(label.str().c_str(),label.str().c_str(),1200,800);
    myc1D[iV]->Divide(2,1);
  }

  TFile *output = TFile::Open((plotdir+".root").c_str(),"RECREATE");
  output->cd();
  TH1F *h1Dqcd[nF][nV];
  TH1F *h1Dsig[nF][nP][nV];
  TH1F *htransfer[nF+1][nP][nV];
  gStyle->SetOptStat("eiuo");

  unsigned Eht[11] = {125,200,250,300,350,400,475,600,650,800,900};

  TLegend *leg = new TLegend(0.5,0.5,0.8,0.9);
  leg->SetFillStyle(0);

  for (unsigned iF(0); iF<nF;++iF){//loop on files
    std::string lcutbase = "weight_nolepnotrig*(ht30>500 && dijet_dphi<1.3 && metnomuons>=100 && dijet_deta>1 &&dijet_M>600&&jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7&&jet1_pt>80&& jet2_pt>40 && nvetomuons==0 && nvetoelectrons==0 && nvetotaus==0 && n_jets_csv2medium==0";
    if (iF==0 || iF==3) lcutbase = "ht30>500 && dijet_dphi<1.3 && metnomuons>=100 && dijet_M>600&&dijet_deta>1 &&jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7&&jet1_pt>80&& jet2_pt>40 && nvetomuons==0 && nvetoelectrons==0 && nvetotaus==0 && n_jets_csv2medium==0";
    
    //f[iF]->cd();
    //TTree *tree = (TTree*)f[iF]->Get("LightTree");
      
    if (!chain[iF]) return 1;

    for (unsigned iV(0); iV<nV;++iV){
      //Z and W+jets MC - no trigger req
      if (iF!=0 && iF!=3){
	std::ostringstream label;
	label << "h1Dqcd_" << iF << "_" << iV;
	h1Dqcd[iF][iV] = new TH1F(label.str().c_str(),(";"+var[iV]+";events").c_str(),nbins[iV],min[iV],max[iV]);
	myc1D[iF]->cd(1);
	gPad->SetLogy(1);
	chain[iF]->Draw((var[iV]+">>"+label.str()).c_str(),(lcutbase+"&& alljetsmetnomu_mindphi<0.5)").c_str(),"");
	h1Dqcd[iF][iV]->SetLineColor(1);      
	h1Dqcd[iF][iV]->Sumw2();

	for (unsigned iP(0); iP<nP; ++iP){
	  label.str("");
	  label << "h1Dsig_" << iF << "_" << iP << "_" << iV;
	  h1Dsig[iF][iP][iV] = new TH1F(label.str().c_str(),(";"+var[iV]+";events").c_str(),nbins[iV],min[iV],max[iV]);
	  h1Dsig[iF][iP][iV]->Sumw2();
	  myc1D[iF]->cd(2);
	  gPad->SetLogy(1);
	  std::ostringstream mindphicut;
	  mindphicut << "&& alljetsmetnomu_mindphi>" << cutval[iP];
	  chain[iF]->Draw((var[iV]+">>"+label.str()).c_str(),(lcutbase+mindphicut.str()+")").c_str(),iP==0?"":"same");
	  h1Dsig[iF][iP][iV]->SetLineColor(iP+1);
	  //h1Dsig[iF][iV]->Scale(1./h1Dsig[iF][iV]->GetEntries());
	  //h1Dqcd[iF][iV]->Scale(1./h1Dqcd[iF][iV]->GetEntries());
	}
      }
      else {
	std::ostringstream label;
	label.str("");
	label << "h1Dqcd_" << iF << "_" << iV;
	h1Dqcd[iF][iV] = new TH1F(label.str().c_str(),(";"+var[iV]+";events").c_str(),nbins[iV],min[iV],max[iV]);
	myc1D[iF]->cd(1);
	gPad->SetLogy(1);
	for (unsigned iB(0); iB<11;++iB){
	  std::ostringstream lbase;
	  lbase << "pass_httrigger_" << Eht[iB];
	  if (iF==3) lbase << "*weight_nolepnotrig";
	  lbase << "*(" << lcutbase << " && pass_httrigger_" << Eht[iB] << ">0";
	  chain[iF]->Draw((var[iV]+">>+"+label.str()).c_str(),(lbase.str()+"&& alljetsmetnomu_mindphi<0.5)").c_str(),"");
	}//loop on triggers
	h1Dqcd[iF][iV]->SetLineColor(1);      
	h1Dqcd[iF][iV]->Sumw2();
	h1Dqcd[iF][iV]->Draw();
	for (unsigned iP(0); iP<nP; ++iP){
	  label.str("");
	  label << "h1Dsig_" << iF << "_" << iP << "_" << iV;
	  h1Dsig[iF][iP][iV] = new TH1F(label.str().c_str(),(";"+var[iV]+";events").c_str(),nbins[iV],min[iV],max[iV]);
	  h1Dsig[iF][iP][iV]->Sumw2();
	  h1Dsig[iF][iP][iV]->SetLineColor(iP+1);      
	  myc1D[iF]->cd(2);
	  std::ostringstream mindphicut;
	  mindphicut << "&& alljetsmetnomu_mindphi>" << cutval[iP];
	  for (unsigned iB(0); iB<11;++iB){
	    std::ostringstream lbase;
	    lbase << "pass_httrigger_" << Eht[iB];
	    if (iF==3) lbase << "*weight_nolepnotrig";
	    lbase<< "*(" << lcutbase << " && pass_httrigger_" << Eht[iB] << ">0";
	    chain[iF]->Draw((var[iV]+">>+"+label.str()).c_str(),(lbase.str()+mindphicut.str()+")").c_str());
	  }//loop on triggers
	}
	myc1D[iF]->cd(2);
	gPad->SetLogy(1);
	for (unsigned iP(0); iP<nP; ++iP){
	  h1Dsig[iF][iP][iV]->Draw(iP==0?"":"same");      
	}
      }
      mycR[iF]->cd();
      gPad->SetGridx(1);      
      gPad->SetGridy(1);      

      for (unsigned iP(0); iP<nP; ++iP){
	htransfer[iF][iP][iV]=(TH1F*)h1Dsig[iF][iP][iV]->Clone();
	htransfer[iF][iP][iV]->SetStats(0);
	htransfer[iF][iP][iV]->GetYaxis()->SetTitle("trueMET/fakeMET");
	htransfer[iF][iP][iV]->Divide(h1Dsig[iF][iP][iV],h1Dqcd[iF][iV],1,1);
	htransfer[iF][iP][iV]->SetMarkerColor(1+iP);
	htransfer[iF][iP][iV]->SetLineColor(1+iP);
	htransfer[iF][iP][iV]->SetMarkerStyle(20+iP);
	htransfer[iF][iP][iV]->Draw(iP==0?"PE":"PEsame");
	std::ostringstream mindphicut;
	mindphicut << "min#Delta#phi(j,MET)>" << cutval[iP];
	if (iF==0 && iV==0) leg->AddEntry(htransfer[iF][iP][iV],mindphicut.str().c_str(),"P");
	//if (iV>0) htransfer[iF][iP][iV-1]->Draw("PEsame");
	if (iF!=1) htransfer[iF][iP][iV]->GetYaxis()->SetRangeUser(0.000001,1.5);
	//if (iF==0) htransfer[iF][iP][iV]->GetYaxis()->SetRangeUser(0,1.5);
      }
      //if (iF>0) {
      //htransfer[0][0][iV]->SetLineStyle(2);
      //htransfer[0][0][iV]->Draw("Lsame");
      //if (iF==1) leg->AddEntry(htransfer[0][0][iV],"JetHT Data","L");;
      //}
      leg->Draw("same");
      lat.DrawLatexNDC(0.1,0.92,fileLabel[iF].c_str());
      mycR[iF]->Update();
      mycR[iF]->Print((plotdir+"/"+fileLabel[iF]+"_"+var[iV]+"_transfer.pdf").c_str());
      mycR[iF]->Print((plotdir+"/"+fileLabel[iF]+"_"+var[iV]+"_transfer.C").c_str());

      myc1D[iF]->cd(1);
      lat.DrawLatexNDC(0.1,0.92,fileLabel[iF].c_str());
      myc1D[iF]->cd(2);
      leg->Draw("same");
      lat.DrawLatexNDC(0.1,0.92,fileLabel[iF].c_str());
      myc1D[iF]->Update();
      myc1D[iF]->Print((plotdir+"/"+fileLabel[iF]+"_"+var[iV]+".pdf").c_str());
      myc1D[iF]->Print((plotdir+"/"+fileLabel[iF]+"_"+var[iV]+".C").c_str());

    }//loop on vars
  }//loop on files


  for (unsigned iV(0); iV<nV;++iV){
    mycR[0]->cd();
    TH1F *htmp = (TH1F*)h1Dqcd[0][iV]->Clone();
    htmp->Add(h1Dqcd[1][iV],-1);
    htmp->Add(h1Dqcd[2][iV],-1);
    for (unsigned iP(0); iP<nP; ++iP){
      htransfer[nF][iP][iV]=(TH1F*)h1Dsig[0][iP][iV]->Clone();
      htransfer[nF][iP][iV]->Add(h1Dsig[1][iP][iV],-1);
      htransfer[nF][iP][iV]->Add(h1Dsig[2][iP][iV],-1);

      htransfer[nF][iP][iV]->SetStats(0);
      htransfer[nF][iP][iV]->GetYaxis()->SetTitle("trueMET/fakeMET");
      htransfer[nF][iP][iV]->Divide(htransfer[nF][iP][iV],htmp,1,1);
      htransfer[nF][iP][iV]->SetMarkerColor(1+iP);
      htransfer[nF][iP][iV]->SetLineColor(1+iP);
      htransfer[nF][iP][iV]->SetMarkerStyle(20+iP);
      htransfer[nF][iP][iV]->Draw(iP==0?"PE":"PEsame");
      std::ostringstream mindphicut;
      mindphicut << "min#Delta#phi(j,MET)>" << cutval[iP];
      //if (iV>0) htransfer[iF][iP][iV-1]->Draw("PEsame");
      htransfer[nF][iP][iV]->GetYaxis()->SetRangeUser(0.000001,1.5);
      //if (iF==0) htransfer[iF][iP][iV]->GetYaxis()->SetRangeUser(0,1.5);
    }
    //if (iF>0) {
    //htransfer[0][0][iV]->SetLineStyle(2);
    //htransfer[0][0][iV]->Draw("Lsame");
    //if (iF==1) leg->AddEntry(htransfer[0][0][iV],"JetHT Data","L");;
    //}
    leg->Draw("same");
    mycR[0]->Update();
    mycR[0]->Print((plotdir+"/Final_"+var[iV]+"_transfer.pdf").c_str());
    mycR[0]->Print((plotdir+"/Final_"+var[iV]+"_transfer.C").c_str());
    


  }//loop on vars



  return 0;
}//main

