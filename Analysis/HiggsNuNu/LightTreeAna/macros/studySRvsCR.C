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


int studySRvsCR(){//main

  const unsigned nF = 4;
  TFile *f[nF];
  f[0] = TFile::Open("/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170327_loosencuts/ZJetsNuNu.root");
  f[1] = TFile::Open("/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170327_loosencuts/DYJetsToLL.root");
  f[2] = TFile::Open("/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170327_loosencuts/Wenu.root");
  f[3] = TFile::Open("/vols/cms/magnan/Hinvisible/RunIILT/output_lighttree_170327_loosencuts/Wmunu.root");

  std::string filename[nF] = {"Znunu","DY","Wenu","Wmunu"};

  const unsigned nP = 8;
  std::string fstate[nP] = {"ZeeSR","ZmumuSR","Zee","Zmumu","WeSR","WmuSR","We","Wmu"};

  const unsigned file[nP] = {0,0,1,1,2,3,2,3};

  std::string lweight[nP] = {
    "total_weight_lepveto",
    "total_weight_lepveto",
    "weight_nolepnotrig*weight_leptight",
    "total_weight_leptight",
    "total_weight_lepveto",
    "total_weight_lepveto",
    "weight_nolepnotrig*weight_leptight",
    "total_weight_leptight"
  };
  std::string lcutbase = "jet1_eta*jet2_eta<0 && abs(jet1_eta)<4.7 && abs(jet2_eta)<4.7 && dijet_M>1100&&jet1_pt>80&&dijet_deta>3.6&& jet2_pt>40 && nvetotaus==0 && nloosephotons==0 && dijet_dphi<1.5";

  std::string lcutlep[nP] = {
    " && nvetomuons==0 && nvetoelectrons==0 && (abs(calomet-met)/metnomuons)<0.5 && alljetsmetnomu_mindphi>0.5 && metnomuons>200",
    " && nvetomuons==0 && nvetoelectrons==0 && (abs(calomet-met)/metnomuons)<0.5 && alljetsmetnomu_mindphi>0.5 && metnomuons>200",
    " && nvetomuons==0 && nselelectrons>=1 && nvetoelectrons==2 && m_ee>60 && m_ee<120 && oppsign_ee && ele1_pt>40 && (abs(calomet-met)/metnoelectrons)<0.5 && alljetsmetnoel_mindphi>0.5 && m_ee_gen>60 && m_ee_gen<120 && metnoelectrons>200",
    " && nvetomuons==2 && nselmuons>=1 && nvetoelectrons==0 && m_mumu>60 && m_mumu<120 && oppsign_mumu && (abs(calomet-met)/metnomuons)<0.5 && alljetsmetnomu_mindphi>0.5 && m_mumu_gen>60 && m_mumu_gen<120 && metnomuons>200",
    " && nvetomuons==0 && nvetoelectrons==0 && (abs(calomet-met)/metnomuons)<0.5 && alljetsmetnomu_mindphi>0.5 && metnomuons>200",
    " && nvetomuons==0 && nvetoelectrons==0 && (abs(calomet-met)/metnomuons)<0.5 && alljetsmetnomu_mindphi>0.5 && metnomuons>200",
    " && nselelectrons==1 &&nvetomuons==0 &&nvetoelectrons==1 &&ele1_pt>40  && (abs(calomet-met)/metnoelectrons)<0.5 && alljetsmetnoel_mindphi>0.5 && metnoelectrons>200",
    " && nselmuons==1&&nvetomuons==1&&nvetoelectrons==0 && (abs(calomet-met)/metnomuons)<0.5 && alljetsmetnomu_mindphi>0.5 && metnomuons>200"
  };

  const unsigned nV = 1;
  std::string var[nV] = {
    /*    "jet1_pt",
    "jet2_pt",
    "jet1_eta",
    "jet2_eta",
    "dijet_M",
    "dijet_deta",
    "dijet_dphi",
    "metnomuons",
    "metnoelectrons"
    "jet1_E",                                                                                                                                  
    "jet2_E",
    "boson_pt",
    "gen_ele1_pt",
    "gen_ele1_eta",
    "gen_ele1_mindR_j1",
    "gen_ele1_mindR_j2",
    "gen_mu1_pt",
    "gen_mu1_eta",
    "gen_mu1_mindR_j1",
    "gen_mu1_mindR_j2",*/
    "n_vertices"
  };

  std::string title[nV] = {
    /*    ";p_{T}^{j1} (GeV);events",
    ";p_{T}^{j2} (GeV);events",
    ";#eta^{j1} (GeV);events",
    ";#eta^{j2} (GeV);events",
    ";M_{jj} (GeV);events",
    ";#Delta#eta_{jj};events",
    ";#Delta#phi_{jj};events",
    ";METno#mu;events",
    ";METnoe;events"
    ";E^{j1} (GeV);events",                                                                                                                 
    ";E^{j2} (GeV);events",
    ";boson p_{T}(GeV);events",
    ";p_{T}^{e_{gen}} (GeV); events",
    ";#eta^{e_{gen}}; events",
    ";min#Delta R(e_{gen},j1); events",
    ";min#Delta R(e_{gen},j2); events",
    ";p_{T}^{#mu_{gen}} (GeV); events",
    ";#eta^{#mu_{gen}}; events",
    ";min#DeltaR(#mu_{gen},j1); events",
    ";min#DeltaR(#mu_{gen},j2); events"*/
    ";n_{vtx};events"
  };

  //unsigned nbins[nV] = {50,50,50,50,50,50,50,50,50};
  //double binmin[nV] = {80,40,-5,-5,1100,3.6,0,200,200};
  //double binmax[nV] = {700,400,5,5,5000,9,1.5,600,600};
  //unsigned nbins[nV] = {50,50,50,50,50,50,50,50,50,50,50};
  //double binmin[nV] = {0,0,0,0,-5,0,0,0,-5,0,0};
  //double binmax[nV] = {5000,5000,600,400,5,10,10,400,5,10,10};
  unsigned nbins[nV] = {30};
  double binmin[nV] = {0};
  double binmax[nV] = {60};

  TCanvas *myc[nV];
  for (unsigned iC(0); iC<nV; ++iC){
    std::ostringstream label;
    label << "myc" << iC;
    myc[iC] = new TCanvas(label.str().c_str(),label.str().c_str(),1000,1000);
  }

  gStyle->SetOptStat(0);
  TH1F *hvar[nP][nV];

  TLatex lat;

  for (unsigned iP(0); iP<nP; ++iP){//loop on final state
      if (!f[file[iP]]) return 1;
      f[file[iP]]->cd();
      TTree *tree = (TTree*)gDirectory->Get("LightTree");
      if (!tree) return 1;

      for (unsigned iV(0); iV<nV; ++iV){//loop on vars
	myc[iV]->cd();
	std::ostringstream label;
	label.str("");
	label << "h_" << fstate[iP] << "_" << iV;
	hvar[iP][iV] = new TH1F(label.str().c_str(),title[iV].c_str(),nbins[iV],binmin[iV],binmax[iV]);
	hvar[iP][iV]->Sumw2();
	tree->Draw((var[iV]+">>"+label.str()).c_str(),(lweight[iP]+"*("+lcutbase+lcutlep[iP]+")").c_str());
	std::cout << fstate[iP] << " " << var[iV] << " " << hvar[iP][iV]->Integral() << std::endl;
	hvar[iP][iV]->Scale(1./hvar[iP][iV]->Integral());
	myc[iV]->Update();
      }//loop on vars
    
  }//loop on fs

  TLegend *leg = new TLegend(0.8,0.8,0.9,0.9);
  for (unsigned iV(0); iV<nV; ++iV){//loop on vars
    myc[iV]->Clear();
    myc[iV]->Divide(2,2);
    myc[iV]->cd(1);
    gPad->SetLogy(0);
    hvar[0][iV]->SetLineColor(1);
    hvar[0][iV]->SetMarkerColor(1);
    hvar[0][iV]->SetMarkerStyle(20);      
    hvar[0][iV]->Draw("PEL");
    if (iV==0) leg->AddEntry(hvar[0][iV],"SR","P");
    hvar[2][iV]->SetLineColor(2);
    hvar[2][iV]->SetMarkerColor(2);
    hvar[2][iV]->SetMarkerStyle(22);      
    hvar[2][iV]->Draw("PELsame");
    if (iV==0) leg->AddEntry(hvar[2][iV],"CR","P");
    leg->Draw("same");
    lat.DrawLatexNDC(0.5,0.93,"Zee");

    myc[iV]->cd(2);
    gPad->SetLogy(0);
    hvar[1][iV]->SetLineColor(1);
    hvar[1][iV]->SetMarkerColor(1);
    hvar[1][iV]->SetMarkerStyle(20);      
    hvar[1][iV]->Draw("PEL");
    hvar[3][iV]->SetLineColor(2);
    hvar[3][iV]->SetMarkerColor(2);
    hvar[3][iV]->SetMarkerStyle(22);      
    hvar[3][iV]->Draw("PELsame");
    leg->Draw("same");
    lat.DrawLatexNDC(0.5,0.93,"Z#mu#mu");

    myc[iV]->cd(3);
    gPad->SetLogy(0);
    hvar[4][iV]->SetLineColor(1);
    hvar[4][iV]->SetMarkerColor(1);
    hvar[4][iV]->SetMarkerStyle(20);      
    hvar[4][iV]->Draw("PEL");
    hvar[6][iV]->SetLineColor(2);
    hvar[6][iV]->SetMarkerColor(2);
    hvar[6][iV]->SetMarkerStyle(22);      
    hvar[6][iV]->Draw("PELsame");
    leg->Draw("same");
    lat.DrawLatexNDC(0.5,0.93,"We#nu");

    myc[iV]->cd(4);
    gPad->SetLogy(0);
    hvar[5][iV]->SetLineColor(1);
    hvar[5][iV]->SetMarkerColor(1);
    hvar[5][iV]->SetMarkerStyle(20);      
    hvar[5][iV]->Draw("PEL");
    hvar[7][iV]->SetLineColor(2);
    hvar[7][iV]->SetMarkerColor(2);
    hvar[7][iV]->SetMarkerStyle(22);      
    hvar[7][iV]->Draw("PELsame");
    leg->Draw("same");
    lat.DrawLatexNDC(0.5,0.93,"W#mu#nu");

    myc[iV]->Update();
    myc[iV]->Print(("SRCRPLOTS/"+var[iV]+".pdf").c_str());

    for (unsigned iC(0); iC<4; ++iC){
      myc[iV]->cd(iC+1);
      gPad->SetLogy(1);
    }
    myc[iV]->Update();
    myc[iV]->Print(("SRCRPLOTS/"+var[iV]+"_log.pdf").c_str());
    
  }

  return 1;
  //  
    

  /*  for (unsigned iP(0); iP<nP; ++iP){//loop on processes
    for (unsigned iS(0); iS<2; ++iS){//loop on updown
      
      TH1Fhjes[iP][iS] = new TH1F(label.str().c_str(),";p_{T}^{j2} (GeV);JESvar/central",100,0,500);
      hjes[iP][iS]->Sumw2();
      hjes[iP][iS]->Divide(hvar[iP][iS+1],hvar[iP][0],1,1);
      hjes[iP][iS]->SetLineColor(iS+2);
      hjes[iP][iS]->SetMarkerColor(iS+2);      
      hjes[iP][iS]->SetMarkerStyle(iP+20);      
      hjes[iP][iS]->Draw(iP==0&&iS==0?"PE":"PEsame");
    }
  }
  leg->Draw("same");
  myc[0]->Update();
  myc[0]->Print("JESPLOTS/jet2pt_jesunc.pdf");
  */

  return 0;

}//main
