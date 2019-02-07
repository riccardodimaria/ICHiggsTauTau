#include <iostream>

#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTree.h"
#include "TStyle.h"

int extractKFactors(){

  bool doZ = false;

  TFile *oldF = TFile::Open("kfactors.root");
  TFile *newF = TFile::Open(doZ?"kfactor_VBF_zjets_v2.root":"kfactor_VBF_wjets_v2.root");

  TFile *output = TFile::Open(doZ?"nloCorrections_Z.root":"nloCorrections_W.root","RECREATE");

  TH1F *kfactor_monojet_qcd_ewk = 0;
  TH1F *vbfQCD = 0;
  TH1F *inclQCD = 0;
  TH1F *vbf_cnc_modifier = 0;

  oldF->cd();
  kfactor_monojet_qcd_ewk = (TH1F*)(gDirectory->Get(doZ?"EWKcorr/Z":"EWKcorr/W")->Clone("kfactor_monojet_qcd_ewk"));
  kfactor_monojet_qcd_ewk->Divide(kfactor_monojet_qcd_ewk,(TH1F*)gDirectory->Get(doZ?"ZJets_LO/inv_pt":"WJets_LO/inv_pt"),1,1);

  newF->cd();
  vbfQCD = (TH1F*)gDirectory->Get("bosonPt_NLO_vbf");
  inclQCD = (TH1F*)gDirectory->Get("bosonPt_NLO_monojet");
  vbfQCD->Divide(vbfQCD,(TH1F*)gDirectory->Get("bosonPt_LO_vbf"),1,1);
  inclQCD->Divide(inclQCD,(TH1F*)gDirectory->Get("bosonPt_LO_monojet"),1,1);

  vbf_cnc_modifier = (TH1F*)vbfQCD->Clone("vbf_cnc_modifier");
  vbf_cnc_modifier->Divide(vbfQCD,inclQCD,1,1);

  TH1F *vbf_cnc_total = new TH1F("vbf_cnc_total",";boson pT (GeV); VBF/Monojet*QCDEWKMonojet",110,150,1250);


  for (int ibin(1); ibin<vbf_cnc_total->GetNbinsX()+1;++ibin){
    double ptval = vbf_cnc_total->GetXaxis()->GetBinCenter(ibin);
    int idx1 = vbf_cnc_modifier->FindBin(ptval);
    double check1l = vbf_cnc_modifier->GetXaxis()->GetBinLowEdge(idx1);
    double check1h = vbf_cnc_modifier->GetXaxis()->GetBinLowEdge(idx1+1);
    int idx2 = kfactor_monojet_qcd_ewk->FindBin(ptval);
    double check2l = kfactor_monojet_qcd_ewk->GetXaxis()->GetBinLowEdge(idx2);
    double check2h = kfactor_monojet_qcd_ewk->GetXaxis()->GetBinLowEdge(idx2+1);
    std::cout << ptval << " " << check1l << "-" << check1h << " " << check2l << "-" << check2h << std::endl;
    vbf_cnc_total->SetBinContent(ibin,vbf_cnc_modifier->GetBinContent(idx1)*kfactor_monojet_qcd_ewk->GetBinContent(idx2));
  }
  
  
  
  output->cd();
  kfactor_monojet_qcd_ewk->Write();
  vbf_cnc_modifier->Write();
  vbf_cnc_total->Write();

  return 0;
}
