#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/HistPlotter.h"
#include <iostream>
#include "TH1F.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include <map>
#include "TVectorD.h"
//#include "CommonTools/Utils/interface/TFileDirectory.h"
#include <boost/algorithm/string.hpp>
#include "TDirectory.h"
#include "TFile.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TASImage.h"
#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TROOT.h"

namespace ic{

  void setTDRStyle (){

    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasDefH(600);
    gStyle->SetCanvasDefW(600);
    gStyle->SetCanvasDefX(0);
    gStyle->SetCanvasDefY(0);

    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(0); 
    gStyle->SetPadGridX(0);
    gStyle->SetPadGridY(0);
    gStyle->SetGridColor(0);
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(1);

    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameBorderSize(1);
    gStyle->SetFrameFillColor(0);
    gStyle->SetFrameFillStyle(0);
    gStyle->SetFrameLineColor(1);
    gStyle->SetFrameLineStyle(1);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetHistLineColor(1);
    gStyle->SetHistLineStyle(0);
    gStyle->SetHistLineWidth(1);

    gStyle->SetEndErrorSize(2);
    gStyle->SetFuncColor(2);
    gStyle->SetFuncStyle(1);
    gStyle->SetFuncWidth(1);
    gStyle->SetOptDate(0);

    gStyle->SetOptFile(0);
    gStyle->SetOptStat(0);
    gStyle->SetStatColor(0); 
    gStyle->SetStatFont(42);
    gStyle->SetStatFontSize(0.04);
    gStyle->SetStatTextColor(1);
    gStyle->SetStatFormat("6.4g");
    gStyle->SetStatBorderSize(1);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.15);

    gStyle->SetPadTopMargin(0.07);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.05);

    gStyle->SetOptTitle(0);
    gStyle->SetTitleFont(42);
    gStyle->SetTitleColor(1);
    gStyle->SetTitleTextColor(1);
    gStyle->SetTitleFillColor(10);
    gStyle->SetTitleFontSize(0.05);

    gStyle->SetTitleColor(1, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetTitleYOffset(1.05);

    gStyle->SetLabelColor(1, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetLabelOffset(0.007, "XYZ");
    gStyle->SetLabelSize(0.04, "XYZ");

    gStyle->SetAxisColor(1, "XYZ");
    gStyle->SetStripDecimals(1); 
    gStyle->SetTickLength(0.025, "XYZ");
    gStyle->SetNdivisions(510, "XYZ");
    gStyle->SetPadTickX(1); 
    gStyle->SetPadTickY(1);

    gStyle->SetOptLogx(0);
    gStyle->SetOptLogy(0);
    gStyle->SetOptLogz(0);

    gStyle->SetPaperSize(20.,20.);
    gStyle->SetPaintTextFormat(".2f");

  }

  void HistPlotter::SetMCStackStyle(ic::LTPlotElement* ele) {
    ele->set_fill_color(ele->color());
    ele->set_fill_style(1001);
    ele->set_draw_fill(true);
    ele->set_draw_fill_in_legend(true);
    ele->set_draw_marker(false);
    ele->set_draw_line(false);
    ele->set_line_width(2);
    ele->set_draw_stat_error_y(false);
    ele->set_line_color(1);
    return;
  }
  void HistPlotter::SetSignalStyle(ic::LTPlotElement* ele) {
    ele->set_fill_style(1001);
    ele->set_draw_fill(true);
    ele->set_draw_marker(false);
    ele->set_draw_line(true);
    ele->set_draw_stat_error_y(false);
    ele->set_draw_fill_in_legend(false);
    ele->set_line_style(ele->line_style());
    ele->set_fill_color(0);
    ele->set_line_color(ele->color());
    ele->set_line_width(3);
    return;
  }
  void HistPlotter::SetDataStyle(ic::LTPlotElement* ele) {
    ele->set_line_color(kBlack);
    ele->set_marker_color(kBlack);
    ele->set_marker_size(1.2);
    ele->set_marker_style(20);
    ele->set_line_style(1);
    ele->set_fill_color(1);
    ele->set_fill_style(0);
    ele->set_draw_fill(false);
    ele->set_draw_fill_in_legend(false);
    ele->set_line_width(2);
    ele->set_draw_marker(true);
    ele->set_draw_line(true);
    ele->set_draw_stat_error_y(true);
    return;
  }

  HistPlotter::HistPlotter(std::string name) : LTModule(name){
    do_debug_=true;
    do_norm_=false;
    do_ratio_=false;
    do_ratio_line_=false;
    do_ratio_fitline_=false;
    add_underflows_=false;
    add_overflows_=false;
    add_preliminary_=true;
    outsuffix_="";
    toterror_ = 0;
  };

  HistPlotter::~HistPlotter(){ ;};

  int HistPlotter::Init(TFile* fs){
    fs_=fs;
    std::cout<<"Initialisation info for "<<module_name_<<":"<<std::endl;

    return 0;
  };

  int HistPlotter::Run(LTFiles* ){
    std::cout<<module_name_<<":"<<std::endl;
    TFile* file=fs_;

    gROOT->SetBatch(kTRUE);
    gROOT->ForceStyle(kTRUE);
    setTDRStyle ();
    gStyle->SetOptStat(0);

    //GET DIRECTORY TO WRITE TO
    TDirectory* stacksdir=0;
    TDirectory* numsdir=0;
    TDirectory* densdir=0;
    if(!fs_->GetDirectory("stacksdir")){
      stacksdir=file->mkdir("stacksdir");
      numsdir=file->mkdir("numsdir");
      densdir=file->mkdir("densdir");
    }
    else{
      int idir=2;
      while(fs_->GetDirectory(("stacksdir"+idir))){
        idir++;
      }
      stacksdir=file->mkdir(("stacksdir"+idir));
      numsdir=file->mkdir(("numsdir"+idir));
      densdir=file->mkdir(("densdir"+idir));
    }

    TDirectory* writedir=0;
    if(dirname_==""){
      writedir=file->mkdir("controlplots");
    }
    else if(!fs_->GetDirectory(dirname_.c_str())){
      writedir=file->mkdir(dirname_.c_str());
    }
    else{
      writedir=file->GetDirectory(dirname_.c_str());
    }
    writedir->cd();

    TCanvas* c1= new TCanvas("canvas","canvas", 600, 700);
    c1->SetTickx(1);
    c1->SetTicky(1);
    c1->cd();
    c1->SetRightMargin(0.06);


    double uppermin=0;
    if(do_ratio_) uppermin=0.3;
    else uppermin=0.;
    TPad* upper = new TPad("upper","pad",0, uppermin ,1 ,1);
    TPad* lower = new TPad("lower","pad",0, 0   ,1 ,uppermin);
    if(do_ratio_){
      upper->SetBottomMargin(0.03);
      lower->SetTopMargin(0.08);
      lower->SetBottomMargin(0.38);
      lower->Draw();
      lower->SetGridy();
    }
    else{
      upper->SetBottomMargin(0.17);
    }
    upper->Draw();
    c1->cd();
    //LOOP OVER ALL THE VARIABLES TO PLOT
    for(unsigned iShape=0;iShape<shapes_.size();iShape++){
      std::cout<<"  Drawing plot for "<<shapes_[iShape].name()<<std::endl;
      std::cout<<iShape<<std::endl;

      c1->SetTitle(shapes_[iShape].name().c_str());
      c1->SetName(shapes_[iShape].name().c_str());
      upper->SetLogy(shapes_[iShape].dology());

      THStack *stack=new THStack(("stack"+shapes_[iShape].name()).c_str(),("stacked "+shapes_[iShape].name()+" plots").c_str());
      bool stackempty=true;
      bool first=true;

      TH1F *tothisto=0;

      if(!do_ratio_) stack->SetTitle(shapes_[iShape].histtitle().c_str());

      //EXTRACT ALL THE HISTOS AND PUT THEM IN STACKED OR UNSTACKED GROUPS
      std::cout<<"    Getting histograms using "<< elements_.size() << " elements" << std::endl;
      for(unsigned iElement=0;iElement<elements_.size();iElement++){
        if(elements_[iElement].sample()==""){
          std::cout<<"HistPlotter ERROR: Element " << iElement << " with empty name, exiting with status 1"<<std::endl;
          std::cout<<"Ele before: " << elements_[iElement>0?iElement-1 : 0].sample() << " after " << elements_[iElement<(elements_.size()-1)?iElement+1:elements_.size()-1].sample() << std::endl;
          return 1;
        }
        if(!fs_->GetDirectory(elements_[iElement].sample().c_str())){
          std::cout<<"ERROR: No directory with name: "<<elements_[iElement].sample()<<std::endl;
          std::cout<<"Exiting with status 1"<<std::endl;
          return 1;
        }
        writedir->cd();
        TH1F* histo=0;
        if(file->Get((elements_[iElement].sample()+"/"+shapes_[iShape].name()).c_str())){
          histo =dynamic_cast<TH1F*>(file->Get((elements_[iElement].sample()+"/"+shapes_[iShape].name()).c_str()));
        }
        else{
          std::cout<<"histo " << elements_[iElement].sample() << "/" << shapes_[iShape].name() << " didn't open properly"<<std::endl;
          return 1;
        }
        writedir->cd();

        //Do blinding of shapes
        for(unsigned iblindshape=0;iblindshape<elements_[iElement].blindvar().size();iblindshape++){
          if(shapes_[iShape].name()==elements_[iElement].blindvar()[iblindshape]){
            for (int j = 0; j < histo->GetNbinsX()+2; ++j) {
              double low_edge = histo->GetBinLowEdge(j);
              double high_edge = histo->GetBinWidth(j)+histo->GetBinLowEdge(j);
              if (j==0 || j==histo->GetNbinsX()+1) std::cout << " Low edge - high edge " << low_edge << " " << high_edge << std::endl;
              double x_blind_min_=elements_[iElement].blindrange()[iblindshape].first;
              double x_blind_max_=elements_[iElement].blindrange()[iblindshape].second;
              if ((low_edge > x_blind_min_ && low_edge < x_blind_max_) || (high_edge > x_blind_min_ && high_edge < x_blind_max_)) {
          histo->SetBinContent(j,0);
          histo->SetBinError(j,0);
              }
            }
          }
        }

        if(do_norm_){
          if(Integral(histo)!=0)histo->Scale(1/(Integral(histo)));
        }

        if (add_underflows_ || add_overflows_){
          int tmpbins = histo->GetNbinsX();
          double tmpmin = histo->GetXaxis()->GetBinLowEdge(1);
          double tmpmax = histo->GetXaxis()->GetBinLowEdge(tmpbins+1);
          double tmpsize = (tmpmax-tmpmin)/tmpbins;
          if (add_underflows_) {
            tmpbins+=1;
            tmpmin = tmpmin-tmpsize;
          }
          if (add_overflows_&&shapes_[iShape].name()!="fourjetsmetnomu_mindphi"){
            tmpbins+=1;
            tmpmax = tmpmax+tmpsize;
          }
          std::ostringstream tmpname;
          tmpname << histo->GetName() << "_mod";
          TH1F *histTmp = new TH1F(tmpname.str().c_str(),histo->GetTitle(),tmpbins,tmpmin,tmpmax);
          histTmp->Sumw2();
          if (add_underflows_) {
            histTmp->SetBinContent(1,histo->GetBinContent(0));
            histTmp->SetBinError(1,histo->GetBinError(0));
            for (int ibin(1); ibin<histo->GetNbinsX()+1;++ibin){
              histTmp->SetBinContent(ibin+1,histo->GetBinContent(ibin));
              histTmp->SetBinError(ibin+1,histo->GetBinError(ibin));
            }
          }
          else {
            for (int ibin(1); ibin<histo->GetNbinsX()+1;++ibin){
              histTmp->SetBinContent(ibin,histo->GetBinContent(ibin));
              histTmp->SetBinError(ibin,histo->GetBinError(ibin));
            }
          }
          if (add_overflows_&&shapes_[iShape].name()!="fourjetsmetnomu_mindphi") {
            histTmp->SetBinContent(tmpbins,histo->GetBinContent(histo->GetNbinsX()+1));
            histTmp->SetBinError(tmpbins,histo->GetBinError(histo->GetNbinsX()+1));
          }

          elements_[iElement].set_hist_ptr(histTmp);
        }
        else elements_[iElement].set_hist_ptr(histo);


        //SETUP STYLE
        if(elements_[iElement].is_data()){
          SetDataStyle(&(elements_[iElement]));
        }
        else if(elements_[iElement].in_stack()){
          SetMCStackStyle(&(elements_[iElement]));
          stackempty=false;
        }
        else{
          SetSignalStyle(&(elements_[iElement]));
        }

        elements_[iElement].hist_ptr()->SetName(elements_[iElement].legname().c_str());
        elements_[iElement].hist_ptr()->Scale(elements_[iElement].scale());
        if(elements_[iElement].draw_marker()){
          elements_[iElement].hist_ptr()->SetMarkerColor(elements_[iElement].marker_color());
          elements_[iElement].hist_ptr()->SetMarkerStyle(elements_[iElement].marker_style());
          elements_[iElement].hist_ptr()->SetMarkerSize(elements_[iElement].marker_size());
          elements_[iElement].set_drawopts(elements_[iElement].drawopts()+"P");
          elements_[iElement].set_legopts(elements_[iElement].legopts()+"P");
        }
        //!!

        if(elements_[iElement].draw_fill()){
          elements_[iElement].hist_ptr()->SetFillColor(elements_[iElement].fill_color());
          elements_[iElement].hist_ptr()->SetFillStyle(elements_[iElement].fill_style());
          elements_[iElement].hist_ptr()->SetLineColor(elements_[iElement].line_color());
          elements_[iElement].set_drawopts(elements_[iElement].drawopts()+"hist");
          if (elements_[iElement].draw_fill_in_legend()) elements_[iElement].set_legopts(elements_[iElement].legopts()+"f");
          else elements_[iElement].set_legopts(elements_[iElement].legopts()+"l");
        }

        if(elements_[iElement].draw_line()){
          elements_[iElement].hist_ptr()->SetLineColor(elements_[iElement].line_color());
          elements_[iElement].hist_ptr()->SetLineStyle(elements_[iElement].line_style());
          elements_[iElement].hist_ptr()->SetLineWidth(elements_[iElement].line_width());
          elements_[iElement].set_legopts(elements_[iElement].legopts()+"l");
        }

        if(elements_[iElement].draw_stat_error_y()==true){
          elements_[iElement].set_drawopts(elements_[iElement].drawopts()+"E1");
          elements_[iElement].set_legopts(elements_[iElement].legopts()+"e");
        }

        //ADD STACKED HISTOS TO STACK
        if((!elements_[iElement].is_data())&&elements_[iElement].in_stack()){
          stack->Add(elements_[iElement].hist_ptr());
          if (first) {
            tothisto=(TH1F*)(elements_[iElement].hist_ptr())->Clone(("totalerror"+shapes_[iShape].name()).c_str());
            first=false;
          }
          else tothisto->Add(elements_[iElement].hist_ptr());
        }
      }
      stacksdir->cd();
      stack->Write();

      //set total error // There was a bug here, producing error bars on MC on the same size of the yield
      if (tothisto && toterror_>0){
        for (int iB(1);iB<tothisto->GetNbinsX()+1;++iB){
          tothisto->SetBinError(iB,toterror_*tothisto->GetBinContent(iB));
                std::cout << " ---- DEBUG bin " << iB << " bin content: " << tothisto->GetBinContent(iB) << " error: " << tothisto->GetBinError(iB) << std::endl;
        }
      }

      //Get maximum
      first=true;
      double ymax=0;
      if(!stackempty) ymax=stack->GetMaximum();
      for(unsigned iElement=0;iElement<elements_.size();iElement++){
        if(elements_[iElement].hist_ptr()->GetMaximum()>ymax) ymax=elements_[iElement].hist_ptr()->GetMaximum();
      }
      upper->cd();
      //DRAW THE STACK AND ALL THE UNSTACKED HISTOS INCLUDING THE DATA
      if(!stackempty){ 
        std::cout<<"    Drawing Stack.."<<std::endl;
        if(first){
          if(ymax>=1){
            stack->SetMaximum(shapes_[iShape].axisrangemultiplier()*(ymax+sqrt(ymax)));
          }else{
            stack->SetMaximum(shapes_[iShape].axisrangemultiplier()*(ymax));
          }
          upper->cd();
          stack->Draw("hist");

          first=false;
          stack->GetYaxis()->SetLabelSize(0.04);

          stack->GetYaxis()->SetTitleSize(0.05);
          stack->GetYaxis()->SetTitleOffset(1.15);
          if(do_ratio_){
            stack->GetXaxis()->SetLabelOffset(999);
            stack->GetXaxis()->SetLabelSize(0);
          }else{
           stack->GetXaxis()->SetLabelSize(0.04);
           stack->GetXaxis()->SetTitleSize(0.05);
           stack->GetXaxis()->SetTitleOffset(0.9);

           std::string xtitle;
           xtitle=shapes_[iShape].histtitle().substr(shapes_[iShape].histtitle().find(";")+1);
           xtitle=xtitle.substr(0,xtitle.find(";"));
           stack->GetXaxis()->SetTitle(xtitle.c_str());
           stack->GetYaxis()->SetTitle("Events");
          }
          std::string ytitle;
          ytitle=shapes_[iShape].histtitle().substr(shapes_[iShape].histtitle().find(";")+1);
          ytitle=ytitle.substr(ytitle.find(";")+1);
          ytitle=ytitle.substr(0,ytitle.find(";"));
          if(do_ratio_) stack->SetTitle((";;"+ytitle).c_str());
          else{
            stack->SetTitle(shapes_[iShape].histtitle().c_str());
          }
        }else stack->Draw("histsame");
      }

      //drawing total stack with error band
      if (tothisto){
        tothisto->SetFillStyle(3004);
        tothisto->SetFillColor(1);
        tothisto->SetLineColor(1);
        tothisto->Draw("E2same");

        TH1F *tothisto2 = (TH1F*)tothisto->Clone(("totline"+shapes_[iShape].name()).c_str());
        tothisto2->SetFillStyle(0);
        tothisto2->Draw("histsame");
      }

      std::cout<<"    Drawing Unstacked.."<<std::endl;
      TGraphAsymmErrors * g = 0;
      for(unsigned iElement=0;iElement<elements_.size();iElement++){
        if(do_debug_)std::cout<<"  "<<elements_[iElement].hist_ptr()->GetName()<<std::endl;
        if(!(elements_[iElement].in_stack())){
          if(!(elements_[iElement].is_data()&&!(elements_[iElement].is_trigeff()))){
            if(first){
              if(!do_ratio_) elements_[iElement].hist_ptr()->SetTitle(shapes_[iShape].histtitle().c_str());
              elements_[iElement].hist_ptr()->Draw(elements_[iElement].drawopts().c_str());
              if(ymax>=1){
                elements_[iElement].hist_ptr()->GetYaxis()->SetRangeUser(0,shapes_[iShape].axisrangemultiplier()*(ymax+sqrt(ymax)));
              }else{
                elements_[iElement].hist_ptr()->GetYaxis()->SetRangeUser(0,shapes_[iShape].axisrangemultiplier()*(ymax));
              }
              //!!	      elements_[iElement].hist_ptr()->Draw(elements_[iElement].drawopts().c_str());
              //c1->Update();
              first=false;
              elements_[iElement].hist_ptr()->GetYaxis()->SetLabelSize(0.04);
              elements_[iElement].hist_ptr()->GetYaxis()->SetTitleSize(0.05);

              if(do_ratio_){
                elements_[iElement].hist_ptr()->GetXaxis()->SetLabelOffset(999);
                elements_[iElement].hist_ptr()->GetXaxis()->SetLabelSize(0);
              }else{
                elements_[iElement].hist_ptr()->GetXaxis()->SetLabelSize(0.04);
                elements_[iElement].hist_ptr()->GetXaxis()->SetTitleSize(0.05);
                elements_[iElement].hist_ptr()->GetXaxis()->SetTitleOffset(0.9);
                std::string xtitle;
                xtitle=shapes_[iShape].histtitle().substr(shapes_[iShape].histtitle().find(";")+1);
                xtitle=xtitle.substr(0,xtitle.find(";"));
                elements_[iElement].hist_ptr()->GetXaxis()->SetTitle(xtitle.c_str());
              }
              std::string ytitle;
              ytitle=shapes_[iShape].histtitle().substr(shapes_[iShape].histtitle().find(";")+1);
              ytitle=ytitle.substr(ytitle.find(";")+1);
              ytitle=ytitle.substr(0,ytitle.find(";"));
              elements_[iElement].hist_ptr()->GetYaxis()->SetTitle(ytitle.c_str());
              elements_[iElement].hist_ptr()->SetTitle("");
            }else{
              elements_[iElement].hist_ptr()->Draw(("same"+elements_[iElement].drawopts()).c_str());
            }
          }
          //!!
          else{//For data do asymmetric data errors
            const double alpha = 1 - 0.6827;
            g = new TGraphAsymmErrors(elements_[iElement].hist_ptr());
            g->SetMarkerSize(1.1);
            g->SetMarkerStyle (20);

            for (int i = 0; i < g->GetN(); ++i) {
              int N = g->GetY()[i];
              double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
              double U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
              g->SetPointEYlow(i, N-L);
              g->SetPointEYhigh(i, U-N);
            }
            if(first){
              if(!do_ratio_) g->SetTitle(shapes_[iShape].histtitle().c_str());
              g->Draw("AP");
              g->Draw(elements_[iElement].drawopts().c_str());
              g->GetYaxis()->SetRangeUser(0.,shapes_[iShape].axisrangemultiplier()*(ymax+sqrt(ymax)+1));
              g->Draw(elements_[iElement].drawopts().c_str());

              first=false;
              g->GetYaxis()->SetLabelSize(0.04);
              g->GetYaxis()->SetTitleSize(0.05);

              if(do_ratio_){
                g->GetXaxis()->SetLabelOffset(999);
                g->GetXaxis()->SetLabelSize(0);
              }else{
                g->GetXaxis()->SetLabelSize(0.04);
                g->GetXaxis()->SetTitleSize(0.05);
                g->GetXaxis()->SetTitleOffset(0.9);
                std::string xtitle;
                xtitle=shapes_[iShape].histtitle().substr(shapes_[iShape].histtitle().find(";")+1);
                xtitle=xtitle.substr(0,xtitle.find(";"));
                g->GetXaxis()->SetTitle(xtitle.c_str());
              }
              std::string ytitle;
              ytitle=shapes_[iShape].histtitle().substr(shapes_[iShape].histtitle().find(";")+1);
              ytitle=ytitle.substr(ytitle.find(";")+1);
              ytitle=ytitle.substr(0,ytitle.find(";"));
              g->SetTitle(ytitle.c_str());
            }
            else g->Draw("sameP");
          }
        }
      }
      if ( g != 0 ){g->Draw("sameP");}
      else {
        std::cout<<" -- g not defined! -- ISSUE "<<std::endl;
        std::cout<<" -- The element size is: " << elements_.size() <<std::endl;
      }
      if(do_debug_)std::cout<<"  Drawing legend"<<std::endl;

      //SETUP AND DRAW THE LEGEND
      double legleft=shapes_[iShape].legleft();
      double legright=shapes_[iShape].legright();
//       double legbottom=0.3+((10-elements_.size())*(0.89-0.3)/10);
//       if(legbottom<0.3||legbottom>0.89)legbottom=0.3;

      TLegend* leg=new TLegend(legleft,0.55,legright,0.9);
      leg->SetName("thelegend");
      leg->SetFillStyle(0);
      leg->SetFillColor(0);
      leg->SetBorderSize(0);
      for(unsigned iElement=0;iElement<elements_.size();iElement++){
        if(do_debug_)std::cout<<"  Adding legend entry for: "<<elements_[iElement].hist_ptr()->GetName()<<std::endl;
        leg->AddEntry(elements_[iElement].hist_ptr(),elements_[iElement].hist_ptr()->GetName(),elements_[iElement].legopts().c_str());
      }
      if(do_debug_)std::cout<<"  All entries added"<<std::endl;

      leg->Draw("SAME");


      if(do_debug_)std::cout<<"  Drawing CMS Logo"<<std::endl;
      if(add_preliminary_){
        DrawCMSLogoTest(upper,"CMS","Preliminary",10);
      } else{
        DrawCMSLogoTest(upper,"CMS","",10);
      }
      upper->cd();

      //DRAW RATIO PLOT
      //!!
      if(do_ratio_){
        std::cout<<"  Drawing ratio plot"<<std::endl;
        c1->cd();
        lower->cd();

        bool firstnum=true;
        bool firstden=true;
        TH1F* num = 0;
        TH1F* den = 0;
        double dentoterr=0;

        for(unsigned iElement=0;iElement<elements_.size();iElement++){
          std::cout<<elements_[iElement].hist_ptr()->GetName()<<std::endl;


          if(elements_[iElement].is_inrationum()){
            //ADD TO num HIST
            if(firstnum){
              num=(TH1F*)(elements_[iElement].hist_ptr()->Clone(("num"+shapes_[iShape].name()).c_str()));
              firstnum=false;
            }else{
              num->Add(elements_[iElement].hist_ptr());
            }
          }

          if(elements_[iElement].is_inratioden()){
            //ADD TO den HIST
            if(firstden){
              den=(TH1F*)(elements_[iElement].hist_ptr()->Clone(("den"+shapes_[iShape].name()).c_str()));
              firstden=false;
            }else{
              den->Add(elements_[iElement].hist_ptr());
            }
            //get error on this contribution
            double thiselementfracerr;
            double thiselementintegral=Integral(elements_[iElement].hist_ptr());
            if(elements_[iElement].has_dderrors()==1){
              TVectorD* dennormerr =(TVectorD*)file->Get((elements_[iElement].sample()+"/normerrs").c_str());//!!
              thiselementfracerr=sqrt(pow((*dennormerr)[0],2)+pow((*dennormerr)[1],2));
            }else{
              thiselementfracerr=Error(elements_[iElement].hist_ptr())/thiselementintegral;
            }
            if(thiselementintegral!=0){
              dentoterr=sqrt(pow(dentoterr,2)+pow(thiselementfracerr*thiselementintegral,2));
            }
            //std::cout<<dentoterr<<std::endl
          }
        }
        numsdir->cd();
        num->Write();
        densdir->cd();
        den->Write();
        lower->cd();

        if(firstnum||firstden)std::cout<<"To draw ratio plot you must specify elements to be numerator and denominator! Ratio plot will be missing."<<std::endl;
        else{
          //Set den error to zero will take den error into account in error band!!
          double denfracerr=dentoterr/Integral(den);
          double numfracerr=1/sqrt(Integral(num));
          //std::cout<<Integral(den);
          //std::cout<<denfracerr<<std::endl;
          TH1F* errorband=(TH1F*)(den->Clone("errorband"));
          for(int bin=0;bin<=den->GetNbinsX()+1;bin++){
            den->SetBinError(bin,0);
            errorband->SetBinContent(bin,Integral(num)/Integral(den));
            errorband->SetBinError(bin,sqrt(pow(denfracerr,2)+pow(numfracerr,2)));
          }

          //DIVIDE NUM BY DEN and put in ratio
          TH1F* ratio=0;
          ratio=(TH1F*)(num->Clone("ratio"));
          ratio->GetYaxis()->SetNdivisions(505);
          ratio->GetYaxis()->SetNdivisions(505);

          ratio->GetXaxis()->SetLabelSize(0.1);
          ratio->GetXaxis()->SetLabelOffset(0.005);

          ratio->GetXaxis()->SetTitleOffset(1.20);
          ratio->GetXaxis()->SetTitleSize(0.12);
          ratio->GetYaxis()->SetTitleSize(0.1);
          ratio->GetYaxis()->SetTitleOffset(0.5);
          ratio->GetYaxis()->SetLabelSize(0.09);

          errorband->GetYaxis()->SetNdivisions(505);
          errorband->GetXaxis()->SetLabelSize(0.1);

          errorband->GetXaxis()->SetTitleOffset(1.20);
          errorband->GetXaxis()->SetTitleSize(0.12);
          errorband->GetYaxis()->SetTitleSize(0.1);
          errorband->GetYaxis()->SetTitleOffset(0.5);
          errorband->GetYaxis()->SetLabelSize(0.09);


          std::string xtitle;
          xtitle=shapes_[iShape].histtitle().substr(shapes_[iShape].histtitle().find(";")+1);
          xtitle=xtitle.substr(0,xtitle.find(";"));
          ratio->GetXaxis()->SetTitle(xtitle.c_str());//!!GET TITLE FOR X AXIS
          //	  ratio->GetXaxis()->SetTitleSize(0.1);
          //ratio->GetXaxis()->SetTitleOffset(0.8);
          //std::cout<<ratio->GetMaximum();
          double ratiomax=ratio->GetMaximum()+0.1;
          if(ratiomax<1.5 || ratiomax>3) {
            std::cout << " Found ratio max = " << ratiomax << " setting to 2." << std::endl;
            ratiomax=1.5;
          }
          ratio->GetYaxis()->SetRangeUser(0.5,ratiomax);
          ratio->SetTitle("");
          ratio->GetYaxis()->SetTitle("Data/Bkg   ");
          ratio->GetYaxis()->SetTitleSize(0.1);

          ratio->Divide(den);

          ratio->SetStats(0);

          errorband->SetMarkerSize(0);
          errorband->SetFillColor(16);
          errorband->SetFillStyle(1000);//3013);
          errorband->SetLineWidth(1);
          errorband->GetXaxis()->SetTitle(xtitle.c_str());//!!GET TITLE FOR X AXIS
          //	  errorband->GetXaxis()->SetTitleSize(0.1);
          //errorband->GetXaxis()->SetTitleOffset(0.8);
          errorband->GetYaxis()->SetRangeUser(0.5,1.5);
          errorband->SetTitle("");
          errorband->GetYaxis()->SetTitle("data/MC");
          //errorband->GetYaxis()->SetTitleSize(0.1);
          //errorband->GetYaxis()->SetTitleOffset(0.3);
          //errorband->GetXaxis()->SetLabelSize(0.1);
          //errorband->GetYaxis()->SetLabelSize(0.1);
          errorband->SetStats(0);



      //	  TF1* fiterrup=0;
      //	  TF1* fiterrdown=0;
          TLine *averageLine=0;
          TLine *averageupLine=0;
          TLine *averagedownLine=0;
          if(do_ratio_fitline_){
            //ratio->Fit("pol0","E");
            averageLine = new TLine(ratio->GetXaxis()->GetBinLowEdge(1),Integral(num)/Integral(den),ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1),Integral(num)/Integral(den));
            averageupLine = new TLine(ratio->GetXaxis()->GetBinLowEdge(1),(Integral(num)+sqrt(Integral(num)))/Integral(den),ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1),(Integral(num)+sqrt(Integral(num)))/Integral(den));
                  averagedownLine = new TLine(ratio->GetXaxis()->GetBinLowEdge(1),(Integral(num)-sqrt(Integral(num)))/Integral(den),ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1),(Integral(num)-sqrt(Integral(num)))/Integral(den));
            averageLine->SetLineColor(1);
            averageLine->SetLineStyle(2);
            averageupLine->SetLineColor(3);
            averagedownLine->SetLineColor(3);
            averageupLine->SetLineStyle(3);
            averagedownLine->SetLineStyle(3);
            // fiterrup=new TF1("pol0up","pol0");
      // 	    TF1 *fitresult=ratio->GetFunction("pol0");
      // 	    //fiterrup->SetParameters(fitresult->GetParameter(0)+fitresult->GetParError(0),fitresult->GetParameter(1)+fitresult->GetParError(1));
      // 	    fiterrup->SetParameter(0,fitresult->GetParameter(0)+fitresult->GetParError(0));                                                             
      // 	    fiterrdown=new TF1("pol0up","pol0");
      // 	    //fiterrdown->SetParameters(fitresult->GetParameter(0)-fitresult->GetParError(0),fitresult->GetParameter(1)-fitresult->GetParError(1));
      // 	    fiterrdown->SetParameter(0,fitresult->GetParameter(0)-fitresult->GetParError(0));                                                            
      // 	    fiterrup->SetLineStyle(2);
      // 	    fiterrdown->SetLineStyle(2);
      // 	    fiterrup->SetLineColor(2);
      // 	    fiterrdown->SetLineColor(2);
          }
          if(do_ratio_fitline_){
            //fiterrup->Draw("same");
            //fiterrdown->Draw("same");
            errorband->Draw("E2");
            averageLine->Draw("same");
            ratio->Draw("E1same");
            //averageupLine->Draw("same");
            //averagedownLine->Draw("same");
          }else ratio->Draw("E1");

          TLine *centerLine = new TLine(ratio->GetXaxis()->GetBinLowEdge(1),1.0,ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1),1.0);
          centerLine->SetLineWidth(1);
          if(!do_ratio_fitline_){
//             centerLine->SetLineColor(2);
            centerLine->SetLineColor(1);
            centerLine->SetLineStyle(2);
          }
          else centerLine->SetLineColor(1);
          centerLine->Draw("same");
          if(do_ratio_line_){
            TLine *lowerLine = new TLine(ratio->GetXaxis()->GetBinLowEdge(1),0.9,ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1),0.9);
            TLine *upperLine = new TLine(ratio->GetXaxis()->GetBinLowEdge(1),1.1,ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1),1.1);
            lowerLine->SetLineWidth(2);
            upperLine->SetLineWidth(2);
            lowerLine->SetLineColor(7);
            upperLine->SetLineColor(7);
            lowerLine->Draw();
            upperLine->Draw();
          }
        }
      }
      //save as PDF
      //c1->Update();
      std::ostringstream lsave;
//       std::ostringstream lsavepng;
//       std::ostringstream lsavec;
      std::string tmpstr = file->GetName();
      tmpstr.erase(std::string(file->GetName()).find(".root"),5);
      tmpstr=tmpstr+outsuffix_;
      lsave << tmpstr ;
      if (shapes_[iShape].dology()) lsave << "_logy";
      lsave << ".pdf" ;
      if (iShape==0) {
      	lsave << "[";
      	c1->Print( lsave.str().c_str() );//open the file
      	lsave.str("");//reset for adding the first plot
      	lsave << tmpstr ;
        if (shapes_[iShape].dology()) lsave << "_logy";
      	lsave << ".pdf" ;
      }

      std::cout << " -- iShape " << iShape << std::endl;
      std::cout << " -- Address of c1 "<< &c1 << std::endl;
      std::cout << " -- Address of lsave "<< lsave.str() << std::endl;
      c1->Print( lsave.str().c_str() );
      if (iShape==shapes_.size()-1) {
      	lsave << "]";
      	c1->Print(lsave.str().c_str());//close the file
      }

      lsave.str("");
      lsave << tmpstr << "_" << c1->GetName() ;
      if (shapes_[iShape].dology()) lsave << "_logy";
      lsave << ".pdf" ;
      c1->Print((lsave.str()).c_str());
//       lsavepng.str("");
//       lsavepng << tmpstr << "_" << c1->GetName();
//       if (shapes_[iShape].dology()) lsavepng << "_logy";
//       lsavepng << ".png" ;
//       c1->Print((lsavepng.str()).c_str());

//       lsavec.str("");
//       lsavec << tmpstr << "_" << c1->GetName();
//       if (shapes_[iShape].dology()) lsavec << "_logy";
//       lsavec << ".C" ;
//       c1->Print((lsavec.str()).c_str());

      //WRITE TO FILE
      writedir->cd();
      c1->Write();
      c1->Clear("D");
    }
    writedir->Close();

    return 0;
  };

}
