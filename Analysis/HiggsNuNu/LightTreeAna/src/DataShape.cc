#include "UserCode/ICHiggsTauTau/Analysis/HiggsNuNu/LightTreeAna/interface/DataShape.h"
#include <iostream>
#include "TH1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include <map>

namespace ic{

  DataShape::DataShape(std::string name) : LTModule(name){
    dataweight_="1";//"weight_nolepnotrig*weight_trig_0";
    std::vector<std::string> shapes;
    shapes.push_back("jet2_pt(200,0.,1000.)");
    shape_=shapes;
    dirname_="";
  };

  DataShape::~DataShape(){ ;};

  int DataShape::Init(TFile* fs){
    fs_=fs;
    std::cout<<"Initialisation info for "<<module_name_<<":"<<std::endl;
    std::cout<<"  Data sets are: "<<std::endl;
      for(unsigned iset=0;iset<dataset_.size();iset++){
	std::cout<<"    "<<dataset_[iset]<<std::endl;
      }
    std::cout<<"  Base selection is: "<<basesel_<<std::endl;
    std::cout<<"  Extra selection is: "<<cat_<<std::endl;
    std::cout<<"  Data weight is: "<<dataweight_<<std::endl;
    if((shapename_.size()!=shape_.size())&&shapename_.size()!=0){
      std::cout<<"  WARNING: different numbers of shape names " << shapename_.size() << " and shapes" << shape_.size() << "! Exiting..."<<std::endl;
      return 1;
    }
    if (!fs_) {
      std::cout << " Error in getting pointer to output file in DataShape, exiting! " << std::endl;
      return 1;
    }
    return 0;
  };

  int DataShape::Run(LTFiles* filemanager){
    std::cout<<"DataShape " << module_name_<<":"<<std::endl;

    TFile *file=fs_;
    TDirectory* dir;
    if(dirname_==""){
      dir=file->mkdir(this->module_name().c_str());
    }
    else if(!fs_->GetDirectory(dirname_.c_str())){
      dir=file->mkdir(dirname_.c_str());
    }
    else{
      dir=file->GetDirectory(dirname_.c_str());
    }
    dir->cd();
    //Get Shapes for Data
    std::cout<<"  Getting shape"<<std::endl;

    for(unsigned iShape=0;iShape<shape_.size();iShape++){
      std::string histname;
      histname=shape_[iShape];
      //access original special characters to recognise 2D hists
      size_t n = std::count(histname.begin(), histname.end(), ':');
      if (histname.find("::")!=histname.npos) n-=2;
      //replace special characters for root to store hist properly
      if(shapename_.size()==0){
	histname=extractShapeName(shape_[iShape]);
      }
      else{
	histname=shapename_[iShape];
      }
      TH1F datashape;
      TH2F datashape2D;
      TH3F datashape3D;
      bool is2D = false;
      bool is3D = false;
      std::cout << " Identifying histo's dimension: " << histname << ": " << n+1 << std::endl;
      if (n==1) {
	datashape2D=filemanager->GetSetsShape2D(dataset_,shape_[iShape],basesel_,cat_,dataweight_,false);
	is2D=true;
      } else if (n==2){
	datashape3D=filemanager->GetSetsShape3D(dataset_,shape_[iShape],basesel_,cat_,dataweight_,false);
	is3D=true;
      }
      else if (n==0) datashape = filemanager->GetSetsShape(dataset_,shape_[iShape],basesel_,cat_,dataweight_,false);
      else {
	std::cout << "Warning, more than 2 : found in histname " << histname << ": " << n << " file " << __FILE__ << " line " << __LINE__ << std::endl;
	return 1;
      }
      dir->cd();
      if(iShape==0){
	if (!is2D && !is3D) std::cout<<"  nevents: "<<Integral(&datashape)<<"+-"<<Error(&datashape)<<std::endl;
	else if (is2D) std::cout<<"  nevents: "<<Integral(&datashape2D)<<"+-"<<Error(&datashape2D)<<std::endl;
	else if (is3D) std::cout<<"  nevents: "<<Integral(&datashape3D)<<"+-"<<Error(&datashape3D)<<std::endl;
      }
      if (is2D){
	datashape2D.SetName(histname.c_str());
	datashape2D.Write();
      } else if (is3D){
	datashape3D.SetName(histname.c_str());
	datashape3D.Write();
      } else {
	datashape.SetName(histname.c_str());
	datashape.Write();
      }
    }
    return 0;
  };

}
