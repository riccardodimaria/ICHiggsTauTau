#include <iostream>
#include <vector>
#include <string>
#include "TSystem.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "Core/interface/AnalysisBase.h"
#include "Utilities/interface/JsonTools.h"
#include "Utilities/interface/FnRootTools.h"
#include "HiggsTauTau/interface/HTTConfig.h"
#include "HiggsTauTau/interface/QCDStudy.h"

using std::string;
using std::vector;

using ic::AnalysisBase;
using std::vector;
using std::string;

int main(int argc, char* argv[]) {

  Json::Value const js = ic::MergedJson(argc, argv);

  fwlite::TFileService* fs =
      new fwlite::TFileService((js["output"].asString()));

  vector<string> files = ic::ParseFileLines(js["filelist"].asString());
  for (auto & f : files) f = js.get("file_prefix", std::string("")).asString() + f;

  AnalysisBase analysis("QCDStudy", files, "icEventProducer/EventTree",
                        js["max_events"].asInt64());
  analysis.SetTTreeCaching(true);
  analysis.StopOnFileFailure(true);
  analysis.RetryFileAfterFailure(7, 3);
  analysis.CalculateTimings(js.get("timings", false).asBool());


  auto qcd_module = ic::QCDStudy("QCDStudy").set_fs(fs)
                      .set_channel(ic::channel::mt);

  analysis.AddModule(&qcd_module);

  analysis.RunAnalysis();

  delete fs;

  return 0;
}
