#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <chrono>
#include <fstream>
#include <iomanip>
#include <string>

#include <Rtypes.h>
#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "Framework/Logger.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"
#include "MFTAlignment/AlignPointHelper.h"
#include "MFTAlignment/AlignPointControl.h"

#include "tracksToAlignControl.h"

#endif

// root -l
// .L ~/cernbox/alice/enigma/macros/runTracksToAlignInfoTree.C++
// runTracksToAlignInfoTree()

//__________________________________________________________________________
void runTracksToAlignInfoTree(const int fileStop = 4315,
                              const int minNumberClustersPerTrack = 4,
                              const bool preferAlignedFile = true,
                              const int nEntriesAutoSave = 10000)
{

  ROOT::EnableImplicitMT(0);
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  // geometry

  const bool applyMisalignment = false;
  o2::base::GeometryManager::loadGeometry("", applyMisalignment, preferAlignedFile);
  o2::mft::GeometryTGeo* geom = o2::mft::GeometryTGeo::Instance();
  geom->fillMatrixCache(
    o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                             o2::math_utils::TransformType::L2G));

  // dictionary

  std::string dictFileName = "MFTdictionary.bin";
  std::ifstream dictFile(dictFileName.c_str());
  if (!dictFile.good()) {
    std::cout << "Error: MFT dictionnary file " << dictFileName << " not found!" << std::endl;
    return;
  }
  o2::itsmft::TopologyDictionary* dict = new o2::itsmft::TopologyDictionary();
  dict->readBinaryFile(dictFileName);

  // cluster and track chains

  const int runN = 505713;
  std::string basePath = "/Users/andry/cernbox/alice/mft/pilotbeam";
  std::string alignStatus = "";
  if (preferAlignedFile || applyMisalignment) {
    alignStatus = "prealigned";
  } else {
    alignStatus = "idealgeo";
  }
  std::stringstream generalPathSs;
  generalPathSs << basePath << "/" << runN << "/" << alignStatus;
  std::string generalPath = generalPathSs.str();

  const Int_t fileStart = 1;
  TChain* mftclusterChain = new TChain("o2sim");
  TChain* mfttrackChain = new TChain("o2sim");
  for (Int_t ii = fileStart; ii <= fileStop; ii++) {
    std::stringstream ss;
    if (ii < 100) {
      ss << generalPath << "/"
         << std::setw(3) << std::setfill('0') << ii;
    } else {
      ss << generalPath << "/" << ii;
    }
    std::string filePath = ss.str();
    mftclusterChain->Add(Form("%s/mftclusters.root", filePath.c_str()));
    mfttrackChain->Add(Form("%s/mfttracks.root", filePath.c_str()));
  }
  std::cout << "Number of files per chain = " << fileStop << std::endl;

  // Alignment point info tree handler

  TracksToAlignControl* mftAlignController = new TracksToAlignControl();
  mftAlignController->setClusterDictionary(dict);
  mftAlignController->setMinNumberClusterCut(minNumberClustersPerTrack);
  mftAlignController->setNEntriesAutoSave(nEntriesAutoSave);
  mftAlignController->setOutFileName(Form("%s/mft_align_point.root", generalPath.c_str()));
  mftAlignController->init();

  // process cluster and track chains

  mftAlignController->processROFs(mfttrackChain, mftclusterChain);
  mftAlignController->printProcessTrackSummary();

  // the end

  std::chrono::steady_clock::time_point stop_time = std::chrono::steady_clock::now();

  std::cout << "----------------------------------- " << endl;
  std::cout << "Total Execution time: \t\t"
            << std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time).count()
            << " seconds" << endl;
}
