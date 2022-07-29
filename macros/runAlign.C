#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <chrono>
#include <fstream>
#include <iomanip>

#include <TFile.h>
#include <TTree.h>
#include <Rtypes.h>
#include <ROOT/RDataFrame.hxx>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"
#include "DetectorsCommonDataFormats/AlignParam.h"

#include "alignHelper.h"

#endif

// root -l
// .L ~/cernbox/alice/enigma/macros/runAlign.C++
// runAlign()

void runAlign(const Int_t fileStop = 1, // 4315,
              const bool preferAlignedFile = true)
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

  std::string generalPath = "/Users/andry/cernbox/alice/mft/pilotbeam/505713/";
  std::string alignStatus = "";
  if (preferAlignedFile || applyMisalignment) {
    alignStatus = "prealigned";
  } else {
    alignStatus = "idealgeo";
  }
  generalPath += alignStatus;

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

  // instantiate and configure the aligner

  AlignHelper aligner;

  aligner.setClusterDictionary(dict);

  AlignHelper::AlignConfig alignConfigParam;
  aligner.setChi2CutNStdDev(alignConfigParam.chi2CutNStdDev);
  aligner.setResidualCutInitial(alignConfigParam.residualCutInitial);
  aligner.setResidualCut(alignConfigParam.residualCut);
  aligner.setAllowedVariationDeltaX(alignConfigParam.allowedVarDeltaX);
  aligner.setAllowedVariationDeltaY(alignConfigParam.allowedVarDeltaY);
  aligner.setAllowedVariationDeltaZ(alignConfigParam.allowedVarDeltaZ);
  aligner.setAllowedVariationDeltaRz(alignConfigParam.allowedVarDeltaRz);
  aligner.setMinNumberClusterCut(alignConfigParam.minPoints);
  aligner.setChi2CutFactor(alignConfigParam.chi2CutFactor);

  // TODO: fix det. elements here

  // init Millipede
  aligner.init();

  // tree to record local measurements and residuals
  aligner.initTree();

  // compute Mille records

  int nRof = aligner.connectToTChains(mfttrackChain, mftclusterChain);
  for (int irof = 0; irof < nRof; irof++) { // loop on ROFs
    mftclusterChain->GetEntry(irof);
    mfttrackChain->GetEntry(irof);
    aligner.processROF();
    aligner.processRecoTracks();
  }
  aligner.printProcessTrackSummary();

  // compute alignment parameters

  // aligner.globalFit();

  // save alignment parameters to file

  std::vector<o2::detectors::AlignParam> alignParams;
  aligner.getAlignParams(alignParams);
  LOGF(info, "Storing MFT alignment params in local file %s/mft_alignment.root", generalPath.c_str());
  TFile afile(Form("%s/mft_alignment.root", generalPath.c_str()), "recreate", "", 505);
  afile.WriteObjectAny(&alignParams, "std::vector<o2::detectors::AlignParam>", "alignment");
  afile.Close();

  // save tree with to record local measurements and residuals and close related file
  aligner.closeTree();

  // the end

  std::chrono::steady_clock::time_point stop_time = std::chrono::steady_clock::now();

  std::cout << "----------------------------------- " << endl;
  std::cout << "Total Execution time: \t\t"
            << std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time).count()
            << " seconds" << endl;
}