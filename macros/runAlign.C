#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <chrono>
#include <fstream>
#include <iomanip>
#include <string>

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

#include "MFTAlignment/Alignment.h"

#endif

struct AlignConfigHelper {
  int minPoints = 6;                ///< mininum number of clusters in a track used for alignment
  int chi2CutNStdDev = 3;           ///< Number of standard deviations for chi2 cut
  double residualCutInitial = 100.; ///< Cut on residual on first iteration
  double residualCut = 100.;        ///< Cut on residual for other iterations
  double allowedVarDeltaX = 0.5;    ///< allowed max delta in x-translation (cm)
  double allowedVarDeltaY = 0.5;    ///< allowed max delta in y-translation (cm)
  double allowedVarDeltaZ = 0.5;    ///< allowed max delta in z-translation (cm)
  double allowedVarDeltaRz = 0.01;  ///< allowed max delta in rotation around z-axis (rad)
  double chi2CutFactor = 256.;      ///< used to reject outliers i.e. bad tracks with sum(chi2) > Chi2DoFLim(fNStdDev, nDoF) * fChi2CutFactor
};

// alienv setenv O2Physics/latest -c root -l
// .L ~/cernbox/alice/enigma/macros/runAlign.C++
// runAlign()

void runAlign(const Int_t fileStop = 2,           // 4315,
              const double chi2CutFactor = 65536, // 256
              const bool doWriteRecords = true,
              const bool doGlobalFit = true,
              const bool preferAlignedFile = true,
              const bool doControl = true,
              const int nEntriesAutoSave = 500)
{
  if (!doWriteRecords && !doGlobalFit)
    return;

  ROOT::EnableImplicitMT(0);
  ROOT::EnableThreadSafety();
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  AlignConfigHelper alignConfigParam;
  alignConfigParam.chi2CutFactor = chi2CutFactor;

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

  TChain* mftclusterChain = nullptr;
  TChain* mfttrackChain = nullptr;
  std::string generalPath = "/Users/andry/cernbox/alice/mft/pilotbeam/505713/";
  std::string alignStatus = "";
  if (preferAlignedFile || applyMisalignment) {
    alignStatus = "prealigned";
  } else {
    alignStatus = "idealgeo";
  }
  generalPath += alignStatus;

  if (doWriteRecords) {
    const Int_t fileStart = 1;
    mftclusterChain = new TChain("o2sim");
    mfttrackChain = new TChain("o2sim");
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
  }

  // instantiate and configure the aligner

  o2::mft::Alignment aligner;

  aligner.setClusterDictionary(dict);

  aligner.setChi2CutNStdDev(alignConfigParam.chi2CutNStdDev);
  aligner.setResidualCutInitial(alignConfigParam.residualCutInitial);
  aligner.setResidualCut(alignConfigParam.residualCut);
  aligner.setAllowedVariationDeltaX(alignConfigParam.allowedVarDeltaX);
  aligner.setAllowedVariationDeltaY(alignConfigParam.allowedVarDeltaY);
  aligner.setAllowedVariationDeltaZ(alignConfigParam.allowedVarDeltaZ);
  aligner.setAllowedVariationDeltaRz(alignConfigParam.allowedVarDeltaRz);
  aligner.setMinNumberClusterCut(alignConfigParam.minPoints);
  aligner.setChi2CutFactor(alignConfigParam.chi2CutFactor);

  aligner.setWithControl(doControl);
  aligner.setNEntriesAutoSave(nEntriesAutoSave);
  aligner.setWithRecordWriter(doWriteRecords);
  aligner.setWithRecordReader(doGlobalFit);

  // TODO: fix det. elements here

  // init Millipede
  aligner.init();

  if (doWriteRecords) {
    // compute Mille records
    aligner.startRecordWriter();
    aligner.processROFs(mfttrackChain, mftclusterChain);
    aligner.printProcessTrackSummary();
    aligner.endRecordWriter();
  }

  TChain* recordChain = nullptr;

  if (doGlobalFit) {
    // data records
    recordChain = new TChain("milleRecords");
    mftclusterChain->Add("mft_mille_records.root");

    // compute alignment parameters
    aligner.connectRecordReaderToChain(recordChain);
    aligner.globalFit();

    // save alignment parameters to file
    std::vector<o2::detectors::AlignParam> alignParams;
    aligner.getAlignParams(alignParams);
    LOGF(info, "Storing MFT alignment params in local file %s/mft_alignment.root", generalPath.c_str());
    TFile afile(Form("%s/mft_alignment.root", generalPath.c_str()), "recreate", "", 505);
    afile.WriteObjectAny(&alignParams, "std::vector<o2::detectors::AlignParam>", "alignment");
    afile.Close();
  }

  // the end

  std::chrono::steady_clock::time_point stop_time = std::chrono::steady_clock::now();

  std::cout << "----------------------------------- " << endl;
  std::cout << "Total Execution time: \t\t"
            << std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time).count()
            << " seconds" << endl;
}