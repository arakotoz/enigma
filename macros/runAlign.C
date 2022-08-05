#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <chrono>
#include <fstream>
#include <iomanip>
#include <string>

#include <Rtypes.h>
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TChain.h>

#include "DetectorsCommonDataFormats/AlignParam.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"

#include "MFTAlignment/RecordsToAlignParams.h"

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

void runAlign(const double chi2CutFactor = 65536, // 256
              const bool preferAlignedFile = true,
              const bool doControl = true,
              const int nEntriesAutoSave = 500)
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

  // instantiate and configure the aligner

  AlignConfigHelper alignConfigParam;
  alignConfigParam.chi2CutFactor = chi2CutFactor;

  o2::mft::RecordsToAlignParams aligner;

  aligner.setChi2CutNStdDev(alignConfigParam.chi2CutNStdDev);
  aligner.setResidualCutInitial(alignConfigParam.residualCutInitial);
  aligner.setResidualCut(alignConfigParam.residualCut);
  aligner.setAllowedVariationDeltaX(alignConfigParam.allowedVarDeltaX);
  aligner.setAllowedVariationDeltaY(alignConfigParam.allowedVarDeltaY);
  aligner.setAllowedVariationDeltaZ(alignConfigParam.allowedVarDeltaZ);
  aligner.setAllowedVariationDeltaRz(alignConfigParam.allowedVarDeltaRz);
  aligner.setChi2CutFactor(alignConfigParam.chi2CutFactor);

  aligner.setWithControl(doControl);
  aligner.setNEntriesAutoSave(nEntriesAutoSave);

  // TODO: fix det. elements here

  // init Millipede

  aligner.init();

  // records chains

  std::string generalPath = "/Users/andry/cernbox/alice/mft/pilotbeam/505713/";
  std::string alignStatus = "";
  if (preferAlignedFile || applyMisalignment) {
    alignStatus = "prealigned";
  } else {
    alignStatus = "idealgeo";
  }
  generalPath += alignStatus;

  TChain* recordChain = new TChain("milleRecords");
  recordChain->Add("mft_mille_records.root");

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

  // the end

  std::chrono::steady_clock::time_point stop_time = std::chrono::steady_clock::now();

  std::cout << "----------------------------------- " << endl;
  std::cout << "Total Execution time: \t\t"
            << std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time).count()
            << " seconds" << endl;
}