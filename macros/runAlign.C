#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <chrono>
#include <fstream>
#include <iomanip>
#include <string>

#include <Rtypes.h>
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>

#include "DetectorsCommonDataFormats/AlignParam.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
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
  TFile outAlignParamFile(Form("%s/mft_alignment.root", generalPath.c_str()), "recreate", "", 505);
  outAlignParamFile.WriteObjectAny(&alignParams, "std::vector<o2::detectors::AlignParam>", "alignment");
  outAlignParamFile.Close();

  // save Pede outputs to a file

  const int nDofPerSensor = aligner.getNDofPerSensor();
  std::vector<double> vPedeOutParams;
  std::vector<double> vPedeOutParamsErrors;
  std::vector<double> vPedeOutParamsPulls;
  aligner.getPedeOutParams(vPedeOutParams);
  aligner.getPedeOutParamsErrors(vPedeOutParamsErrors);
  aligner.getPedeOutParamsPulls(vPedeOutParamsPulls);
  TFile outPedefile(Form("%s/mft_pede_outputs.root", generalPath.c_str()), "recreate", "", 505);
  double dx = 0., dy = 0., dz = 0., dRz = 0.;
  double dxErr = 0., dyErr = 0., dzErr = 0., dRzErr = 0.;
  double dxPull = 0., dyPull = 0., dzPull = 0., dRzPull = 0.;
  TTree* tree = new TTree("pede", "the Pede output tree");
  tree->Branch("dx", &dx, "dx/D");
  tree->Branch("dy", &dy, "dy/D");
  tree->Branch("dz", &dz, "dz/D");
  tree->Branch("dRz", &dRz, "dRz/D");
  tree->Branch("dxErr", &dxErr, "dxErr/D");
  tree->Branch("dyErr", &dyErr, "dyErr/D");
  tree->Branch("dzErr", &dzErr, "dzErr/D");
  tree->Branch("dRzErr", &dRzErr, "dRzErr/D");
  tree->Branch("dxPull", &dxPull, "dxPull/D");
  tree->Branch("dyPull", &dyPull, "dyPull/D");
  tree->Branch("dzPull", &dzPull, "dzPull/D");
  tree->Branch("dRzPull", &dRzPull, "dRzPull/D");
  o2::itsmft::ChipMappingMFT chipMappingMFT;
  int NChips = o2::itsmft::ChipMappingMFT::NChips;
  for (int iChip = 0; iChip < NChips; iChip++) {
    dx = vPedeOutParams[iChip * nDofPerSensor + 0];
    dy = vPedeOutParams[iChip * nDofPerSensor + 1];
    dz = vPedeOutParams[iChip * nDofPerSensor + 3];
    dRz = vPedeOutParams[iChip * nDofPerSensor + 2];
    dxErr = vPedeOutParamsErrors[iChip * nDofPerSensor + 0];
    dyErr = vPedeOutParamsErrors[iChip * nDofPerSensor + 1];
    dzErr = vPedeOutParamsErrors[iChip * nDofPerSensor + 3];
    dRzErr = vPedeOutParamsErrors[iChip * nDofPerSensor + 2];
    dxPull = vPedeOutParamsPulls[iChip * nDofPerSensor + 0];
    dyPull = vPedeOutParamsPulls[iChip * nDofPerSensor + 1];
    dzPull = vPedeOutParamsPulls[iChip * nDofPerSensor + 3];
    dRzPull = vPedeOutParamsPulls[iChip * nDofPerSensor + 2];
    tree->Fill();
  }
  tree->Write();
  outPedefile.Close();

  // the end

  std::chrono::steady_clock::time_point stop_time = std::chrono::steady_clock::now();

  std::cout << "----------------------------------- " << endl;
  std::cout << "Total Execution time: \t\t"
            << std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time).count()
            << " seconds" << endl;
}