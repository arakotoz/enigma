#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

#include <Rtypes.h>
#include <ROOT/RDataFrame.hxx>
#include <TChain.h>

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"

#include "MFTAlignment/TracksToRecords.h"

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
// .L ~/cernbox/alice/enigma/macros/runTracksToRecords.C++
// runTracksToRecords()

// fileStop = 44 for reco-with-mille/old-ctf/pass1
// fileStop = 4315 for prealigned/old-ctf
// fileStop = 173 for reco-with-mille/old-ctf/pass2
// fileStop = 819 for reco-with-mille/new-ctf/new-pass1
// fileStop = 819 for prealigned/new-ctf
void runTracksToRecords(const Int_t fileStop = 10,
                        const int minPoints = 6,
                        const bool preferAlignedFile = true,
                        const bool useMilleAlignment = true,
                        const bool useNewCTFs = true,
                        const bool doControl = true,
                        const int nEntriesAutoSave = 10000)
{

  ROOT::EnableImplicitMT(0);
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  // geometry

  ///< load geometry from file
  ///< When applyMisalignedment == false --> read from unaligned file
  ///< When preferAlignedFile == true and applyMisalignment == true : Prefer reading from existing aligned file

  const bool applyMisalignment = false;
  o2::base::GeometryManager::loadGeometry("", applyMisalignment, preferAlignedFile);
  o2::mft::GeometryTGeo* geom = o2::mft::GeometryTGeo::Instance();
  geom->fillMatrixCache(
    o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                             o2::math_utils::TransformType::L2G));

  // dictionary

  std::string dictFileName = "MFTdictionary.bin";
  if (useNewCTFs) {
    dictFileName = "o2_mft_dictionary.root";
  }
  o2::itsmft::TopologyDictionary* dict = nullptr;
  try {
    dict = o2::itsmft::TopologyDictionary::loadFrom(dictFileName);
  } catch (std::exception e) {
    std::cout << "Error " << e.what() << std::endl;
    return;
  }
  dict->readFromFile(dictFileName);

  // cluster and track chains

  const int runN = 505713;
  std::string basePath = "/Users/andry/cernbox/alice/mft/pilotbeam";
  std::string alignStatus = "";
  if (preferAlignedFile || applyMisalignment) {
    if (useMilleAlignment) {
      if (useNewCTFs) {
        alignStatus = "reco-with-mille/new-ctf/new-pass1";
      } else {
        alignStatus = "reco-with-mille/old-ctf/pass1";
      }
    } else {
      if (useNewCTFs) {
        alignStatus = "prealigned/new-ctf";
      } else {
        alignStatus = "prealigned/old-ctf";
      }
    }
  } else {
    if (useMilleAlignment) {
      alignStatus = "reco-with-mille/old-ctf/pass2";
    } else {
      alignStatus = "idealgeo/old-ctf";
    }
  }
  std::stringstream generalPathSs;
  generalPathSs << basePath << "/" << runN << "/" << alignStatus;
  std::string generalPath = generalPathSs.str();

  const Int_t fileStart = 1;
  TChain* mftclusterChain = new TChain("o2sim");
  TChain* mfttrackChain = new TChain("o2sim");

  /*
    // This is for reco-with-mille/old-ctf/pass2
    static constexpr int nFiles = 28;
    static constexpr std::array<int, nFiles> gridSubJob{
      4, 6, 7, 10, 11, 12, 13, 14, 15, 16,
      17, 18, 19, 20, 21, 22, 25, 26, 30, 32,
      34, 35, 37, 38, 39, 42, 43, 44};
    Int_t countFiles = 0;
    for (const auto& ii : gridSubJob) {
      if (countFiles > fileStop) {
        break;
      }
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
      countFiles++;
    }
    */

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
    LOG(debug) << "Add " << Form("%s/mftclusters.root", filePath.c_str());
    LOG(debug) << "Add " << Form("%s/mfttracks.root", filePath.c_str());
  }

  std::cout << "Number of files per chain = " << fileStop << std::endl;

  // instantiate and configure the aligner

  AlignConfigHelper alignConfigParam;
  alignConfigParam.minPoints = minPoints;

  o2::mft::TracksToRecords aligner;

  aligner.setRunNumber(runN);
  aligner.setBz(0.);

  aligner.setClusterDictionary(dict);
  aligner.setMinNumberClusterCut(alignConfigParam.minPoints);

  aligner.setWithControl(doControl);
  aligner.setNEntriesAutoSave(nEntriesAutoSave);

  // TODO: fix det. elements here

  // init Millipede

  aligner.init();

  // compute Mille records

  aligner.startRecordWriter();
  aligner.processROFs(mfttrackChain, mftclusterChain);
  aligner.printProcessTrackSummary();
  aligner.endRecordWriter();

  // the end

  std::chrono::steady_clock::time_point stop_time = std::chrono::steady_clock::now();

  std::cout << "----------------------------------- " << endl;
  std::cout << "Total Execution time: \t\t"
            << std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time).count()
            << " seconds" << endl;
}