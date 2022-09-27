#include <vector>

#include <Rtypes.h>
#include <ROOT/RDataFrame.hxx>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "Framework/Logger.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "ReconstructionDataFormats/BaseCluster.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"
#include "MFTAlignment/AlignPointHelper.h"
#include "MFTAlignment/AlignPointControl.h"

class TracksToAlignControl
{
 public:
  TracksToAlignControl();
  ~TracksToAlignControl();

  void init();

  void setClusterDictionary(const o2::itsmft::TopologyDictionary* d) { mDictionary = d; }
  void setMinNumberClusterCut(const int value) { mMinNumberClusterCut = value; }
  void setNEntriesAutoSave(const int value) { mNEntriesAutoSave = value; }
  void setOutFileName(TString fname) { mPointControl.setOutFileName(fname); }

  void processROFs(TChain* mfttrackChain, TChain* mftclusterChain);
  void printProcessTrackSummary();

 protected:
  bool mIsInitDone;
  int mNumberOfClusterChainROFs;                           ///< number of ROFs in the cluster chain
  int mNumberOfTrackChainROFs;                             ///< number of ROFs in the track chain
  int mCounterSkippedTracks;                               ///< count how many tracks did not met the cut on the min. nb of clusters
  int mCounterUsedTracks;                                  ///< count how many tracks were used to make Mille records
  int mMinNumberClusterCut;                                ///< Minimum number of clusters in the track to be used for alignment
  const o2::itsmft::TopologyDictionary* mDictionary;       ///< cluster patterns dictionary
  o2::mft::AlignPointHelper* mAlignPoint;                  ///< Alignment point helper
  long mNEntriesAutoSave;                                  ///< number of entries needed to call AutoSave for the output TTrees
  o2::mft::AlignPointControl mPointControl;                ///< AlignPointControl handles the control tree
  std::vector<o2::BaseCluster<double>> mMFTClustersLocal;  ///< MFT clusters in local coordinate system
  std::vector<o2::BaseCluster<double>> mMFTClustersGlobal; ///< MFT clusters in global coordinate system

  ClassDefNV(TracksToAlignControl, 1);
};

//__________________________________________________________________________
TracksToAlignControl::TracksToAlignControl()
  : mIsInitDone(false),
    mNumberOfClusterChainROFs(0),
    mNumberOfTrackChainROFs(0),
    mCounterSkippedTracks(0),
    mCounterUsedTracks(0),
    mMinNumberClusterCut(6),
    mDictionary(nullptr),
    mAlignPoint(new o2::mft::AlignPointHelper()),
    mNEntriesAutoSave(10000)
{
  LOGF(debug, "TracksToAlignControl instantiated");
}

//__________________________________________________________________________
TracksToAlignControl::~TracksToAlignControl()
{
  if (mAlignPoint) {
    delete mAlignPoint;
  }
  if (mDictionary) {
    mDictionary = nullptr;
  }
  LOGF(debug, "TracksToAlignControl destroyed");
}

//__________________________________________________________________________
void TracksToAlignControl::init()
{
  if (mIsInitDone) {
    return;
  }
  if (mDictionary == nullptr) {
    LOGF(fatal, "TracksToAlignControl::init() failed because no cluster dictionary is defined");
    mIsInitDone = false;
    return;
  }
  mAlignPoint->setClusterDictionary(mDictionary);
  mPointControl.setCyclicAutoSave(mNEntriesAutoSave);

  mIsInitDone = true;
  LOGF(info, "TracksToAlignControl init done");
}

//__________________________________________________________________________
void TracksToAlignControl::processROFs(TChain* mfttrackChain, TChain* mftclusterChain)
{
  if (!mIsInitDone) {
    LOGF(fatal, "TracksToAlignControl::processROFs() aborted because init was not done !");
    return;
  }
  mPointControl.init();

  LOG(info) << "TracksToAlignControl::processROFs() - start";

  TTreeReader mftTrackChainReader(mfttrackChain);
  TTreeReader mftClusterChainReader(mftclusterChain);
  std::vector<unsigned char>::iterator pattIterator;

  TTreeReaderValue<std::vector<o2::mft::TrackMFT>> mftTracks =
    {mftTrackChainReader, "MFTTrack"};
  TTreeReaderValue<std::vector<o2::itsmft::ROFRecord>> mftTracksROF =
    {mftTrackChainReader, "MFTTracksROF"};
  TTreeReaderValue<std::vector<int>> mftTrackClusIdx =
    {mftTrackChainReader, "MFTTrackClusIdx"};

  TTreeReaderValue<std::vector<o2::itsmft::CompClusterExt>> mftClusters =
    {mftClusterChainReader, "MFTClusterComp"};
  TTreeReaderValue<std::vector<o2::itsmft::ROFRecord>> mftClustersROF =
    {mftClusterChainReader, "MFTClustersROF"};
  TTreeReaderValue<std::vector<unsigned char>> mftClusterPatterns =
    {mftClusterChainReader, "MFTClusterPatt"};

  int nCounterAllTracks = 0;

  while (mftTrackChainReader.Next() && mftClusterChainReader.Next()) {

    mNumberOfTrackChainROFs += (*mftTracksROF).size();
    mNumberOfClusterChainROFs += (*mftClustersROF).size();
    assert(mNumberOfTrackChainROFs == mNumberOfClusterChainROFs);

    pattIterator = (*mftClusterPatterns).begin();
    mAlignPoint->convertCompactClusters(
      *mftClusters, pattIterator, mMFTClustersLocal, mMFTClustersGlobal);

    //______________________________________________________
    for (auto& oneTrack : *mftTracks) { // track loop

      LOGF(debug, "Processing track # %5d", nCounterAllTracks);

      // Skip the track if not enough clusters
      auto ncls = oneTrack.getNumberOfPoints();
      if (ncls < mMinNumberClusterCut) {
        nCounterAllTracks++;
        mCounterSkippedTracks++;
        continue;
      }

      // Skip presumably quite low momentum track
      if (!oneTrack.isLTF()) {
        nCounterAllTracks++;
        mCounterSkippedTracks++;
        continue;
      }

      auto offset = oneTrack.getExternalClusterIndexOffset();

      // Store the initial track parameters
      mAlignPoint->resetTrackInitialParam();
      mAlignPoint->recordTrackInitialParam(oneTrack);

      bool isTrackUsed = true;

      for (int icls = 0; icls < ncls; ++icls) { // cluster loop

        mAlignPoint->resetAlignPoint();

        // Store measured positions
        auto clsEntry = (*mftTrackClusIdx)[offset + icls];
        auto localCluster = mMFTClustersLocal[clsEntry];
        auto globalCluster = mMFTClustersGlobal[clsEntry];
        mAlignPoint->setMeasuredPosition(localCluster, globalCluster);
        if (!mAlignPoint->isClusterOk()) {
          LOGF(warning, "TracksToAlignControl::processROFs() - will not use track # %5d with at least a bad cluster", nCounterAllTracks);
          mCounterSkippedTracks++;
          isTrackUsed = false;
          break;
        }

        // Propagate track to the current z plane of this cluster
        oneTrack.propagateParamToZlinear(mAlignPoint->getGlobalMeasuredPosition().Z());

        // Store reco positions
        mAlignPoint->setGlobalRecoPosition(oneTrack);

        // compute residuals
        mAlignPoint->setLocalResidual();
        mAlignPoint->setGlobalResidual();

        mPointControl.fill(mAlignPoint, mCounterUsedTracks);
        LOGF(debug, "TracksToAlignControl::processROFs() - track %i h %d d %d l %d s %4d lMpos x %.2e y %.2e z %.2e gMpos x %.2e y %.2e z %.2e gRpos x %.2e y %.2e z %.2e",
             mCounterUsedTracks, mAlignPoint->half(), mAlignPoint->disk(), mAlignPoint->layer(), mAlignPoint->getSensorId(),
             mAlignPoint->getLocalMeasuredPosition().X(), mAlignPoint->getLocalMeasuredPosition().Y(), mAlignPoint->getLocalMeasuredPosition().Z(),
             mAlignPoint->getGlobalMeasuredPosition().X(), mAlignPoint->getGlobalMeasuredPosition().Y(), mAlignPoint->getGlobalMeasuredPosition().Z(),
             mAlignPoint->getGlobalRecoPosition().X(), mAlignPoint->getGlobalRecoPosition().Y(), mAlignPoint->getGlobalRecoPosition().Z());

      } // end of loop on clusters

      if (isTrackUsed) {
        mCounterUsedTracks++;
      }
      nCounterAllTracks++;
    } // end of loop on tracks

  } // end of loop on TChain reader
  mPointControl.terminate();
  LOG(info) << "TracksToAlignControl::processROFs() - end";
}

//__________________________________________________________________________
void TracksToAlignControl::printProcessTrackSummary()
{
  LOGF(info, "TracksToAlignControl processRecoTracks() summary: ");
  LOGF(info,
       "n ROFs = %d, used tracks = %d, skipped tracks = %d",
       mNumberOfTrackChainROFs, mCounterUsedTracks, mCounterSkippedTracks);
}
