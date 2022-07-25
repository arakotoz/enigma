#include <array>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <Rtypes.h>
#include <TString.h>
#include <TChain.h>

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "ReconstructionDataFormats/BaseCluster.h"
#include "MFTAlignment/MillePedeRecord.h"
#include "MFTAlignment/MillePede2.h"
#include "MFTAlignment/AlignPointHelper.h"
#include "MFTBase/GeometryTGeo.h"

#include "Framework/Logger.h"
#include "MFTAlignment/AlignSensorHelper.h"
#include "MFTTracking/IOUtils.h"
#include "MFTBase/Geometry.h"
#include "MFTAlignment/MillePede2.h"

using namespace o2::mft;

class AlignHelper
{
 public:
  struct AlignConfig {
    int minPoints = 6;                  ///< mininum number of clusters in a track used for alignment
    Int_t chi2CutNStdDev = 3;           ///< Number of standard deviations for chi2 cut
    Double_t residualCutInitial = 100.; ///< Cut on residual on first iteration
    Double_t residualCut = 100.;        ///< Cut on residual for other iterations
    Double_t allowedVarDeltaX = 0.5;    ///< allowed max delta in x-translation (cm)
    Double_t allowedVarDeltaY = 0.5;    ///< allowed max delta in y-translation (cm)
    Double_t allowedVarDeltaZ = 0.5;    ///< allowed max delta in z-translation (cm)
    Double_t allowedVarDeltaRz = 0.01;  ///< allowed max delta in rotation around z-axis (rad)
  };

 public:
  /// \brief construtor
  AlignHelper();

  /// \brief destructor
  ~AlignHelper() = default;

  /// \brief init Millipede and AlignPointHelper
  void init();

  // simple setters

  void setClusterDictionary(const o2::itsmft::TopologyDictionary* d) { mDictionary = d; }
  void setRunNumber(const int value) { mRunNumber = value; }
  void setBz(const float bz) { mBz = bz; }
  void setSaveTrackRecordToFile(const bool choice) { mSaveTrackRecordToFile = choice; }
  void setChi2CutNStdDev(const Int_t value) { mChi2CutNStdDev = value; }
  void setResidualCutInitial(const Double_t value) { mResCutInitial = value; }
  void setResidualCut(const Double_t value) { mResCut = value; }
  void setMinNumberClusterCut(const int value) { mMinNumberClusterCut = value; }
  void setAllowedVariationDeltaX(const double value) { mAllowVar[0] = value; }
  void setAllowedVariationDeltaY(const double value) { mAllowVar[1] = value; }
  void setAllowedVariationDeltaZ(const double value) { mAllowVar[3] = value; }
  void setAllowedVariationDeltaRz(const double value) { mAllowVar[2] = value; }

  /// \brief set pointer to geometry prepared outside of the class i.e. already had fillMatrixCache()
  void setGeometry(const o2::mft::GeometryTGeo* geom) { mGeometry = geom; }

  /// \brief access mft tracks and clusters in the ROOT files, return total number of ROFs
  int connectToTChains(TChain* mfttrackChain, TChain* mftclusterChain);

  /// \brief process a given ROF
  void processROF();

  /// \brief use valid tracks to build Mille records
  void processRecoTracks();

  /// \brief perform the simultaneous fit of track and alignement parameters
  void globalFit();

  /// \brief print a summary status of what happened in processRecoTracks()
  void printProcessTrackSummary();

  /// \brief provide access to the AlignParam vector
  void getAlignParams(std::vector<o2::detectors::AlignParam>& alignParams) { alignParams = mAlignParams; }

 protected:
  int mRunNumber = 0;                                                            ///< run number
  float mBz = 0;                                                                 ///< magnetic field status
  int mNumberOfClusterChainROFs = 0;                                             ///< number of ROFs in the cluster chain
  int mNumberOfTrackChainROFs = 0;                                               ///< number of ROFs in the track chain
  int mCounterLocalEquationFailed = 0;                                           ///< count how many times we failed to set a local equation
  int mCounterSkippedTracks = 0;                                                 ///< count how many tracks did not met the cut on the min. nb of clusters
  int mCounterUsedTracks = 0;                                                    ///< count how many tracks were used to make Mille records
  static constexpr int mNumberOfTrackParam = 4;                                  ///< Number of track (= local) parameters (X0, Tx, Y0, Ty)
  static constexpr int mNDofPerSensor = 4;                                       ///< translation in global x, y, z, and rotation Rz around global z-axis
  static o2::itsmft::ChipMappingMFT mChipMapping;                                ///< MFT chip <-> ladder, layer, disk, half mapping
  static constexpr int mNumberOfSensors = mChipMapping.getNChips();              ///< Total number of sensors (detection elements) in the MFT
  static constexpr int mNumberOfGlobalParam = mNDofPerSensor * mNumberOfSensors; ///< Number of alignment (= global) parameters
  Double_t mGlobalDerivatives[mNumberOfGlobalParam];                             ///< Array of global derivatives {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}
  Double_t mLocalDerivatives[mNumberOfTrackParam];                               ///< Array of local derivatives {dX0, dTx, dY0, dTz}
  std::array<Double_t, mNDofPerSensor> mAllowVar;                                ///< "Encouraged" variation for degrees of freedom {dx, dy, dRz, dz}
  double mStartFac = 256;                                                        ///< Initial value for chi2 cut (if > 1, iterations in Millepede are turned on)
  Int_t mChi2CutNStdDev = 3;                                                     ///< Number of standard deviations for chi2 cut
  Double_t mResCutInitial = 100.;                                                ///< Cut on residual on first iteration
  Double_t mResCut = 100.;                                                       ///< Cut on residual for other iterations
  int mMinNumberClusterCut = 6;                                                  ///< Minimum number of clusters in the track to be used for alignment
  o2::mft::MillePedeRecord mTrackRecord;                                         ///< running MillePede Track record
  double mWeightRecord = 1.;                                                     ///< the weight given to a single Mille record in Millepede algorithm
  bool mSaveTrackRecordToFile = false;                                           ///< true if we want to save Mille records in a ROOT file
  TString mMilleRecordsFileName;                                                 ///< output file name when saving the Mille records
  TString mMilleConstraintsRecFileName;                                          ///< output file name when saving the records of the constraints
  std::unique_ptr<o2::mft::MillePede2> mMillepede = nullptr;                     ///< Millepede2 implementation copied from AliROOT
  const o2::itsmft::TopologyDictionary* mDictionary = nullptr;                   ///< cluster patterns dictionary
  std::unique_ptr<o2::mft::AlignPointHelper> mAlignPoint = nullptr;              ///< AlignHelper point helper
  std::vector<o2::detectors::AlignParam> mAlignParams;                           ///< vector of alignment parameters computed by Millepede global fit
  const o2::mft::GeometryTGeo* mGeometry = nullptr;                              ///< geometry that must be initialised outside of Alignment
  bool mIsInitDone = false;                                                      ///< boolean to follow the initialisation status
  Int_t mGlobalParameterStatus[mNumberOfGlobalParam];                            ///< Array of effective degrees of freedom, used to fix detectors, parameters, etc.

  // used to fix some degrees of freedom

  static constexpr Int_t mFixedParId = -1;
  static constexpr Int_t mFreeParId = mFixedParId - 1;

  // access these data from mfttracks.root and mftclusters.root

  std::vector<o2::mft::TrackMFT> mMFTTracks, *mMFTTracksVecP = &mMFTTracks;
  std::vector<o2::itsmft::ROFRecord> mMFTTracksROF, *mMFTTracksROFVecP = &mMFTTracksROF;
  std::vector<int> mMFTTrackClusIdx, *mMFTTrackClusIdxVecP = &mMFTTrackClusIdx;
  std::vector<o2::itsmft::CompClusterExt> mMFTClusters, *mMFTClustersVecP = &mMFTClusters;
  std::vector<o2::itsmft::ROFRecord> mMFTClustersROF, *mMFTClustersROFVecP = &mMFTClustersROF;
  std::vector<unsigned char> mMFTClusterPatterns, *mMFTClusterPatternsVecP = &mMFTClusterPatterns;
  std::vector<o2::BaseCluster<float>> mMFTClustersGlobal;

  /// \brief set array of local derivatives
  bool setLocalDerivative(Int_t index, Double_t value);

  /// \brief set array of global derivatives
  bool setGlobalDerivative(Int_t index, Double_t value);

  /// \brief reset the array of the Local derivative
  void resetLocalDerivative();

  /// \brief reset the array of the Global derivative
  void resetGlocalDerivative();

  /// \brief set the first component of the local equation vector for a given alignment point
  bool setLocalEquationX();

  /// \brief set the 2nd component of the local equation vector for a given alignment point
  bool setLocalEquationY();

  /// \brief set the last component of the local equation vector for a given alignment point
  bool setLocalEquationZ();

  ClassDefNV(AlignHelper, 1);
};

//__________________________________________________________________________
AlignHelper::AlignHelper()
  : mRunNumber(0),
    mBz(0),
    mNumberOfClusterChainROFs(0),
    mNumberOfTrackChainROFs(0),
    mCounterLocalEquationFailed(0),
    mCounterSkippedTracks(0),
    mCounterUsedTracks(0),
    mStartFac(256),
    mChi2CutNStdDev(3),
    mResCutInitial(100.),
    mResCut(100.),
    mMinNumberClusterCut(6),
    mWeightRecord(1.),
    mSaveTrackRecordToFile(false),
    mMilleRecordsFileName("mft_mille_records.root"),
    mMilleConstraintsRecFileName("mft_mille_constraints.root"),
    mIsInitDone(false)
{
  mMillepede = std::make_unique<MillePede2>();
  // default allowed variations w.r.t. global system coordinates
  mAllowVar[0] = 0.5;  // delta translation in x (cm)
  mAllowVar[1] = 0.5;  // delta translation in y (cm)
  mAllowVar[2] = 0.01; // rotation angle Rz around z-axis (rad)
  mAllowVar[3] = 0.5;  // delta translation in z (cm)

  for (int iPar = 0; iPar < mNumberOfGlobalParam; iPar++) {
    mGlobalParameterStatus[iPar] = mFreeParId;
  }
}

//__________________________________________________________________________
void AlignHelper::init()
{
  if (mIsInitDone)
    return;
  if (mGeometry == nullptr) {
    LOGF(fatal, "AlignHelper::init() failed because no geometry is defined");
    mIsInitDone = false;
    return;
  }
  mAlignPoint = std::make_unique<AlignPointHelper>(mGeometry);
  mMillepede->InitMille(mNumberOfGlobalParam,
                        mNumberOfTrackParam,
                        mChi2CutNStdDev,
                        mResCut,
                        mResCutInitial);
  // filenames for the processed tracks and constraints records
  mMillepede->SetDataRecFName(mMilleRecordsFileName.Data());
  mMillepede->SetConsRecFName(mMilleConstraintsRecFileName.Data());

  if (mSaveTrackRecordToFile) {
    mMillepede->InitDataRecStorage(kFALSE);
  }

  LOG(info) << "-------------- Alignment configured with -----------------";
  LOGF(info, "Chi2CutNStdDev = %d", mChi2CutNStdDev);
  LOGF(info, "ResidualCutInitial = %.3f", mResCutInitial);
  LOGF(info, "ResidualCut = %.3f", mResCut);
  LOGF(info, "MinNumberClusterCut = %d", mMinNumberClusterCut);
  LOGF(info, "mStartFac = %.3f", mStartFac);
  LOGF(info,
       "Allowed variation: dx = %.3f, dy = %.3f, dz = %.3f, dRz = %.4f",
       mAllowVar[0], mAllowVar[1], mAllowVar[3], mAllowVar[2]);
  LOG(info) << "-----------------------------------------------------------";

  // set allowed variations for all parameters
  for (int chipId = 0; chipId < mNumberOfSensors; ++chipId) {
    for (Int_t iPar = 0; iPar < mNDofPerSensor; ++iPar) {
      mMillepede->SetParSigma(chipId * mNDofPerSensor + iPar, mAllowVar[iPar]);
    }
  }

  // set iterations
  if (mStartFac > 1) {
    mMillepede->SetIterations(mStartFac);
  }

  LOGF(info, "AlignHelper init done");
  mIsInitDone = true;
}

//__________________________________________________________________________
int AlignHelper::connectToTChains(TChain* mfttrackChain, TChain* mftclusterChain)
{
  // get tracks
  mfttrackChain->SetBranchAddress("MFTTrack", &mMFTTracksVecP);
  mfttrackChain->SetBranchAddress("MFTTracksROF", &mMFTTracksROFVecP);
  mfttrackChain->SetBranchAddress("MFTTrackClusIdx", &mMFTTrackClusIdxVecP);
  mNumberOfTrackChainROFs = mfttrackChain->GetEntries();
  LOG(info) << "Number of track chain entries = " << mNumberOfTrackChainROFs;

  // get clusters
  mftclusterChain->SetBranchAddress("MFTClusterComp", &mMFTClustersVecP);
  mftclusterChain->SetBranchAddress("MFTClustersROF", &mMFTClustersROFVecP);
  mftclusterChain->SetBranchAddress("MFTClusterPatt", &mMFTClusterPatternsVecP);
  mNumberOfClusterChainROFs = mftclusterChain->GetEntries();
  LOG(info) << "Number of cluster chain entries = " << mNumberOfClusterChainROFs;

  assert(mNumberOfTrackChainROFs == mNumberOfClusterChainROFs);
  return mNumberOfTrackChainROFs;
}

//__________________________________________________________________________
void AlignHelper::processROF()
{
  std::vector<unsigned char>::iterator pattIt = mMFTClusterPatterns.begin();
  mMFTClustersGlobal.clear();
  mMFTClustersGlobal.reserve(mMFTClusters.size());

  for (const auto& compCluster : mMFTClusters) {

    auto chipID = compCluster.getChipID();
    auto pattID = compCluster.getPatternID();

    o2::math_utils::Point3D<float> locXYZ;
    // Dummy COG errors (about half pixel size)
    float sigmaX2 = o2::mft::ioutils::DefClusError2Row;
    float sigmaY2 = o2::mft::ioutils::DefClusError2Col;
    if (pattID != o2::itsmft::CompCluster::InvalidPatternID) {
      // ALPIDE local Y coordinate => MFT global X coordinate (ALPIDE rows)
      sigmaX2 = mDictionary->getErr2X(pattID);
      // ALPIDE local Z coordinate => MFT global Y coordinate (ALPIDE columns)
      sigmaY2 = mDictionary->getErr2Z(pattID);
      if (!mDictionary->isGroup(pattID)) {
        locXYZ = mDictionary->getClusterCoordinates(compCluster);
      } else {
        o2::itsmft::ClusterPattern patt(pattIt);
        locXYZ = mDictionary->getClusterCoordinates(compCluster, patt);
      }
    } else {
      o2::itsmft::ClusterPattern patt(pattIt);
      locXYZ = mDictionary->getClusterCoordinates(compCluster, patt, false);
    }
    mMFTClustersGlobal.emplace_back(chipID, locXYZ);
  }
}

//__________________________________________________________________________
void AlignHelper::processRecoTracks()
{
  if (!mIsInitDone) {
    LOGF(fatal, "AlignHelper::processRecoTracks() aborted because init was not done!");
    return;
  }

  for (auto oneTrack : mMFTTracks) { // track loop

    // Skip the track if not enough clusters
    auto ncls = oneTrack.getNumberOfPoints();
    if (ncls < mMinNumberClusterCut) {
      mCounterSkippedTracks++;
      continue;
    }

    auto offset = oneTrack.getExternalClusterIndexOffset();

    mTrackRecord.Reset();
    if (mMillepede->GetRecord()) {
      mMillepede->GetRecord()->Reset();
    }

    // Store the initial track parameters
    mAlignPoint->resetTrackInitialParam();
    mAlignPoint->recordTrackInitialParam(oneTrack);

    bool isTrackUsed = true;

    for (int icls = 0; icls < ncls; ++icls) { // cluster loop

      mAlignPoint->resetDerivatives();
      mAlignPoint->resetAlignPoint();

      auto clsEntry = mMFTTrackClusIdx[offset + icls];
      auto globalCluster = mMFTClustersGlobal[clsEntry];

      // Propagate track to the current z plane of this cluster
      oneTrack.propagateParamToZlinear(globalCluster.getZ());

      // Store reco and measured positions
      mAlignPoint->setGlobalRecoPosition(oneTrack);
      mAlignPoint->setLocalMeasuredPosition(globalCluster);

      // Compute derivatives
      mAlignPoint->computeLocalDerivatives();
      mAlignPoint->computeGlobalDerivatives();

      // Set local equations
      bool success = true;
      success &= setLocalEquationX();
      success &= setLocalEquationY();
      success &= setLocalEquationZ();
      isTrackUsed &= success;

    } // end of loop on clusters

    if (isTrackUsed) {
      // copy track record
      mMillepede->SetRecordRun(mRunNumber);
      mMillepede->SetRecordWeight(mWeightRecord);
      mTrackRecord = *mMillepede->GetRecord();

      // save record data
      if (mSaveTrackRecordToFile) {
        mMillepede->SaveRecordData();
      }

      mCounterUsedTracks++;
    }
  } // end of loop on tracks

  if (mSaveTrackRecordToFile) {
    mMillepede->CloseDataRecStorage();
  }
}

//__________________________________________________________________________
void AlignHelper::globalFit()
{
  if (!mIsInitDone) {
    LOGF(fatal, "AlignHelper::globalFit() aborted because init was not done!");
    return;
  }

  if (!mCounterUsedTracks) {
    LOGF(fatal, "AlignHelper::globalFit() aborted because no reco track was used!");
    return;
  }

  // allocate memory in arrays to temporarily store the results of the global fit

  Double_t* params = new Double_t[mNumberOfGlobalParam];
  Double_t* paramsErrors = new Double_t[mNumberOfGlobalParam];
  Double_t* paramsPulls = new Double_t[mNumberOfGlobalParam];

  // initialise the content of each array

  for (int ii = 0; ii < mNumberOfGlobalParam; ii++) {
    params[ii] = 0.;
    paramsErrors[ii] = 0.;
    paramsPulls[ii] = 0.;
  }

  // perform the simultaneous fit of track and alignement parameters

  mMillepede->GlobalFit(params, paramsErrors, paramsPulls);

  // post-treatment:
  // debug output + save Millepede global fit result in AlignParam vector

  LOGF(info, "AlignHelper: done fitting global parameters");
  LOGF(info, "sensor info, dx (cm), dy (cm), dz (cm), dRz (rad)");

  AlignSensorHelper chipHelper(mGeometry);
  double dRx = 0., dRy = 0., dRz = 0.; // delta rotations
  double dx = 0., dy = 0., dz = 0.;    // delta translations
  bool global = true;                  // delta in global ref. system
  bool withSymName = false;

  for (int chipId = 0; chipId < mNumberOfSensors; chipId++) {
    chipHelper.setSensorOnlyInfo(chipId);
    std::stringstream name = chipHelper.getSensorFullName(withSymName);
    dx = params[chipId * mNDofPerSensor + 0];
    dy = params[chipId * mNDofPerSensor + 1];
    dz = params[chipId * mNDofPerSensor + 3];
    dRz = params[chipId * mNDofPerSensor + 2];
    LOGF(info,
         "%s, %.3e, %.3e, %.3e, %.3e",
         name.str().c_str(), dx, dy, dz, dRz);
    mAlignParams.emplace_back(
      chipHelper.geoSymbolicName(),
      chipHelper.sensorUid(),
      dx, dy, dz,
      dRx, dRy, dRz,
      global);
  }

  // free allocated memory

  delete[] params;
  delete[] paramsErrors;
  delete[] paramsPulls;
}

//__________________________________________________________________________
void AlignHelper::printProcessTrackSummary()
{
  LOGF(info, "AlignHelper processRecoTracks() summary: ");
  LOGF(info,
       "n ROFs = %d, used tracks = %d, skipped tracks = %d, local equations failed = %d",
       mNumberOfTrackChainROFs, mCounterUsedTracks,
       mCounterSkippedTracks, mCounterLocalEquationFailed);
}

//__________________________________________________________________________
bool AlignHelper::setLocalDerivative(Int_t index, Double_t value)
{
  // index [0 .. 3] for {dX0, dTx, dY0, dTz}

  bool success = true;
  if (index < mNumberOfTrackParam) {
    mLocalDerivatives[index] = value;
  } else {
    LOGF(error,
         "AlignHelper::setLocalDerivative() - index %d >= %d",
         index, mNumberOfTrackParam);
    success = false;
  }
  return success;
}

//__________________________________________________________________________
bool AlignHelper::setGlobalDerivative(Int_t index, Double_t value)
{
  // index [0 .. 3] for {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}

  bool success = true;
  if (index < mNumberOfGlobalParam) {
    mGlobalDerivatives[index] = value;
  } else {
    LOGF(error,
         "AlignHelper::setGlobalDerivative() - index %d >= %d",
         index, mNumberOfGlobalParam);
    success = false;
  }
  return success;
}

//__________________________________________________________________________
void AlignHelper::resetLocalDerivative()
{
  for (int i = 0; i < mNumberOfTrackParam; ++i) {
    mLocalDerivatives[i] = 0.0;
  }
}

//__________________________________________________________________________
void AlignHelper::resetGlocalDerivative()
{
  for (int i = 0; i < mNumberOfGlobalParam; ++i) {
    mGlobalDerivatives[i] = 0.0;
  }
}

//__________________________________________________________________________
bool AlignHelper::setLocalEquationX()
{

  if (!mAlignPoint->isAlignPointSet())
    return false;
  if (!mAlignPoint->isGlobalDerivativeDone())
    return false;
  if (!mAlignPoint->isLocalDerivativeDone())
    return false;

  // clean slate for the local equation for this measurement

  resetGlocalDerivative();
  resetLocalDerivative();

  bool success = true;

  // local derivatives
  // index [0 .. 3] for {dX0, dTx, dY0, dTz}

  success &= setLocalDerivative(0, mAlignPoint->localDerivativeX().dX0());
  success &= setLocalDerivative(1, mAlignPoint->localDerivativeX().dTx());
  success &= setLocalDerivative(2, mAlignPoint->localDerivativeX().dY0());
  success &= setLocalDerivative(3, mAlignPoint->localDerivativeX().dTy());

  // global derivatives
  // index [0 .. 3] for {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}

  Int_t chipId = mAlignPoint->getSensorId();
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 0, mAlignPoint->globalDerivativeX().dDeltaX());
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 1, mAlignPoint->globalDerivativeX().dDeltaY());
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 2, mAlignPoint->globalDerivativeX().dDeltaRz());
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 3, mAlignPoint->globalDerivativeX().dDeltaZ());

  if (success) {
    if (mCounterUsedTracks < 5)
      LOGF(info,
           "setLocalEquationX(): track %i local %.3e %.3e %.3e %.3e, global %.3e %.3e %.3e %.3e X %.3e",
           mCounterUsedTracks,
           mLocalDerivatives[0], mLocalDerivatives[1], mLocalDerivatives[2], mLocalDerivatives[3],
           mGlobalDerivatives[0], mGlobalDerivatives[1], mGlobalDerivatives[2], mGlobalDerivatives[3],
           mAlignPoint->getLocalMeasuredPosition().X());

    mMillepede->SetLocalEquation(
      mGlobalDerivatives,
      mLocalDerivatives,
      mAlignPoint->getLocalMeasuredPosition().X(),
      mAlignPoint->getMeasuredPositionSigma().X());
  } else {
    mCounterLocalEquationFailed++;
  }

  return success;
}

//__________________________________________________________________________
bool AlignHelper::setLocalEquationY()
{
  if (!mAlignPoint->isAlignPointSet())
    return false;
  if (!mAlignPoint->isGlobalDerivativeDone())
    return false;
  if (!mAlignPoint->isLocalDerivativeDone())
    return false;

  // clean slate for the local equation for this measurement

  resetGlocalDerivative();
  resetLocalDerivative();

  bool success = true;

  // local derivatives
  // index [0 .. 3] for {dX0, dTx, dY0, dTz}

  success &= setLocalDerivative(0, mAlignPoint->localDerivativeY().dX0());
  success &= setLocalDerivative(1, mAlignPoint->localDerivativeY().dTx());
  success &= setLocalDerivative(2, mAlignPoint->localDerivativeY().dY0());
  success &= setLocalDerivative(3, mAlignPoint->localDerivativeY().dTy());

  // global derivatives
  // index [0 .. 3] for {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}

  Int_t chipId = mAlignPoint->getSensorId();
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 0, mAlignPoint->globalDerivativeY().dDeltaX());
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 1, mAlignPoint->globalDerivativeY().dDeltaY());
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 2, mAlignPoint->globalDerivativeY().dDeltaRz());
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 3, mAlignPoint->globalDerivativeY().dDeltaZ());

  if (success) {
    if (mCounterUsedTracks < 5)
      LOGF(info,
           "setLocalEquationY(): track %i local %.3e %.3e %.3e %.3e, global %.3e %.3e %.3e %.3e Y %.3e",
           mCounterUsedTracks,
           mLocalDerivatives[0], mLocalDerivatives[1], mLocalDerivatives[2], mLocalDerivatives[3],
           mGlobalDerivatives[0], mGlobalDerivatives[1], mGlobalDerivatives[2], mGlobalDerivatives[3],
           mAlignPoint->getLocalMeasuredPosition().Y());

    mMillepede->SetLocalEquation(
      mGlobalDerivatives,
      mLocalDerivatives,
      mAlignPoint->getLocalMeasuredPosition().Y(),
      mAlignPoint->getMeasuredPositionSigma().Y());
  } else {
    mCounterLocalEquationFailed++;
  }

  return success;
}

//__________________________________________________________________________
bool AlignHelper::setLocalEquationZ()
{

  if (!mAlignPoint->isAlignPointSet())
    return false;
  if (!mAlignPoint->isGlobalDerivativeDone())
    return false;
  if (!mAlignPoint->isLocalDerivativeDone())
    return false;

  // clean slate for the local equation for this measurement

  resetGlocalDerivative();
  resetLocalDerivative();

  bool success = true;

  // local derivatives
  // index [0 .. 3] for {dX0, dTx, dY0, dTz}

  success &= setLocalDerivative(0, mAlignPoint->localDerivativeZ().dX0());
  success &= setLocalDerivative(1, mAlignPoint->localDerivativeZ().dTx());
  success &= setLocalDerivative(2, mAlignPoint->localDerivativeZ().dY0());
  success &= setLocalDerivative(3, mAlignPoint->localDerivativeZ().dTy());

  // global derivatives
  // index [0 .. 3] for {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}

  Int_t chipId = mAlignPoint->getSensorId();
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 0, mAlignPoint->globalDerivativeZ().dDeltaX());
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 1, mAlignPoint->globalDerivativeZ().dDeltaY());
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 2, mAlignPoint->globalDerivativeZ().dDeltaRz());
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 3, mAlignPoint->globalDerivativeZ().dDeltaZ());

  if (success) {
    if (mCounterUsedTracks < 5)
      LOGF(info,
           "setLocalEquationZ(): track %i local %.3e %.3e %.3e %.3e, global %.3e %.3e %.3e %.3e Z %.3e",
           mCounterUsedTracks,
           mLocalDerivatives[0], mLocalDerivatives[1], mLocalDerivatives[2], mLocalDerivatives[3],
           mGlobalDerivatives[0], mGlobalDerivatives[1], mGlobalDerivatives[2], mGlobalDerivatives[3],
           mAlignPoint->getLocalMeasuredPosition().Y());
    mMillepede->SetLocalEquation(
      mGlobalDerivatives,
      mLocalDerivatives,
      mAlignPoint->getLocalMeasuredPosition().Z(),
      mAlignPoint->getMeasuredPositionSigma().Z());
  } else {
    mCounterLocalEquationFailed++;
  }

  return success;
}
