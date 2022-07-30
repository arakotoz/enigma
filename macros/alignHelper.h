#include <array>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <Rtypes.h>
#include <TString.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

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

struct AlignPoint {
  UShort_t sensor;          // sensor id
  UShort_t layer;           // layer id
  UShort_t disk;            // disk id
  UShort_t half;            // half id
  Double_t measuredGlobalX; // cluster x, global frame (cm)
  Double_t measuredGlobalY; // cluster y, global frame (cm)
  Double_t measuredGlobalZ; // cluster z, global frame (cm)
  Double_t measuredLocalX;  // cluster x, local frame (cm)
  Double_t measuredLocalY;  // cluster y, local frame (cm)
  Double_t measuredLocalZ;  // cluster z, local frame (cm)
  Double_t residualX;       // track global x - cluster global x (cm)
  Double_t residualY;       // track global y - cluster global y (cm)
  Double_t residualZ;       // track global z - cluster global z (cm)
  Double_t residualLocalX;  // track local x - cluster local x (cm)
  Double_t residualLocalY;  // track local y - cluster local y (cm)
  Double_t residualLocalZ;  // track local z - cluster local z (cm)
  Double_t recoGlobalX;     // track x, global frame (cm)
  Double_t recoGlobalY;     // track y, global frame (cm)
  Double_t recoGlobalZ;     // track z, global frame (cm)
  Double_t recoLocalX;      // track x, local frame (cm)
  Double_t recoLocalY;      // track y, local frame (cm)
  Double_t recoLocalZ;      // track z, local frame (cm)
};

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
    Double_t chi2CutFactor = 256.;      ///< used to reject outliers i.e. bad tracks with sum(chi2) > Chi2DoFLim(fNStdDev, nDoF) * fChi2CutFactor
  };

 public:
  /// \brief construtor
  AlignHelper();

  /// \brief destructor
  ~AlignHelper();

  /// \brief init Millipede and AlignPointHelper
  void init();

  // simple setters

  void setClusterDictionary(const o2::itsmft::TopologyDictionary* d) { mDictionary = d; }
  void setRunNumber(const int value) { mRunNumber = value; }
  void setBz(const float bz) { mBz = bz; }
  void setChi2CutNStdDev(const Int_t value) { mChi2CutNStdDev = value; }
  void setResidualCutInitial(const Double_t value) { mResCutInitial = value; }
  void setResidualCut(const Double_t value) { mResCut = value; }
  void setMinNumberClusterCut(const int value) { mMinNumberClusterCut = value; }
  void setAllowedVariationDeltaX(const double value) { mAllowVar[0] = value; }
  void setAllowedVariationDeltaY(const double value) { mAllowVar[1] = value; }
  void setAllowedVariationDeltaZ(const double value) { mAllowVar[3] = value; }
  void setAllowedVariationDeltaRz(const double value) { mAllowVar[2] = value; }
  void setChi2CutFactor(const double value) { mStartFac = value; }

  /// \brief  access mft tracks and clusters in the ROOT files, process all ROFs
  void processROFs(TChain* mfttrackChain, TChain* mftclusterChain);

  /// \brief use valid tracks to build Mille records
  // void processRecoTracks();

  /// \brief perform the simultaneous fit of track and alignement parameters
  void globalFit();

  /// \brief print a summary status of what happened in processRecoTracks()
  void printProcessTrackSummary();

  /// \brief provide access to the AlignParam vector
  void getAlignParams(std::vector<o2::detectors::AlignParam>& alignParams) { alignParams = mAlignParams; }

  void initControlTree();
  void closeControlTree();

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
  double* mGlobalDerivatives = nullptr;                                          ///< Array of global derivatives {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}
  double* mLocalDerivatives = nullptr;                                           ///< Array of local derivatives {dX0, dTx, dY0, dTz}
  std::array<Double_t, mNDofPerSensor> mAllowVar;                                ///< "Encouraged" variation for degrees of freedom {dx, dy, dRz, dz}
  double mStartFac = 256;                                                        ///< Initial value for chi2 cut (if > 1, iterations in Millepede are turned on)
  Int_t mChi2CutNStdDev = 3;                                                     ///< Number of standard deviations for chi2 cut
  Double_t mResCutInitial = 100.;                                                ///< Cut on residual on first iteration
  Double_t mResCut = 100.;                                                       ///< Cut on residual for other iterations
  int mMinNumberClusterCut = 6;                                                  ///< Minimum number of clusters in the track to be used for alignment
  o2::mft::MillePedeRecord mTrackRecord;                                         ///< running MillePede Track record
  double mWeightRecord = 1.;                                                     ///< the weight given to a single Mille record in Millepede algorithm
  TString mMilleRecordsFileName;                                                 ///< output file name when saving the Mille records
  TString mMilleConstraintsRecFileName;                                          ///< output file name when saving the records of the constraints
  std::unique_ptr<o2::mft::MillePede2> mMillepede = nullptr;                     ///< Millepede2 implementation copied from AliROOT
  const o2::itsmft::TopologyDictionary* mDictionary = nullptr;                   ///< cluster patterns dictionary
  std::unique_ptr<o2::mft::AlignPointHelper> mAlignPoint = nullptr;              ///< AlignHelper point helper
  std::vector<o2::detectors::AlignParam> mAlignParams;                           ///< vector of alignment parameters computed by Millepede global fit
  bool mIsInitDone = false;                                                      ///< boolean to follow the initialisation status
  int* mGlobalParameterStatus = nullptr;                                         ///< Array of effective degrees of freedom, used to fix detectors, parameters, etc.

  // used to fix some degrees of freedom

  static constexpr Int_t mFixedParId = -1;
  static constexpr Int_t mFreeParId = mFixedParId - 1;

  TFile* mFile;
  TTree* mTree;
  AlignPoint mPointInfo;

  void fillTree();

  /// \brief set array of local derivatives
  bool setLocalDerivative(int index, double value);

  /// \brief set array of global derivatives
  bool setGlobalDerivative(int index, double value);

  /// \brief reset the array of the Local derivative
  bool resetLocalDerivative();

  /// \brief reset the array of the Global derivative
  bool resetGlocalDerivative();

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
    mGlobalDerivatives(nullptr),
    mLocalDerivatives(nullptr),
    mStartFac(256),
    mChi2CutNStdDev(3),
    mResCutInitial(100.),
    mResCut(100.),
    mMinNumberClusterCut(6),
    mWeightRecord(1.),
    mMilleRecordsFileName("mft_mille_records.root"),
    mMilleConstraintsRecFileName("mft_mille_constraints.root"),
    mIsInitDone(false),
    mGlobalParameterStatus(nullptr),
    mFile(nullptr),
    mTree(nullptr)
{
  // default allowed variations w.r.t. global system coordinates
  mAllowVar[0] = 0.5;  // delta translation in x (cm)
  mAllowVar[1] = 0.5;  // delta translation in y (cm)
  mAllowVar[2] = 0.01; // rotation angle Rz around z-axis (rad)
  mAllowVar[3] = 0.5;  // delta translation in z (cm)

  mGlobalDerivatives = (double*)malloc(sizeof(double) * mNumberOfGlobalParam);
  mLocalDerivatives = new Double_t[mNumberOfTrackParam];

  // initialise the content of each array
  resetGlocalDerivative();
  resetLocalDerivative();

  mGlobalParameterStatus = (int*)malloc(sizeof(int) * mNumberOfGlobalParam);
  for (int iPar = 0; iPar < mNumberOfGlobalParam; iPar++) {
    mGlobalParameterStatus[iPar] = mFreeParId;
  }
  mPointInfo.sensor = 0;
  mPointInfo.layer = 0;
  mPointInfo.disk = 0;
  mPointInfo.half = 0;
  mPointInfo.measuredGlobalX = 0;
  mPointInfo.measuredGlobalY = 0;
  mPointInfo.measuredGlobalZ = 0;
  mPointInfo.measuredLocalX = 0;
  mPointInfo.measuredLocalY = 0;
  mPointInfo.measuredLocalZ = 0;
  mPointInfo.residualX = 0;
  mPointInfo.residualY = 0;
  mPointInfo.residualZ = 0;
  mPointInfo.residualLocalX = 0;
  mPointInfo.residualLocalY = 0;
  mPointInfo.residualLocalZ = 0;
  mPointInfo.recoGlobalX = 0;
  mPointInfo.recoGlobalY = 0;
  mPointInfo.recoGlobalZ = 0;
  mPointInfo.recoLocalX = 0;
  mPointInfo.recoLocalY = 0;
  mPointInfo.recoLocalZ = 0;
  LOGF(info, "AlignHelper instantiated");
}

//__________________________________________________________________________
AlignHelper::~AlignHelper()
{
  free(mGlobalDerivatives);
  delete[] mLocalDerivatives;
  free(mGlobalParameterStatus);
}

//__________________________________________________________________________
void AlignHelper::init()
{
  if (mIsInitDone)
    return;
  if (mDictionary == nullptr) {
    LOGF(fatal, "Alignment::init() failed because no cluster dictionary is defined");
    mIsInitDone = false;
    return;
  }
  mMillepede = std::make_unique<MillePede2>();
  mAlignPoint = std::make_unique<AlignPointHelper>();
  mAlignPoint->setClusterDictionary(mDictionary);
  mMillepede->InitMille(mNumberOfGlobalParam,
                        mNumberOfTrackParam,
                        mChi2CutNStdDev,
                        mResCut,
                        mResCutInitial);
  // filenames for the processed tracks and constraints records
  mMillepede->SetDataRecFName(mMilleRecordsFileName.Data());
  mMillepede->SetConsRecFName(mMilleConstraintsRecFileName.Data());

  bool read = false;
  mMillepede->InitDataRecStorage(read);
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

  initControlTree();

  LOGF(info, "AlignHelper init done");
  mIsInitDone = true;
}

//__________________________________________________________________________
void AlignHelper::processROFs(TChain* mfttrackChain, TChain* mftclusterChain)
{
  if (!mIsInitDone) {
    LOGF(fatal, "AlignHelper::processROFs() aborted because init was not done!");
    return;
  }

  LOG(info) << "AlignHelper::processROFs() - start";

  TTreeReader mftTrackChainReader(mfttrackChain);
  TTreeReader mftClusterChainReader(mftclusterChain);
  std::vector<unsigned char>::iterator pattIt;

  TTreeReaderValue<std::vector<o2::mft::TrackMFT>> mMFTTracks =
    {mftTrackChainReader, "MFTTrack"};
  TTreeReaderValue<std::vector<o2::itsmft::ROFRecord>> mMFTTracksROF =
    {mftTrackChainReader, "MFTTracksROF"};
  TTreeReaderValue<std::vector<int>> mMFTTrackClusIdx =
    {mftTrackChainReader, "MFTTrackClusIdx"};

  TTreeReaderValue<std::vector<o2::itsmft::CompClusterExt>> mMFTClusters =
    {mftClusterChainReader, "MFTClusterComp"};
  TTreeReaderValue<std::vector<o2::itsmft::ROFRecord>> mMFTClustersROF =
    {mftClusterChainReader, "MFTClustersROF"};
  TTreeReaderValue<std::vector<unsigned char>> mMFTClusterPatterns =
    {mftClusterChainReader, "MFTClusterPatt"};

  bool firstEntry = true;
  while (mftTrackChainReader.Next() && mftClusterChainReader.Next()) {

    if (firstEntry)
      pattIt = (*mMFTClusterPatterns).begin();

    mNumberOfTrackChainROFs += (*mMFTTracksROF).size();
    mNumberOfClusterChainROFs += (*mMFTClustersROF).size();
    assert(mNumberOfTrackChainROFs == mNumberOfClusterChainROFs);

    for (const auto& oneTrack : *mMFTTracks) { // track loop

      // Skip the track if not enough clusters
      auto ncls = oneTrack.getNumberOfPoints();
      if (ncls < mMinNumberClusterCut) {
        mCounterSkippedTracks++;
        continue;
      }

      // Skip presumably quite low momentum track
      if (!oneTrack.isLTF()) {
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

        mAlignPoint->resetAlignPoint();

        // Store measured positions
        auto clsEntry = (*mMFTTrackClusIdx)[offset + icls];
        const auto compCluster = (*mMFTClusters)[clsEntry];
        mAlignPoint->setMeasuredPosition(compCluster, pattIt);

        // Propagate track to the current z plane of this cluster
        auto track = oneTrack;
        track.propagateParamToZlinear(mAlignPoint->getGlobalMeasuredPosition().Z());

        // Store reco positions
        mAlignPoint->setGlobalRecoPosition(track);

        // compute residuals
        mAlignPoint->setLocalResidual();
        mAlignPoint->setGlobalResidual();
        fillTree();

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
        mMillepede->SaveRecordData();
        mCounterUsedTracks++;
      }
    } // end of loop on tracks

    firstEntry = false;

  } // end of loop on TChain reader

  LOG(info) << "AlignHelper::processROFs() - end";
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

  double* params = (double*)malloc(sizeof(double) * mNumberOfGlobalParam);
  double* paramsErrors = (double*)malloc(sizeof(double) * mNumberOfGlobalParam);
  double* paramsPulls = (double*)malloc(sizeof(double) * mNumberOfGlobalParam);

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

  AlignSensorHelper chipHelper;
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

  free(params);
  free(paramsErrors);
  free(paramsPulls);
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
void AlignHelper::initControlTree()
{
  mFile = TFile::Open("align_point.root", "recreate", "", 505);

  mTree = new TTree("point", "the align point info tree");
  mTree->Branch("sensor", &mPointInfo.sensor, "sensor/s");
  mTree->Branch("layer", &mPointInfo.layer, "layer/s");
  mTree->Branch("disk", &mPointInfo.disk, "disk/s");
  mTree->Branch("half", &mPointInfo.half, "half/s");
  mTree->Branch("measuredGlobalX", &mPointInfo.measuredGlobalX, "measuredGlobalX/D");
  mTree->Branch("measuredGlobalY", &mPointInfo.measuredGlobalY, "measuredGlobalY/D");
  mTree->Branch("measuredGlobalZ", &mPointInfo.measuredGlobalZ, "measuredGlobalZ/D");
  mTree->Branch("measuredLocalX", &mPointInfo.measuredLocalX, "measuredLocalX/D");
  mTree->Branch("measuredLocalY", &mPointInfo.measuredLocalY, "measuredLocalY/D");
  mTree->Branch("measuredLocalZ", &mPointInfo.measuredLocalZ, "measuredLocalZ/D");
  mTree->Branch("residualX", &mPointInfo.residualX, "residualX/D");
  mTree->Branch("residualY", &mPointInfo.residualY, "residualY/D");
  mTree->Branch("residualZ", &mPointInfo.residualZ, "residualZ/D");
  mTree->Branch("residualLocalX", &mPointInfo.residualLocalX, "residualLocalX/D");
  mTree->Branch("residualLocalY", &mPointInfo.residualLocalY, "residualLocalY/D");
  mTree->Branch("residualLocalZ", &mPointInfo.residualLocalZ, "residualLocalZ/D");
  mTree->Branch("recoGlobalX", &mPointInfo.recoGlobalX, "recoGlobalX/D");
  mTree->Branch("recoGlobalY", &mPointInfo.recoGlobalY, "recoGlobalY/D");
  mTree->Branch("recoGlobalZ", &mPointInfo.recoGlobalZ, "recoGlobalZ/D");
  mTree->Branch("recoLocalX", &mPointInfo.recoLocalX, "recoLocalX/D");
  mTree->Branch("recoLocalY", &mPointInfo.recoLocalY, "recoLocalY/D");
  mTree->Branch("recoLocalZ", &mPointInfo.recoLocalZ, "recoLocalZ/D");
}

//__________________________________________________________________________
void AlignHelper::closeControlTree()
{
  if (mTree) {
    if (mFile && mFile->IsWritable()) {
      mFile->cd();
      mTree->Write();
    }
    delete mTree;
    if (mFile) {
      mFile->Close();
      delete mFile;
    }
  }
  LOG(info) << "Closed file align_point.root";
}

//__________________________________________________________________________
bool AlignHelper::setLocalDerivative(int index, double value)
{
  // index [0 .. 3] for {dX0, dTx, dY0, dTz}

  bool success = false;
  if (index < mNumberOfTrackParam) {
    mLocalDerivatives[index] = value;
    success = true;
  } else {
    LOGF(error,
         "AlignHelper::setLocalDerivative() - index %d >= %d",
         index, mNumberOfTrackParam);
  }
  return success;
}

//__________________________________________________________________________
bool AlignHelper::setGlobalDerivative(int index, double value)
{
  // index [0 .. 3] for {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}

  bool success = false;
  if (index < mNumberOfGlobalParam) {
    mGlobalDerivatives[index] = value;
    success = true;
  } else {
    LOGF(error,
         "AlignHelper::setGlobalDerivative() - index %d >= %d",
         index, mNumberOfGlobalParam);
  }
  return success;
}

//__________________________________________________________________________
bool AlignHelper::resetLocalDerivative()
{
  bool success = false;
  for (int i = 0; i < mNumberOfTrackParam; ++i) {
    success = false;
    mLocalDerivatives[i] = 0.0;
    success = true;
  }
  return success;
}

//__________________________________________________________________________
bool AlignHelper::resetGlocalDerivative()
{
  bool success = false;
  for (int i = 0; i < mNumberOfGlobalParam; ++i) {
    success = false;
    mGlobalDerivatives[i] = 0.0;
    success = true;
  }
  return success;
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

  bool success = true;

  // clean slate for the local equation for this measurement

  success &= resetGlocalDerivative();
  success &= resetLocalDerivative();

  // local derivatives
  // index [0 .. 3] for {dX0, dTx, dY0, dTz}

  success &= setLocalDerivative(0, mAlignPoint->localDerivativeX().dX0());
  success &= setLocalDerivative(1, mAlignPoint->localDerivativeX().dTx());
  success &= setLocalDerivative(2, mAlignPoint->localDerivativeX().dY0());
  success &= setLocalDerivative(3, mAlignPoint->localDerivativeX().dTy());

  // global derivatives
  // index [0 .. 3] for {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}

  Int_t chipId = mAlignPoint->getSensorId();
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 0, mAlignPoint->globalDerivativeX().dDeltaX());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 1, mAlignPoint->globalDerivativeX().dDeltaY());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 2, mAlignPoint->globalDerivativeX().dDeltaRz());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 3, mAlignPoint->globalDerivativeX().dDeltaZ());

  if (success) {
    bool debugPrint = false;
    if (mCounterUsedTracks < 5) {
      LOGF(info,
           "setLocalEquationX(): track %i sr %4d local %.3e %.3e %.3e %.3e, global %.3e %.3e %.3e %.3e X %.3e sigma %.3e",
           mCounterUsedTracks, chipId,
           mLocalDerivatives[0], mLocalDerivatives[1], mLocalDerivatives[2], mLocalDerivatives[3],
           mGlobalDerivatives[chipId * mNDofPerSensor + 0],
           mGlobalDerivatives[chipId * mNDofPerSensor + 1],
           mGlobalDerivatives[chipId * mNDofPerSensor + 2],
           mGlobalDerivatives[chipId * mNDofPerSensor + 3],
           mAlignPoint->getLocalMeasuredPosition().X(),
           mAlignPoint->getLocalMeasuredPositionSigma().X());
      debugPrint = true;
    }

    mMillepede->SetLocalEquation(
      mGlobalDerivatives,
      mLocalDerivatives,
      mAlignPoint->getLocalMeasuredPosition().X(),
      mAlignPoint->getLocalMeasuredPositionSigma().X(),
      debugPrint);
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

  bool success = true;

  // clean slate for the local equation for this measurement

  success &= resetGlocalDerivative();
  success &= resetLocalDerivative();

  // local derivatives
  // index [0 .. 3] for {dX0, dTx, dY0, dTz}

  success &= setLocalDerivative(0, mAlignPoint->localDerivativeY().dX0());
  success &= setLocalDerivative(1, mAlignPoint->localDerivativeY().dTx());
  success &= setLocalDerivative(2, mAlignPoint->localDerivativeY().dY0());
  success &= setLocalDerivative(3, mAlignPoint->localDerivativeY().dTy());

  // global derivatives
  // index [0 .. 3] for {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}

  Int_t chipId = mAlignPoint->getSensorId();
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 0, mAlignPoint->globalDerivativeY().dDeltaX());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 1, mAlignPoint->globalDerivativeY().dDeltaY());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 2, mAlignPoint->globalDerivativeY().dDeltaRz());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 3, mAlignPoint->globalDerivativeY().dDeltaZ());

  if (success) {
    bool debugPrint = false;
    if (mCounterUsedTracks < 5) {
      LOGF(info,
           "setLocalEquationY(): track %i sr %4d local %.3e %.3e %.3e %.3e, global %.3e %.3e %.3e %.3e Y %.3e sigma %.3e",
           mCounterUsedTracks, chipId,
           mLocalDerivatives[0], mLocalDerivatives[1], mLocalDerivatives[2], mLocalDerivatives[3],
           mGlobalDerivatives[chipId * mNDofPerSensor + 0],
           mGlobalDerivatives[chipId * mNDofPerSensor + 1],
           mGlobalDerivatives[chipId * mNDofPerSensor + 2],
           mGlobalDerivatives[chipId * mNDofPerSensor + 3],
           mAlignPoint->getLocalMeasuredPosition().Y(),
           mAlignPoint->getLocalMeasuredPositionSigma().Y());
      debugPrint = true;
    }

    mMillepede->SetLocalEquation(
      mGlobalDerivatives,
      mLocalDerivatives,
      mAlignPoint->getLocalMeasuredPosition().Y(),
      mAlignPoint->getLocalMeasuredPositionSigma().Y(),
      debugPrint);
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

  bool success = true;

  // clean slate for the local equation for this measurement

  success &= resetGlocalDerivative();
  success &= resetLocalDerivative();

  // local derivatives
  // index [0 .. 3] for {dX0, dTx, dY0, dTz}

  success &= setLocalDerivative(0, mAlignPoint->localDerivativeZ().dX0());
  success &= setLocalDerivative(1, mAlignPoint->localDerivativeZ().dTx());
  success &= setLocalDerivative(2, mAlignPoint->localDerivativeZ().dY0());
  success &= setLocalDerivative(3, mAlignPoint->localDerivativeZ().dTy());

  // global derivatives
  // index [0 .. 3] for {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}

  Int_t chipId = mAlignPoint->getSensorId();
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 0, mAlignPoint->globalDerivativeZ().dDeltaX());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 1, mAlignPoint->globalDerivativeZ().dDeltaY());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 2, mAlignPoint->globalDerivativeZ().dDeltaRz());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 3, mAlignPoint->globalDerivativeZ().dDeltaZ());

  if (success) {
    bool debugPrint = false;
    if (mCounterUsedTracks) {
      LOGF(info,
           "setLocalEquationZ(): track %i sr %4d local %.3e %.3e %.3e %.3e, global %.3e %.3e %.3e %.3e Z %.3e sigma %.3e",
           mCounterUsedTracks, chipId,
           mLocalDerivatives[0], mLocalDerivatives[1], mLocalDerivatives[2], mLocalDerivatives[3],
           mGlobalDerivatives[chipId * mNDofPerSensor + 0],
           mGlobalDerivatives[chipId * mNDofPerSensor + 1],
           mGlobalDerivatives[chipId * mNDofPerSensor + 2],
           mGlobalDerivatives[chipId * mNDofPerSensor + 3],
           mAlignPoint->getLocalMeasuredPosition().Z(),
           mAlignPoint->getLocalMeasuredPositionSigma().Z());
      debugPrint = true;
    }

    mMillepede->SetLocalEquation(
      mGlobalDerivatives,
      mLocalDerivatives,
      mAlignPoint->getLocalMeasuredPosition().Z(),
      mAlignPoint->getLocalMeasuredPositionSigma().Z(),
      debugPrint);
  } else {
    mCounterLocalEquationFailed++;
  }

  return success;
}

//__________________________________________________________________________
void AlignHelper::fillTree()
{
  if (mTree) {
    mPointInfo.sensor = mAlignPoint->getSensorId();
    mPointInfo.layer = mAlignPoint->layer();
    mPointInfo.disk = mAlignPoint->disk();
    mPointInfo.half = mAlignPoint->half();
    mPointInfo.measuredGlobalX = mAlignPoint->getGlobalMeasuredPosition().X();
    mPointInfo.measuredGlobalY = mAlignPoint->getGlobalMeasuredPosition().Y();
    mPointInfo.measuredGlobalZ = mAlignPoint->getGlobalMeasuredPosition().Z();
    mPointInfo.measuredLocalX = mAlignPoint->getLocalMeasuredPosition().X();
    mPointInfo.measuredLocalY = mAlignPoint->getLocalMeasuredPosition().Y();
    mPointInfo.measuredLocalZ = mAlignPoint->getLocalMeasuredPosition().Z();
    mPointInfo.residualX = mAlignPoint->getGlobalResidual().X();
    mPointInfo.residualY = mAlignPoint->getGlobalResidual().Y();
    mPointInfo.residualZ = mAlignPoint->getGlobalResidual().Z();
    mPointInfo.residualLocalX = mAlignPoint->getLocalResidual().X();
    mPointInfo.residualLocalY = mAlignPoint->getLocalResidual().Y();
    mPointInfo.residualLocalZ = mAlignPoint->getLocalResidual().Z();
    mPointInfo.recoGlobalX = mAlignPoint->getGlobalRecoPosition().X();
    mPointInfo.recoGlobalY = mAlignPoint->getGlobalRecoPosition().Y();
    mPointInfo.recoGlobalZ = mAlignPoint->getGlobalRecoPosition().Z();
    mPointInfo.recoLocalX = mAlignPoint->getLocalRecoPosition().X();
    mPointInfo.recoLocalY = mAlignPoint->getLocalRecoPosition().Y();
    mPointInfo.recoLocalZ = mAlignPoint->getLocalRecoPosition().Z();

    if (mCounterUsedTracks < 5) {
      LOGF(info, "track %i h %d d %d l %d s %4d lMpos x %.2e y %.2e z %.2e gMpos x %.2e y %.2e z %.2e",
           mCounterUsedTracks, mPointInfo.half, mPointInfo.disk, mPointInfo.layer, mPointInfo.sensor,
           mPointInfo.measuredLocalX, mPointInfo.measuredLocalY, mPointInfo.measuredLocalZ,
           mPointInfo.measuredGlobalX, mPointInfo.measuredGlobalY, mPointInfo.measuredGlobalZ);
    }
    mTree->Fill();
  }
}