#include <bitset>
#include <iomanip>
#include <iostream>
#include <string>

#include <TObject.h>
#include <Rtypes.h>

#include "DataFormatsITSMFT/ClusterPattern.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "MathUtils/Cartesian.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"
#include "MFTTracking/IOUtils.h"

struct HitStruct {
  UInt_t rofIdx;            // readout frame index
  UInt_t bc;                // BC
  UInt_t orbit;             // orbit
  UShort_t sensor;          // sensor id
  UShort_t layer;           // layer id
  UShort_t disk;            // disk id
  UShort_t half;            // half id
  Int_t trackIdx;           // if attached to a track, track index > -1
  Double_t measuredGlobalX; // cluster x, global frame (cm)
  Double_t measuredGlobalY; // cluster y, global frame (cm)
  Double_t measuredGlobalZ; // cluster z, global frame (cm)
  Double_t measuredLocalX;  // cluster x, local frame (cm)
  Double_t measuredLocalY;  // cluster y, local frame (cm)
  Double_t measuredLocalZ;  // cluster z, local frame (cm)
  Double_t measuredSigmaX2; // cluster, variance x (cm2)
  Double_t measuredSigmaY2; // cluster, variance y (cm2)
  Double_t measuredSigmaZ2; // cluster, variance z (cm2)
  Double_t recoGlobalX;     // track x, global frame (cm)
  Double_t recoGlobalY;     // track y, global frame (cm)
  Double_t recoGlobalZ;     // track z, global frame (cm)
  Double_t recoLocalX;      // track x, local frame (cm)
  Double_t recoLocalY;      // track y, local frame (cm)
  Double_t recoLocalZ;      // track z, local frame (cm)
  Double_t residualX;       // track global x - cluster global x (cm)
  Double_t residualY;       // track global y - cluster global y (cm)
  Double_t residualZ;       // track global z - cluster global z (cm)
  Double_t residualLocalX;  // track local x - cluster local x (cm)
  Double_t residualLocalY;  // track local y - cluster local y (cm)
  Double_t residualLocalZ;  // track local z - cluster local z (cm)
};

class Hit
{
 public:
  enum RefFrame_t {
    isLocal = 0, // Local reference frame
    isGlobal = 1 // global reference frame
  };

  Hit();
  ~Hit() = default;

  // non-default constructors

  Hit(const Int_t sensor,
      o2::itsmft::ChipMappingMFT mapping,
      const o2::math_utils::Point3D<double>& clusterPosition,
      const RefFrame_t refFrame = isGlobal);
  Hit(const Int_t sensor,
      o2::itsmft::ChipMappingMFT mapping,
      const double clusterX,
      const double clusterY,
      const double clusterZ,
      const RefFrame_t refFrame = isGlobal);
  Hit(o2::itsmft::CompClusterExt cluster,
      std::vector<unsigned char>::iterator& pattIt,
      o2::mft::GeometryTGeo* geom,
      o2::itsmft::ChipMappingMFT chipMappingMFT,
      o2::itsmft::TopologyDictionary& dict,
      const UInt_t rof = 0);

  // convert a compact cluster to 3D spacepoint stored into a Hit

  void convertCompactCluster(o2::itsmft::CompClusterExt c,
                             std::vector<unsigned char>::iterator& pattIt,
                             o2::mft::GeometryTGeo* geom,
                             o2::itsmft::ChipMappingMFT chipMappingMFT,
                             o2::itsmft::TopologyDictionary& dict);

  // convert a track 3D spacepoint from glocal to local coordinates

  void convertT2LTrackPosition(UShort_t chipID,
                               o2::mft::GeometryTGeo* geom);

  // short-named getters

  UInt_t bc() const { return mBc; }
  double clusterX(const RefFrame_t refFrame = isGlobal) const
  {
    if (refFrame == isGlobal) {
      return mGlobalMeasuredPosition.X();
    } else {
      return mLocalMeasuredPosition.X();
    }
  }
  double clusterY(const RefFrame_t refFrame = isGlobal) const
  {
    if (refFrame == isGlobal) {
      return mGlobalMeasuredPosition.Y();
    } else {
      return mLocalMeasuredPosition.Y();
    }
  }
  double clusterZ(const RefFrame_t refFrame = isGlobal) const
  {
    if (refFrame == isGlobal) {
      return mGlobalMeasuredPosition.Z();
    } else {
      return mLocalMeasuredPosition.Z();
    }
  }
  double clusterGlobalX() const { return clusterX(isGlobal); }
  double clusterGlobalY() const { return clusterY(isGlobal); }
  double clusterGlobalZ() const { return clusterZ(isGlobal); }
  double clusterSigmaX2() const { return mMeasuredSigmaX2; }
  double clusterSigmaY2() const { return mMeasuredSigmaY2; }
  double clusterSigmaZ2() const { return mMeasuredSigmaZ2; }
  double clusterLocalX() const { return clusterX(isLocal); }
  double clusterLocalY() const { return clusterY(isLocal); }
  double clusterLocalZ() const { return clusterZ(isLocal); }
  UShort_t disk() const { return mDisk; }
  UShort_t half() const { return mHalf; }
  Bool_t isInTrack() const { return (mTrackIdx > -1) ? kTRUE : kFALSE; }
  UShort_t layer() const { return mLayer; }
  UShort_t orbit() const { return mOrbit; }
  double residualX(const RefFrame_t refFrame = isGlobal) const; // (cm)
  double residualY(const RefFrame_t refFrame = isGlobal) const; // (cm)
  double residualZ(const RefFrame_t refFrame = isGlobal) const; // (cm)
  UInt_t rofIdx() const { return mRofIdx; }
  UShort_t sensor() const { return mSensor; }
  double trackX(const RefFrame_t refFrame = isGlobal) const
  {
    if (refFrame == isGlobal) {
      return mGlobalRecoPosition.X();
    } else {
      return mLocalRecoPosition.X();
    }
  }
  double trackY(const RefFrame_t refFrame = isGlobal) const
  {
    if (refFrame == isGlobal) {
      return mGlobalRecoPosition.Y();
    } else {
      return mLocalRecoPosition.Y();
    }
  }
  double trackZ(const RefFrame_t refFrame = isGlobal) const
  {
    if (refFrame == isGlobal) {
      return mGlobalRecoPosition.Z();
    } else {
      return mLocalRecoPosition.Z();
    }
  }
  double trackGlobalX() const { return trackX(isGlobal); }
  double trackGlobalY() const { return trackY(isGlobal); }
  double trackGlobalZ() const { return trackZ(isGlobal); }
  double trackLocalX() const { return trackX(isLocal); }
  double trackLocalY() const { return trackY(isLocal); }
  double trackLocalZ() const { return trackZ(isLocal); }
  Int_t trackIdx() const { return mTrackIdx; }

  // getters

  o2::math_utils::Point3D<double> getClusterXYZ(
    const RefFrame_t refFrame = isGlobal) const
  {
    if (refFrame == isGlobal) {
      return mGlobalMeasuredPosition;
    } else {
      return mLocalMeasuredPosition;
    }
  }
  o2::math_utils::Point3D<double>& getClusterXYZ(
    const RefFrame_t refFrame = isGlobal)
  {
    if (refFrame == isGlobal) {
      return mGlobalMeasuredPosition;
    } else {
      return mLocalMeasuredPosition;
    }
  }
  HitStruct getHitStruct() const;
  o2::math_utils::Point3D<double> getTrackXYZ(
    const RefFrame_t refFrame = isGlobal) const
  {
    if (refFrame == isGlobal) {
      return mGlobalRecoPosition;
    } else {
      return mLocalRecoPosition;
    }
  }
  o2::math_utils::Point3D<double>& getTrackXYZ(
    const RefFrame_t refFrame = isGlobal)
  {
    if (refFrame == isGlobal) {
      return mGlobalRecoPosition;
    } else {
      return mLocalRecoPosition;
    }
  }

  // print

  void print(const RefFrame_t refFrame = isGlobal);

  // setters

  void setBc(const UInt_t bc) { mBc = bc; }
  void setClusterPosition(
    const o2::math_utils::Point3D<double>& position,
    const RefFrame_t refFrame = isGlobal)
  {
    if (refFrame == isGlobal) {
      mGlobalMeasuredPosition = position;
    } else {
      mLocalMeasuredPosition = position;
    }
  }
  void setClusterPosition(const double x,
                          const double y,
                          const double z,
                          const RefFrame_t refFrame = isGlobal)
  {
    setClusterX(x, refFrame);
    setClusterY(y, refFrame);
    setClusterZ(z, refFrame);
  }
  void setClusterX(const double x, const RefFrame_t refFrame = isGlobal)
  {
    if (refFrame == isGlobal) {
      mGlobalMeasuredPosition.SetX(x);
    } else {
      mLocalMeasuredPosition.SetX(x);
    }
  }
  void setClusterY(const double y, const RefFrame_t refFrame = isGlobal)
  {
    if (refFrame == isGlobal) {
      mGlobalMeasuredPosition.SetY(y);
    } else {
      mLocalMeasuredPosition.SetY(y);
    }
  }
  void setClusterZ(const double z, const RefFrame_t refFrame = isGlobal)
  {
    if (refFrame == isGlobal) {
      mGlobalMeasuredPosition.SetZ(z);
    } else {
      mLocalMeasuredPosition.SetZ(z);
    }
  }
  void setMeasuredErrors(const double sx2,
                         const double sy2, const double sz2)
  {
    setMeasuredSigmaX2(sx2);
    setMeasuredSigmaY2(sy2);
    setMeasuredSigmaZ2(sz2);
  }
  void setMeasuredSigmaX2(const double v) { mMeasuredSigmaX2 = v; }
  void setMeasuredSigmaY2(const double v) { mMeasuredSigmaY2 = v; }
  void setMeasuredSigmaZ2(const double v) { mMeasuredSigmaZ2 = v; }
  void setOrbit(const UInt_t v) { mOrbit = v; }
  void setRofIdx(const UInt_t idx) { mRofIdx = idx; }
  void setSensor(Int_t sensor, o2::itsmft::ChipMappingMFT mapping);
  void setTrackPosition(
    const o2::math_utils::Point3D<double>& position,
    const RefFrame_t refFrame = isGlobal)
  {
    if (refFrame == isGlobal) {
      mGlobalRecoPosition = position;
    } else {
      mLocalRecoPosition = position;
    }
  }
  void setTrackPosition(const double x,
                        const double y,
                        const double z,
                        const RefFrame_t refFrame = isGlobal)
  {
    setTrackX(x, refFrame);
    setTrackY(y, refFrame);
    setTrackZ(z, refFrame);
  }
  void setTrackX(const double x, const RefFrame_t refFrame = isGlobal)
  {
    if (refFrame == isGlobal) {
      mGlobalRecoPosition.SetX(x);
    } else {
      mLocalRecoPosition.SetX(x);
    }
  }
  void setTrackY(const double y, const RefFrame_t refFrame = isGlobal)
  {
    if (refFrame == isGlobal) {
      mGlobalRecoPosition.SetY(y);
    } else {
      mLocalRecoPosition.SetY(y);
    }
  }
  void setTrackZ(const double z, const RefFrame_t refFrame = isGlobal)
  {
    if (refFrame == isGlobal) {
      mGlobalRecoPosition.SetZ(z);
    } else {
      mLocalRecoPosition.SetZ(z);
    }
  }
  void setTrackIdx(const Int_t idx) { mTrackIdx = idx; }

 protected:
  // readout frame index
  UInt_t mRofIdx = 0;
  // BC
  UInt_t mBc = 0;
  // orbit
  UInt_t mOrbit = 0;
  // sensor id
  UShort_t mSensor = 0;
  // layer id
  UShort_t mLayer = 0;
  // disk id
  UShort_t mDisk = 0;
  // half id
  UShort_t mHalf = 0;
  // track index
  Int_t mTrackIdx = -1;
  // Cartesian position (cm, in Global frame) of the reconstructed track
  // analytically propagated to the z position of the cluster
  o2::math_utils::Point3D<double> mGlobalRecoPosition;
  // Cartesian position (cm, in Global frame) of the cluster
  o2::math_utils::Point3D<double> mGlobalMeasuredPosition;
  // Errors on the cluster cartesian position (cm2)
  double mMeasuredSigmaX2 = 0.;
  double mMeasuredSigmaY2 = 0.;
  double mMeasuredSigmaZ2 = 0.;
  // Cartesian position (cm, in Local frame) of the reconstructed track
  // analytically propagated to the z position of the cluster
  o2::math_utils::Point3D<double> mLocalRecoPosition;
  // Cartesian position (cm, in Local frame) of the cluster
  o2::math_utils::Point3D<double> mLocalMeasuredPosition;

  ClassDefNV(Hit, 1);
};

//__________________________________________________________________________
Hit::Hit()
  : mRofIdx(0),
    mBc(0),
    mOrbit(0),
    mSensor(0),
    mLayer(0),
    mDisk(0),
    mHalf(0),
    mTrackIdx(-1),
    mMeasuredSigmaX2(0.),
    mMeasuredSigmaY2(0.),
    mMeasuredSigmaZ2(0.)
{
  setClusterPosition(0., 0., 0., isGlobal);
  setTrackPosition(0., 0., 0., isGlobal);
  setClusterPosition(0., 0., 0., isLocal);
  setTrackPosition(0., 0., 0., isLocal);
}

//__________________________________________________________________________
Hit::Hit(const Int_t sensor,
         o2::itsmft::ChipMappingMFT mapping,
         const o2::math_utils::Point3D<double>& clusterPosition,
         const RefFrame_t refFrame)
  : mRofIdx(0),
    mBc(0),
    mOrbit(0),
    mSensor(0),
    mLayer(0),
    mDisk(0),
    mHalf(0),
    mTrackIdx(-1),
    mMeasuredSigmaX2(0.),
    mMeasuredSigmaY2(0.),
    mMeasuredSigmaZ2(0.)
{
  setClusterPosition(0., 0., 0., isGlobal);
  setTrackPosition(0., 0., 0., isGlobal);
  setClusterPosition(0., 0., 0., isLocal);
  setTrackPosition(0., 0., 0., isLocal);

  setSensor(sensor, mapping);
  setClusterPosition(clusterPosition, refFrame);
}

//__________________________________________________________________________
Hit::Hit(const Int_t sensor,
         o2::itsmft::ChipMappingMFT mapping,
         const double clusterX,
         const double clusterY,
         const double clusterZ,
         const RefFrame_t refFrame)
  : mRofIdx(0),
    mBc(0),
    mOrbit(0),
    mSensor(0),
    mLayer(0),
    mDisk(0),
    mHalf(0),
    mTrackIdx(-1),
    mMeasuredSigmaX2(0.),
    mMeasuredSigmaY2(0.),
    mMeasuredSigmaZ2(0.)
{
  setClusterPosition(0., 0., 0., isGlobal);
  setTrackPosition(0., 0., 0., isGlobal);
  setClusterPosition(0., 0., 0., isLocal);
  setTrackPosition(0., 0., 0., isLocal);

  setSensor(sensor, mapping);
  setClusterPosition(clusterX, clusterY, clusterZ, refFrame);
}

//__________________________________________________________________________
Hit::Hit(o2::itsmft::CompClusterExt cluster,
         std::vector<unsigned char>::iterator& pattIt,
         o2::mft::GeometryTGeo* geom,
         o2::itsmft::ChipMappingMFT chipMappingMFT,
         o2::itsmft::TopologyDictionary& dict,
         const UInt_t rof)
  : mRofIdx(rof),
    mBc(0),
    mOrbit(0),
    mSensor(0),
    mLayer(0),
    mDisk(0),
    mHalf(0),
    mTrackIdx(-1),
    mMeasuredSigmaX2(0.),
    mMeasuredSigmaY2(0.),
    mMeasuredSigmaZ2(0.)
{
  setClusterPosition(0., 0., 0., isGlobal);
  setTrackPosition(0., 0., 0., isGlobal);
  setClusterPosition(0., 0., 0., isLocal);
  setTrackPosition(0., 0., 0., isLocal);

  convertCompactCluster(cluster, pattIt, geom, chipMappingMFT, dict);
}

//__________________________________________________________________________
double Hit::residualX(const RefFrame_t refFrame) const
{
  /// return (track x - cluster x) at this layer, x in ref frame (cm)
  double residual = 0;
  if (isInTrack()) {
    residual = trackX(refFrame) - clusterX(refFrame);
  }
  return residual;
}

//__________________________________________________________________________
double Hit::residualY(const RefFrame_t refFrame) const
{
  /// return (track y - cluster y) at this layer, y in ref frame (cm)
  double residual = 0;
  if (isInTrack()) {
    residual = trackY(refFrame) - clusterY(refFrame);
  }
  return residual;
}

//__________________________________________________________________________
double Hit::residualZ(const RefFrame_t refFrame) const
{
  /// return (track z - cluster z) at this layer, z in ref frame (cm)
  double residual = 0;
  if (isInTrack()) {
    residual = trackZ(refFrame) - clusterZ(refFrame);
  }
  return residual;
}

//__________________________________________________________________________
HitStruct Hit::getHitStruct() const
{
  HitStruct hit{};
  hit.rofIdx = rofIdx();
  hit.bc = bc();
  hit.orbit = orbit();
  hit.sensor = sensor();
  hit.layer = layer();
  hit.disk = disk();
  hit.half = half();
  hit.trackIdx = trackIdx();
  hit.measuredGlobalX = clusterGlobalX();
  hit.measuredGlobalY = clusterGlobalY();
  hit.measuredGlobalZ = clusterGlobalZ();
  hit.measuredLocalX = clusterLocalX();
  hit.measuredLocalY = clusterLocalY();
  hit.measuredLocalZ = clusterLocalZ();
  hit.measuredSigmaX2 = clusterSigmaX2();
  hit.measuredSigmaY2 = clusterSigmaY2();
  hit.measuredSigmaZ2 = clusterSigmaZ2();
  hit.recoGlobalX = trackGlobalX();
  hit.recoGlobalY = trackGlobalY();
  hit.recoGlobalZ = trackGlobalZ();
  hit.recoLocalX = trackLocalX();
  hit.recoLocalY = trackLocalY();
  hit.recoLocalZ = trackLocalZ();
  hit.residualX = residualX(isGlobal);
  hit.residualY = residualY(isGlobal);
  hit.residualZ = residualZ(isGlobal);
  hit.residualLocalX = residualX(isLocal);
  hit.residualLocalY = residualY(isLocal);
  hit.residualLocalZ = residualZ(isLocal);
  return hit;
}

//__________________________________________________________________________
void Hit::print(const RefFrame_t refFrame)
{
  std::cout << "BC " << bc()
            << " orbit " << orbit()
            << " rofIdx " << rofIdx()
            << " hit " << std::noshowpos
            << "s " << std::setw(4) << sensor() << " "
            << "l " << std::setw(2) << layer() << " "
            << "d " << std::setw(2) << disk() << " "
            << "h " << std::setw(1) << half() << " ";
  if (refFrame == isGlobal) {
    std::cout << "(global) ";
  } else {
    std::cout << "(local) ";
  }
  std::cout << "cluster ("
            << std::showpos << std::scientific << std::setprecision(3)
            << clusterX(refFrame) << ", "
            << clusterY(refFrame) << ", "
            << clusterZ(refFrame) << ") "
            << std::noshowpos << std::setprecision(3)
            << "err2 ("
            << clusterSigmaX2() << ", "
            << clusterSigmaY2() << ", "
            << clusterSigmaZ2() << ") ";
  if (isInTrack()) {
    std::cout << "track " << trackIdx() << " ("
              << std::showpos << std::setprecision(3)
              << trackX(refFrame) << ", "
              << trackY(refFrame) << ", "
              << trackZ(refFrame) << ") ";
  }
  std::cout << std::noshowpos << std::setprecision(6)
            << std::endl;
}

//__________________________________________________________________________
void Hit::setSensor(Int_t sensor, o2::itsmft::ChipMappingMFT mapping)
{
  if (sensor < mapping.getNChips()) {
    o2::itsmft::MFTChipMappingData chipMapping = (mapping.getChipMappingData())[sensor];
    mSensor = chipMapping.globalChipSWID;
    mLayer = (UShort_t)chipMapping.layer;
    mDisk = (UShort_t)chipMapping.disk;
    mHalf = (UShort_t)chipMapping.half;
  } else {
    std::cout << "Hit::setSensor() - "
              << "ERROR, sensor id " << sensor
              << " >= " << mapping.getNChips() << std::endl;
  }
}

//_________________________________________________________
void Hit::convertCompactCluster(o2::itsmft::CompClusterExt c,
                                std::vector<unsigned char>::iterator& pattIt,
                                o2::mft::GeometryTGeo* geom,
                                o2::itsmft::ChipMappingMFT chipMappingMFT,
                                o2::itsmft::TopologyDictionary& dict)
{
  /// convert a compact cluster to 3D spacepoint stored into a Hit
  // lines from Detectors/ITSMFT/MFT/tracking/src/IOUtils.cxx
  auto chipID = c.getChipID();
  auto pattID = c.getPatternID();
  o2::math_utils::Point3D<double> locXYZ;
  // Dummy COG errors (about half pixel size)
  double sigmaX2 = o2::mft::ioutils::DefClusError2Row;
  double sigmaY2 = o2::mft::ioutils::DefClusError2Col;
  if (pattID != o2::itsmft::CompCluster::InvalidPatternID) {
    // ALPIDE local Y coordinate => MFT global X coordinate (ALPIDE rows)
    sigmaX2 = dict.getErr2X(pattID);
    // ALPIDE local Z coordinate => MFT global Y coordinate (ALPIDE columns)
    sigmaY2 = dict.getErr2Z(pattID);
    if (!dict.isGroup(pattID)) {
      locXYZ = dict.getClusterCoordinates(c);
    } else {
      o2::itsmft::ClusterPattern cPattern(pattIt);
      locXYZ = dict.getClusterCoordinates(c, cPattern);
    }
  } else {
    o2::itsmft::ClusterPattern cPattern(pattIt);
    locXYZ = dict.getClusterCoordinates(c, cPattern, false);
  }
  setClusterPosition(locXYZ, isLocal);
  // Transformation local --> global coordinates
  auto gloXYZ = geom->getMatrixL2G(chipID) * locXYZ;
  setSensor(c.getSensorID(), chipMappingMFT);
  setClusterPosition(gloXYZ, isGlobal);
  double cook = 1.0; // WARNING!! COOKED
  setMeasuredErrors(cook * sigmaX2, cook * sigmaY2, 0.);
}

//_________________________________________________________
void Hit::convertT2LTrackPosition(UShort_t chipID,
                                  o2::mft::GeometryTGeo* geom)
{
  // Transformation tracking (global) --> local coordinates
  auto locXYZ = geom->getMatrixT2L(chipID) * getTrackXYZ(isGlobal);
  setTrackPosition(locXYZ, isLocal);
}