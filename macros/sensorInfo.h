#include <iomanip>
#include <iostream>
#include <string>

#include <TGeoMatrix.h>
#include <TObject.h>
#include <TString.h>
#include <Rtypes.h>

#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "MathUtils/Cartesian.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"
#include "MFTTracking/IOUtils.h"
#include "Framework/Logger.h"

class SensorInfo
{
 public:
  SensorInfo() = delete;
  SensorInfo(const int chipIndex, const o2::mft::GeometryTGeo* geom);
  ~SensorInfo() = default;

  void setSensor(const int chipIndex);

  /// \brief set pointer to geometry that should already have done fillMatrixCache()
  void setGeometry(const o2::mft::GeometryTGeo* geom);

  /// \brief set the ALICE global unique id of the sensor
  void setSensorUid(const int chipIndex);

  /// \brief print sensor info to screen
  void print(bool wSymName = true,
             bool wTranslation = true,
             bool wRotation = true);

 protected:
  Int_t mChipIndexOnLadder = 0;                     ///< sensor index within the ladder [0, 4]
  Int_t mChipIndexInMft = 0;                        ///< sensor sw index within the MFT [0, 935]
  Int_t mLadderInHalfDisk = 0;                      ///< ladder geo index in this half MFT disk [0, 33]
  Int_t mConnector = 0;                             ///< connector index to which the ladder is plugged in the zone [0, 4]
  Int_t mTransceiver = 0;                           ///< transceiver id to which the sensor is connected in the zone [0, 24]
  Int_t mLayer = 0;                                 ///< layer id [0, 9]
  Int_t mZone = 0;                                  ///< zone id [0,3]
  Int_t mDisk = 0;                                  ///< disk id [0, 4]
  Int_t mHalf = 0;                                  ///< half id [0, 1]
  Int_t mChipUniqueId = 0;                          ///< ALICE global unique id of the sensor
  TGeoHMatrix mTransform;                           ///< sensor transformation matrix L2G
  static o2::itsmft::ChipMappingMFT mChipMapping;   ///< MFT chip <-> ladder, layer, disk, half mapping
  const o2::mft::GeometryTGeo* mGeometry = nullptr; ///< MFT geometry
  TString mGeoSymbolicName;

  /// \brief set the symbolic name of this sensor in the geometry
  void setSymName();

  /// \brief set the matrix that stores the sensor transform L2G
  void setSensorTransform();

  ClassDefNV(SensorInfo, 1);
};

//__________________________________________________________________________
SensorInfo::SensorInfo(const int chipIndex, const o2::mft::GeometryTGeo* geom)
  : mChipIndexOnLadder(0),
    mChipIndexInMft(0),
    mLadderInHalfDisk(0),
    mConnector(0),
    mTransceiver(0),
    mLayer(0),
    mZone(0),
    mDisk(0),
    mHalf(0),
    mChipUniqueId(0)
{
  setGeometry(geom);
  setSensor(chipIndex);
  setSensorUid(chipIndex);
  setSymName();
}

//__________________________________________________________________________
void SensorInfo::setGeometry(const o2::mft::GeometryTGeo* geom)
{
  if (mGeometry == nullptr) {
    mGeometry = geom;
    mGeoSymbolicName = mGeometry->composeSymNameMFT();
  }
}

//__________________________________________________________________________
void SensorInfo::setSensor(const int chipIndex)
{
  if (chipIndex < mChipMapping.getNChips()) {
    o2::itsmft::MFTChipMappingData chipMappingData = (mChipMapping.getChipMappingData())[chipIndex];
    mChipIndexOnLadder = (UShort_t)chipMappingData.chipOnModule;
    mChipIndexInMft = chipMappingData.globalChipSWID;
    mLayer = (UShort_t)chipMappingData.layer;
    mConnector = (UShort_t)chipMappingData.connector;
    mTransceiver = (UShort_t)chipMappingData.cable;
    mZone = (UShort_t)chipMappingData.zone;
    mDisk = (UShort_t)chipMappingData.disk;
    mHalf = (UShort_t)chipMappingData.half;
  } else {
    LOG(error) << "SensorInfo::setSensor() - "
               << "chip index " << chipIndex
               << " >= " << mChipMapping.getNChips() << std::endl;
  }
}

//__________________________________________________________________________
void SensorInfo::setSensorUid(const int chipIndex)
{
  if (chipIndex < mChipMapping.getNChips()) {
    mChipUniqueId = o2::base::GeometryManager::getSensID(o2::detectors::DetID::MFT,
                                                         chipIndex);
  } else {
    LOG(error) << "SensorInfo::setSensorUid() - "
               << "chip index " << chipIndex
               << " >= " << mChipMapping.getNChips() << std::endl;
    mChipUniqueId = o2::base::GeometryManager::getSensID(o2::detectors::DetID::MFT, 0);
  }
}

//__________________________________________________________________________
void SensorInfo::setSymName()
{
  int hf = 0, dk = 0, sr = 0;
  if (mGeometry) {
    mGeometry->getSensorID(mChipIndexInMft, hf, dk, mLadderInHalfDisk, sr);
    bool isIdVerified = true;
    isIdVerified &= (hf == (int)mHalf);
    isIdVerified &= (dk == (int)mDisk);
    isIdVerified &= (sr == (int)mChipIndexOnLadder);
    if (isIdVerified) {
      mGeoSymbolicName = mGeometry->composeSymNameChip(mHalf,
                                                       mDisk,
                                                       mLadderInHalfDisk,
                                                       mChipIndexOnLadder);
    } else {
      LOG(error) << "SensorInfo::setSymName() - mismatch in some index"
                 << std::endl;
    }
  } else {
    LOG(error) << "SensorInfo::setSymName() - nullptr to geometry"
               << std::endl;
  }
}

//__________________________________________________________________________
void SensorInfo::setSensorTransform()
{
  if (mGeometry) {
    mTransform = mGeometry->getMatrixL2G(mChipIndexInMft);
  } else {
    LOG(error) << "SensorInfo::setSensorTransform() - nullptr to geometry"
               << std::endl;
  }
}

//__________________________________________________________________________
void SensorInfo::print(bool wSymName,
                       bool wTranslation,
                       bool wRotation)
{
  if (mGeometry == nullptr) {
    wSymName = false;
    wTranslation = false;
    wRotation = false;
  }
  std::streamsize ss = std::cout.precision();
  std::cout << "h " << mHalf << " d " << mDisk << " layer " << mLayer
            << " z " << mZone << " lr " << std::setw(3) << mLadderInHalfDisk
            << " con " << std::setw(1) << mConnector
            << " tr " << std::setw(2) << mTransceiver
            << " sr " << std::setw(1) << mChipIndexOnLadder
            << " iChip " << std::setw(3) << mChipIndexInMft
            << " uid " << mChipUniqueId;
  if (wSymName) {
    std::cout << " " << mGeoSymbolicName;
  }
  if (wTranslation) {
    Double_t* tra = mTransform.GetTranslation();
    std::cout << std::scientific << std::setprecision(2)
              << " tra " << tra[0] << " " << tra[1] << " " << tra[2];
  }
  if (wRotation) {
    constexpr double rad2deg = 180.0 / 3.14159265358979323846;
    Double_t* rot = mTransform.GetRotationMatrix();
    double roll = rad2deg * std::atan2(-rot[5], rot[8]);
    double pitch = rad2deg * std::asin(rot[2]);
    double yaw = rad2deg * std::atan2(-rot[1], rot[0]);
    std::cout << std::scientific << std::setprecision(2)
              << " roll " << roll << " pitch " << pitch << " yaw " << yaw;
  }
  std::cout << std::setprecision(ss) << std::endl;
}