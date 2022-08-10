#include <iomanip>
#include <iostream>
#include <array>
#include <algorithm>

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
  SensorInfo(o2::mft::GeometryTGeo* geom);
  ~SensorInfo() = default;

  void setSensor(const int chipIndex);

  /// \brief print sensor info to screen
  void print(bool wSymName = true,
             bool wTranslation = true,
             bool wRotation = true,
             bool wDeg = true);

 protected:
  UShort_t mChipIndexOnLadder = 0;                                  ///< sensor index within the ladder [0, 4]
  UShort_t mChipIndexInMft = 0;                                     ///< sensor sw index within the MFT [0, 935]
  UShort_t mConnector = 0;                                          ///< connector index to which the ladder is plugged in the zone [0, 4]
  UShort_t mTransceiver = 0;                                        ///< transceiver id to which the sensor is connected in the zone [0, 24]
  UShort_t mLayer = 0;                                              ///< layer id [0, 9]
  UShort_t mZone = 0;                                               ///< zone id [0,3]
  UShort_t mDisk = 0;                                               ///< disk id [0, 4]
  UShort_t mHalf = 0;                                               ///< half id [0, 1]
  Int_t mChipUniqueId = 0;                                          ///< ALICE global unique id of the sensor
  TGeoHMatrix mTransform;                                           ///< sensor transformation matrix L2G
  static o2::itsmft::ChipMappingMFT mChipMapping;                   ///< MFT chip <-> ladder, layer, disk, half mapping
  o2::mft::GeometryTGeo* mGeometry = nullptr;                       ///< MFT geometry
  static constexpr int mNumberOfSensors = mChipMapping.getNChips(); ///< total number of sensors in the MFT
  TString mGeoSymbolicName;                                         ///< symbolic name in the geometry for that sensor
  std::array<TString, mNumberOfSensors> mSymNames;                  ///> array of symbolic names of all sensors

 protected:
  /// \brief set pointer to geometry that should already have done fillMatrixCache()
  void setGeometry(o2::mft::GeometryTGeo* geom);

  /// \brief set the ALICE global unique id of the sensor
  void setSensorUid(const int chipIndex);

  /// \brief build array of symbolic names of all sensors
  void builSymNames();

  /// \brief set the symbolic name of this sensor in the geometry
  void setSymName();

  /// \brief set the matrix that stores the sensor transform L2G
  void setSensorTransform();

  /// \brief return chip id Geo from a given chip id RO
  int getChipIdGeoFromRO(const int chipIdRO);

  ClassDefNV(SensorInfo, 1);
};

//__________________________________________________________________________
SensorInfo::SensorInfo(o2::mft::GeometryTGeo* geom)
  : mChipIndexOnLadder(0),
    mChipIndexInMft(0),
    mConnector(0),
    mTransceiver(0),
    mLayer(0),
    mZone(0),
    mDisk(0),
    mHalf(0),
    mChipUniqueId(0)

{
  setGeometry(geom);
  builSymNames();
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
    setSensorUid(chipIndex);
    setSymName();
    setSensorTransform();
  } else {
    LOG(error) << "SensorInfo::setSensor() - "
               << "chip index " << chipIndex
               << " >= " << mChipMapping.getNChips();
  }
}

//__________________________________________________________________________
void SensorInfo::setGeometry(o2::mft::GeometryTGeo* geom)
{
  if (mGeometry == nullptr) {
    mGeometry = geom;
    mGeoSymbolicName = mGeometry->composeSymNameMFT();
    mSymNames.fill(mGeoSymbolicName);
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
               << " >= " << mChipMapping.getNChips();
    mChipUniqueId = o2::base::GeometryManager::getSensID(o2::detectors::DetID::MFT, 0);
  }
}

//__________________________________________________________________________
void SensorInfo::builSymNames()
{
  if (mGeometry) {
    Int_t iChip = 0;
    Int_t nHalf = mGeometry->getNumberOfHalfs();
    TString sname = mGeometry->composeSymNameMFT();

    for (Int_t hf = 0; hf < nHalf; hf++) {
      Int_t nDisks = mGeometry->getNumberOfDisksPerHalf(hf);
      sname = mGeometry->composeSymNameHalf(hf);

      for (Int_t dk = 0; dk < nDisks; dk++) {
        sname = mGeometry->composeSymNameDisk(hf, dk);

        Int_t nLadders = 0;
        for (Int_t sensor = mGeometry->getMinSensorsPerLadder();
             sensor < mGeometry->getMaxSensorsPerLadder() + 1; sensor++) {
          nLadders += mGeometry->getNumberOfLaddersPerDisk(hf, dk, sensor);
        }

        for (Int_t lr = 0; lr < nLadders; lr++) { // nLadders
          sname = mGeometry->composeSymNameLadder(hf, dk, lr);
          Int_t nSensorsPerLadder = mGeometry->getNumberOfSensorsPerLadder(hf, dk, lr);

          for (Int_t sr = 0; sr < nSensorsPerLadder; sr++) {
            sname = mGeometry->composeSymNameChip(hf, dk, lr, sr);
            mSymNames[iChip] = sname;
            iChip++;
          }
        }
      }
    }
  } else {
    LOG(error) << "SensorInfo::builSymNames() - nullptr to geometry";
  }
}

//__________________________________________________________________________
void SensorInfo::setSymName()
{
  mGeoSymbolicName = mSymNames[getChipIdGeoFromRO(mChipIndexInMft)];
}

//__________________________________________________________________________
void SensorInfo::setSensorTransform()
{
  if (mGeometry) {
    /*
    o2::math_utils::Transform3D matrix = mGeometry->getMatrixL2G(mChipIndexInMft);
    double rot[9], tra[3];
    matrix.GetComponents(rot[0], rot[1], rot[2], tra[0], rot[3], rot[4], rot[5], tra[1], rot[6], rot[7], rot[8], tra[2]);
    mTransform.SetRotation(rot);
    mTransform.SetTranslation(tra);
    std::cout << std::scientific << std::setprecision(2)
              << " tra " << tra[0] << " " << tra[1] << " " << tra[2];
    std::cout << std::scientific << std::setprecision(2)
              << " rot "
              << rot[0] << " " << rot[1] << " " << rot[2] << " "
              << rot[3] << " " << rot[4] << " " << rot[5] << " "
              << rot[6] << " " << rot[7] << " " << rot[8] << std::endl;
    */
    mTransform = mGeometry->getMatrixL2G(mChipIndexInMft);
  } else {
    LOG(error) << "SensorInfo::setSensorTransform() - nullptr to geometry";
  }
}

//__________________________________________________________________________
void SensorInfo::print(bool wSymName,
                       bool wTranslation,
                       bool wRotation,
                       bool wDeg)
{
  if (mGeometry == nullptr) {
    wSymName = false;
    wTranslation = false;
    wRotation = false;
  }
  std::streamsize ss = std::cout.precision();
  std::cout << "h " << mHalf << " d " << mDisk << " layer " << mLayer
            << " z " << mZone << " con " << std::setw(1) << mConnector
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
    double rotX = std::atan2(-rot[5], rot[8]);
    double rotY = std::asin(rot[2]);
    double rotZ = std::atan2(-rot[1], rot[0]);
    if (wDeg) {
      rotX *= rad2deg;
      rotY *= rad2deg;
      rotZ *= rad2deg;
      std::cout << " (deg)";
    }
    std::cout << std::scientific << std::setprecision(2)
              << " Rx " << rotX << " Ry " << rotY << " Rz " << rotZ;
  }
  std::cout << std::setprecision(ss) << std::endl;
}

//__________________________________________________________________________
int SensorInfo::getChipIdGeoFromRO(const int chipIdRO)
{
  int chipIdGeo = 0;
  for (int index = 0; index < mChipMapping.getNChips(); index++) {
    if (o2::itsmft::ChipMappingMFT::mChipIDGeoToRO[index] == chipIdRO) {
      chipIdGeo = index;
      break;
    }
  }
  return chipIdGeo;
}