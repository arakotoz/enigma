#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <TFile.h>
#include <TString.h>
#include <Rtypes.h>

#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"
#include "MFTAlignment/AlignSensorHelper.h"

#endif

// root -l
// .L ~/cernbox/alice/enigma/macros/printSensorInfo.C++
// printSensorInfo()

void print(o2::mft::AlignSensorHelper chipHelper,
           bool wSymName,
           bool wTranslation,
           bool wRotation,
           bool wDeg);

void printSensorInfo(
  const bool wAllSensors = true,
  const bool wSymName = true,
  const bool wTranslation = true,
  const bool wRotation = true,
  const bool wDeg = true,
  const bool preferAlignedFile = true)
{
  // geometry

  const bool applyMisalignment = false;
  o2::base::GeometryManager::loadGeometry("", applyMisalignment, preferAlignedFile);
  o2::mft::GeometryTGeo* geom = o2::mft::GeometryTGeo::Instance();
  geom->fillMatrixCache(
    o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                             o2::math_utils::TransformType::L2G));

  // mapping

  o2::itsmft::ChipMappingMFT chipMappingMFT;

  int NChips = o2::itsmft::ChipMappingMFT::NChips;
  if (!wAllSensors) {
    NChips = 10;
  }

  // print info for a given sensor

  o2::mft::AlignSensorHelper chipHelper;
  for (int iChip = 0; iChip < NChips; iChip++) {
    chipHelper.setSensor(iChip);
    print(chipHelper, wSymName, wTranslation, wRotation, wDeg);
  }
}

void print(o2::mft::AlignSensorHelper chipHelper,
           bool wSymName,
           bool wTranslation,
           bool wRotation,
           bool wDeg)
{
  std::streamsize ss = std::cout.precision();
  std::stringstream name = chipHelper.getSensorFullName(wSymName);
  std::cout << name.str().c_str();
  if (wTranslation) {
    std::cout << std::scientific << std::setprecision(2)
              << " (cm) dx " << chipHelper.translateX()
              << " dy " << chipHelper.translateY()
              << " dz " << chipHelper.translateZ();
  }
  if (wRotation) {
    constexpr double rad2deg = 180.0 / 3.14159265358979323846;
    double rotX = chipHelper.angleRx();
    double rotY = chipHelper.angleRy();
    double rotZ = chipHelper.angleRz();
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