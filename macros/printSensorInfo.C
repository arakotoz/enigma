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

#include "helperGeom.h"

#endif

// root -l
// .L ~/cernbox/alice/enigma/macros/printSensorInfo.C++
// printSensorInfo()

void printSensorInfo(
  const bool wAllSensors = true,
  const bool wSymName = true,
  const bool wTranslation = true,
  const bool wRotation = true,
  const bool wDeg = true,
  const bool preferAlignedFile = false)
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
    printSensorGlobalTransform(chipHelper, wSymName, wTranslation, wRotation, wDeg);
  }
}
