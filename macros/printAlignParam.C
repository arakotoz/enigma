#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <iomanip>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <TFile.h>
#include <TObject.h>
#include <TString.h>
#include <Rtypes.h>

#include "DetectorsCommonDataFormats/AlignParam.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"

#include "MFTAlignment/AlignSensorHelper.h"

#include "helperGeom.h"

#endif

// root -l
// .L ~/cernbox/alice/enigma/macros/printAlignParam.C++
// printAlignParam()

void printAlignParam(std::string generalPath = "/Users/andry/cernbox/alice/enigma/common-input/lhc2022h/", //"/Users/andry/cernbox/alice/mft/pilotbeam/505713/out-mille/pass1",
                     std::string alignParamFileName = "pass2_mft_alignment_half_mft_shifted",              //"pass2_mft_alignment_sensor_shifted", // "pass2_mft_alignment",
                     const bool wTranslation = true,
                     const bool wRotation = true,
                     const bool wDeg = false)
{
  // geometry

  const bool preferAlignedFile = false;
  const bool applyMisalignment = false;
  o2::base::GeometryManager::loadGeometry("", applyMisalignment, preferAlignedFile);
  o2::mft::GeometryTGeo* geom = o2::mft::GeometryTGeo::Instance();
  geom->fillMatrixCache(
    o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                             o2::math_utils::TransformType::L2G));

  // load alignement parameters from alignParamFileName
  std::vector<o2::detectors::AlignParam> alignParameters;
  try {
    alignParameters = loadAlignParam(
      Form("%s/%s", generalPath.c_str(), alignParamFileName.c_str()));
  } catch (std::exception e) {
    LOG(fatal) << "Abort, " << e.what();
    return;
  }

  // print alignment parameters

  const bool printScreen = true;
  printAlignParam(
    Form("%s/%s", generalPath.c_str(), alignParamFileName.c_str()),
    alignParameters,
    printScreen, wTranslation, wRotation, wDeg);
}
