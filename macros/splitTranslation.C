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

#include "CommonUtils/NameConf.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"

#include "MFTAlignment/AlignSensorHelper.h"

#include "helperGeom.h"

#endif

// root -l
// .L ~/cernbox/alice/enigma/macros/splitTranslation.C++
// splitTranslation()

void splitTranslation(std::string generalPath = "/Users/andry/cernbox/alice/enigma/common-input/lhc2022h/", //"/Users/andry/cernbox/alice/mft/pilotbeam/505713/out-mille/pass1",
                      std::string alignParamFileName = "pass2_mft_alignment")
{
  // geometry

  const bool preferAlignedFile = true;
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

  // create vectors to store align params

  std::vector<o2::detectors::AlignParam> alignParamsXY;
  std::vector<o2::detectors::AlignParam> alignParamsZ;

  // fill the two vectors of align param using loaded align params:
  // - one vector (XY) will have dz = 0
  // - one vector (Z) will have dx = dy = 0

  o2::itsmft::ChipMappingMFT chipMappingMFT;
  int NChips = o2::itsmft::ChipMappingMFT::NChips;

  o2::mft::AlignSensorHelper chipHelper;

  double dx = 0., dy = 0., dz = 0., dRx = 0., dRy = 0., dRz = 0.;

  for (int iChip = 0; iChip < NChips; iChip++) {

    chipHelper.setSensorOnlyInfo(iChip);
    dx = alignParameters[iChip].getX();
    dy = alignParameters[iChip].getY();
    dz = alignParameters[iChip].getZ();
    dRx = alignParameters[iChip].getPsi();
    dRy = alignParameters[iChip].getTheta();
    dRz = alignParameters[iChip].getPhi();

    const bool global = true;
    const double dx0 = 0., dy0 = 0., dz0 = 0.;

    // align param vector (XY) have dz = 0

    alignParamsXY.emplace_back(
      chipHelper.geoSymbolicName(),
      chipHelper.sensorUid(),
      dx, dy, dz0,
      dRx, dRy, dRz,
      global);

    // align param vector (Z) have dx = dy = 0

    alignParamsZ.emplace_back(
      chipHelper.geoSymbolicName(),
      chipHelper.sensorUid(),
      dx0, dy0, dz,
      dRx, dRy, dRz,
      global);
  }

  // save alignment parameters to files

  LOGF(info, "Storing align param vector (XY) in local file %s/%s-xy.root",
       generalPath.c_str(), alignParamFileName.c_str());
  TFile outAlignParamFileXY(
    Form("%s/%s-xy.root", generalPath.c_str(), alignParamFileName.c_str()),
    "recreate", "", 505);
  outAlignParamFileXY.WriteObjectAny(
    &alignParamsXY, "std::vector<o2::detectors::AlignParam>",
    o2::base::NameConf::CCDBOBJECT.data());
  outAlignParamFileXY.Close();

  LOGF(info, "Storing align param vector (Z) in local file %s/%s-z.root",
       generalPath.c_str(), alignParamFileName.c_str());
  TFile outAlignParamFileZ(
    Form("%s/%s-z.root", generalPath.c_str(), alignParamFileName.c_str()),
    "recreate", "", 505);
  outAlignParamFileZ.WriteObjectAny(
    &alignParamsZ, "std::vector<o2::detectors::AlignParam>",
    o2::base::NameConf::CCDBOBJECT.data());
  outAlignParamFileZ.Close();
}