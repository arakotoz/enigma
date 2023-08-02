#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <iomanip>
#include <iostream>
#include <string>

#include <TFile.h>
#include <TObject.h>
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
// .L ~/cernbox/alice/enigma/macros/exportAlignedGeom.C++
// exportAlignedGeom()

void exportAlignedGeom(
  std::string generalPath = "/Users/andry/cernbox/alice/enigma/common-input/lhc2022h/",
  std::string alignParamFileName = "mft_survey_disk", // "pass2_mft_alignment", "pass1_wrt_ideal_mft_alignment", //"pass1_mft_alignment.root", //"mft_alignment.root"
  const bool wAllSensors = true)
{

  // load ideal geometry (o2sim_geometry.root)

  const bool applyMisalignment = false;
  const bool preferAlignedFile = false;
  o2::base::GeometryManager::loadGeometry("", applyMisalignment, preferAlignedFile);
  o2::mft::GeometryTGeo* geom = o2::mft::GeometryTGeo::Instance();

  // load alignement parameters

  std::vector<o2::detectors::AlignParam> alignParameters;
  try {
    alignParameters = loadAlignParam(alignParamFileName);
  } catch (std::exception e) {
    LOG(fatal) << "Abort, " << e.what();
    return;
  }

  // apply alignment

  bool isAlignApplied = o2::base::GeometryManager::applyAlignment(alignParameters);
  if (isAlignApplied) {
    LOG(info) << "Successfully applied alignment parameters from " << alignParamFileName;

    // generate aligned geometry file (o2sim_geometry-aligned.root)
    auto alignedgeomfile = o2::base::NameConf::getAlignedGeomFileName();
    gGeoManager->Export(Form("new-%s", alignedgeomfile.c_str()));
    LOG(info) << "New geometry file generated : "
              << Form("new-%s", alignedgeomfile.c_str());
  }

  // options to be used when printing the global transformations in the new geometry

  const bool wSymName = true;
  const bool wTranslation = true;
  const bool wRotation = true;
  const bool wDeg = true;

  // cache matrix elements in the geometry

  geom->fillMatrixCache(
    o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                             o2::math_utils::TransformType::L2G));

  // mapping

  o2::itsmft::ChipMappingMFT chipMappingMFT;

  int NChips = o2::itsmft::ChipMappingMFT::NChips;
  if (!wAllSensors) {
    NChips = 10;
  }

  // print global transformations in the new geometry for a given sensor

  o2::mft::AlignSensorHelper chipHelper;
  for (int iChip = 0; iChip < NChips; iChip++) {
    chipHelper.setSensor(iChip);
    printSensorGlobalTransform(chipHelper, wSymName, wTranslation, wRotation, wDeg);
  }
}
