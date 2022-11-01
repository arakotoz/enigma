#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <iomanip>
#include <iostream>
#include <string>

#include <TFile.h>
#include <TObject.h>
#include <Rtypes.h>

#include "CommonUtils/NameConf.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"

#endif

// root -l
// .L ~/cernbox/alice/enigma/macros/exportAlignedGeom.C++
// exportAlignedGeom()

void exportAlignedGeom(std::string alignParamFileName = "mft_alignment.root")
{

  // load ideal geometry (o2sim_geometry.root)

  const bool applyMisalignment = false;
  const bool preferAlignedFile = false;
  o2::base::GeometryManager::loadGeometry("", applyMisalignment, preferAlignedFile);
  o2::mft::GeometryTGeo* geom = o2::mft::GeometryTGeo::Instance();

  // load pre-alignement parameters (mftprealignment.root)

  std::vector<o2::detectors::AlignParam> alignParameters;
  if (!alignParamFileName.empty()) {
    std::cout << "Loading alignment parameters from " << alignParamFileName << std::endl;
    TFile algFile(alignParamFileName.c_str());
    auto alignement = algFile.Get<std::vector<o2::detectors::AlignParam>>("alignment");
    if (alignement) {
      alignParameters = *alignement;
    } else {
      throw std::exception();
    }
    algFile.Close();
  }

  // apply alignment

  bool isAlignApplied = o2::base::GeometryManager::applyAlignment(alignParameters);
  if (isAlignApplied) {
    std::cout << "Successfully applied alignment parameters from " << alignParamFileName << std::endl;
  }

  // generate aligned geometry file (o2sim_geometry-aligned.root)

  auto alignedgeomfile = o2::base::NameConf::getAlignedGeomFileName();
  gGeoManager->Export(alignedgeomfile.c_str());
}