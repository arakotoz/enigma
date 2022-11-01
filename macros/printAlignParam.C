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

#endif

// root -l
// .L ~/cernbox/alice/enigma/macros/printAlignParam.C++
// printAlignParam()

void printAlignParam(std::string generalPath = "/Users/andry/cernbox/alice/mft/pilotbeam/505713/out-mille/pass1",
                     std::string alignParamFileName = "mft_alignment",
                     bool wTranslation = true,
                     bool wRotation = true,
                     bool wDeg = false)
{
  // geometry

  const bool preferAlignedFile = true;
  const bool applyMisalignment = false;
  o2::base::GeometryManager::loadGeometry("", applyMisalignment, preferAlignedFile);
  o2::mft::GeometryTGeo* geom = o2::mft::GeometryTGeo::Instance();
  geom->fillMatrixCache(
    o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                             o2::math_utils::TransformType::L2G));

  // load alignement parameters (mft_alignment.root)

  std::vector<o2::detectors::AlignParam> alignParameters;
  if (!alignParamFileName.empty()) {
    std::cout << "Loading alignment parameters from "
              << Form("%s/%s.root", generalPath.c_str(), alignParamFileName.c_str())
              << std::endl;
    TFile algFile(Form("%s/%s.root", generalPath.c_str(), alignParamFileName.c_str()));
    auto alignement = algFile.Get<std::vector<o2::detectors::AlignParam>>("alignment");
    if (alignement) {
      alignParameters = *alignement;
    } else {
      throw std::exception();
    }
    algFile.Close();
  }

  // print alignment parameters

  o2::itsmft::ChipMappingMFT chipMappingMFT;
  int NChips = o2::itsmft::ChipMappingMFT::NChips;

  std::ofstream OutStream;
  OutStream.open(Form("%s/%s.csv", generalPath.c_str(), alignParamFileName.c_str()));
  OutStream << "half,disk,layer,zone,con,tr,chipid,dx,dy,dz,dRx,dRy,dRz" << endl;

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

    OutStream << chipHelper.half() << ","
              << chipHelper.disk() << ","
              << chipHelper.layer() << ","
              << chipHelper.zone() << ","
              << chipHelper.connector() << ","
              << chipHelper.transceiver() << ","
              << iChip << ","
              << dx << "," << dy << "," << dz << ","
              << dRx << "," << dRy << "," << dRz
              << endl;

    std::streamsize ss = std::cout.precision();
    bool wSymName = true;
    std::stringstream name = chipHelper.getSensorFullName(wSymName);
    std::cout << name.str().c_str()
              << " (" << alignParameters[iChip].getSymName() << ") ";
    if (wTranslation) {
      std::cout << std::scientific << std::setprecision(2)
                << " (cm) dx " << dx
                << " dy " << dy
                << " dz " << dz;
    }
    if (wRotation) {
      constexpr double rad2deg = 180.0 / 3.14159265358979323846;
      if (wDeg) {
        dRx *= rad2deg;
        dRy *= rad2deg;
        dRz *= rad2deg;
        std::cout << " (deg)";
      }
      std::cout << std::scientific << std::setprecision(2)
                << " dRx " << dRx << " dRy " << dRy << " dRz " << dRz;
    }
    std::cout << std::setprecision(ss) << std::endl;
  }
  OutStream.close();
}