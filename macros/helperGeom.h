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

#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"

#include "MFTAlignment/AlignSensorHelper.h"

//_________________________________
std::vector<o2::detectors::AlignParam> loadAlignParam(std::string alignParamFileName)
{
  if (alignParamFileName.empty()) {
    LOG(fatal) << "No input align params file name provided !";
    throw std::exception();
  }
  LOG(info) << "Loading alignment parameters from " << alignParamFileName << ".root";
  TFile algFile(Form("%s.root", alignParamFileName.c_str()));
  if (algFile.IsZombie()) {
    LOG(fatal) << "Bad align param file " << alignParamFileName << ".root";
    throw std::exception();
  }
  auto alignment = algFile.Get<std::vector<o2::detectors::AlignParam>>("alignment");
  algFile.Close();
  if (!alignment) {
    LOG(fatal) << "Empty vector of align params !";
    throw std::exception();
  }
  std::vector<o2::detectors::AlignParam> alignParameters = *alignment;
  return alignParameters;
}

//_________________________________
void printAlignParam(std::string alignParamFileName,
                     std::vector<o2::detectors::AlignParam> alignParameters,
                     bool printScreen = false,
                     bool wTranslation = true,
                     bool wRotation = true,
                     bool wDeg = false)
{
  o2::itsmft::ChipMappingMFT chipMappingMFT;
  int NChips = o2::itsmft::ChipMappingMFT::NChips;

  std::ofstream OutStream;
  OutStream.open(Form("%s.csv", alignParamFileName.c_str()));
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

    if (printScreen) {
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
  }
  OutStream.close();
}

//_________________________________
void printSensorGlobalTransform(o2::mft::AlignSensorHelper chipHelper,
                                bool wSymName = true,
                                bool wTranslation = true,
                                bool wRotation = true,
                                bool wDeg = true)
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