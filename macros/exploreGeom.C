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

#endif

// root -l
// .L ~/cernbox/alice/enigma/macros/exploreGeom.C++
// exploreGeom()

std::vector<o2::detectors::AlignParam> loadAlignParam(std::string alignParamFileName);

void printAlignParam(std::string alignParamFileName,
                     std::vector<o2::detectors::AlignParam> alignParameters,
                     bool printScreen = false,
                     bool wTranslation = true,
                     bool wRotation = true,
                     bool wDeg = false);

void printSensorGlobalTransform(o2::mft::AlignSensorHelper chipHelper,
                                bool wSymName = true,
                                bool wTranslation = true,
                                bool wRotation = true,
                                bool wDeg = true);

//_________________________________
void exploreGeom(std::string alignParamFileName1 = "mftprealignment",
                 std::string alignParamFileName2 = "pass1_mft_alignment",
                 const bool wCombinedAlignParams = true,
                 const bool wAllSensors = true)
{

  // load ideal geometry (o2sim_geometry.root)

  const bool applyMisalignment = false;
  const bool preferAlignedFile = false;
  o2::base::GeometryManager::loadGeometry("", applyMisalignment, preferAlignedFile);
  o2::mft::GeometryTGeo* geom = o2::mft::GeometryTGeo::Instance();

  // load alignement parameters from first file

  auto alignParameters1 = loadAlignParam(alignParamFileName1);
  printAlignParam(alignParamFileName1, alignParameters1);

  // load alignement parameters from second file

  auto alignParameters2 = loadAlignParam(alignParamFileName2);
  printAlignParam(alignParamFileName2, alignParameters2);

  if (wCombinedAlignParams) { // combine alignment parameters

    std::vector<o2::detectors::AlignParam> outAlignParams;

    o2::mft::AlignSensorHelper chipHelper;
    double dRx = 0., dRy = 0., dRz = 0.; // delta rotations
    double dx = 0., dy = 0., dz = 0.;    // delta translations
    bool global = true;                  // delta in global ref. system

    // mapping

    o2::itsmft::ChipMappingMFT chipMappingMFT;
    const int numberOfSensors = o2::itsmft::ChipMappingMFT::NChips;

    // add alignment parameters (global delta)

    for (int chipId = 0; chipId < numberOfSensors; chipId++) {
      chipHelper.setSensorOnlyInfo(chipId);
      dx = alignParameters1[chipId].getX() + alignParameters2[chipId].getX();
      dy = alignParameters1[chipId].getY() + alignParameters2[chipId].getY();
      dz = alignParameters1[chipId].getZ() + alignParameters2[chipId].getZ();
      dRx = alignParameters1[chipId].getPsi() + alignParameters2[chipId].getPsi();
      dRy = alignParameters1[chipId].getTheta() + alignParameters2[chipId].getTheta();
      dRz = alignParameters1[chipId].getPhi() + alignParameters2[chipId].getPhi();

      outAlignParams.emplace_back(
        chipHelper.geoSymbolicName(),
        chipHelper.sensorUid(),
        dx, dy, dz,
        dRx, dRy, dRz,
        global);
    }

    // print and save combined alignment parameters

    std::string OutAlignParamFileName = "combined_prealigned_and_pass1_mft_alignment";
    printAlignParam(OutAlignParamFileName, outAlignParams);
    LOGF(info, "Storing *combined* MFT alignment params in local file %s.root",
         OutAlignParamFileName.c_str());
    TFile outAlignParamFile(Form("%s.root", OutAlignParamFileName.c_str()),
                            "recreate", "", 505);
    outAlignParamFile.WriteObjectAny(&outAlignParams, "std::vector<o2::detectors::AlignParam>", "alignment");
    outAlignParamFile.Close();

    // apply the new set of alignment parameters

    bool isAlignApplied = o2::base::GeometryManager::applyAlignment(outAlignParams);
    if (isAlignApplied) {
      LOG(info) << "Successfully applied *combined* alignment parameters from "
                << alignParamFileName1 << ".root + "
                << alignParamFileName2 << ".root";
    }

  } else {

    // apply alignment parameters from first file

    bool isAlignApplied1 = o2::base::GeometryManager::applyAlignment(alignParameters1);
    if (isAlignApplied1) {
      LOG(info) << "Successfully applied alignment parameters from "
                << alignParamFileName1 << ".root";
    }

    // apply alignment parameters from second file

    bool isAlignApplied2 = o2::base::GeometryManager::applyAlignment(alignParameters1);
    if (isAlignApplied2) {
      LOG(info) << "Successfully applied alignment parameters from "
                << alignParamFileName2 << ".root";
    }
  }

  // generate aligned geometry file (o2sim_geometry-aligned.root)

  auto alignedgeomfile = o2::base::NameConf::getAlignedGeomFileName();
  gGeoManager->Export(Form("new-%s", alignedgeomfile.c_str()));
  LOG(info) << "New geometry file generated : "
            << Form("new-%s", alignedgeomfile.c_str());

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
                     bool printScreen,
                     bool wTranslation,
                     bool wRotation,
                     bool wDeg)
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