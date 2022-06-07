// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// @author  Laurent Aphecetche

#define BOOST_TEST_MODULE Test MCHSimulation GeometryMisAligner
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "MCHGeometryCreator/Geometry.h"
#include "MCHGeometryTest/Helpers.h"
#include "MCHGeometryMisAligner/MisAligner.h"
#include "MCHGeometryTransformer/Transformations.h"
#include "TGeoManager.h"
#include "MathUtils/Cartesian.h"
#include "Math/GenVector/Cartesian3D.h"
#include "boost/format.hpp"
#include <boost/test/data/test_case.hpp>
#include <iomanip>
#include <iostream>
#include <fmt/format.h>

namespace but = boost::unit_test;
namespace bdata = boost::unit_test::data;
namespace btools = boost::test_tools;

BOOST_TEST_DONT_PRINT_LOG_VALUE(o2::mch::geo::MisAligner)

struct GEOMETRY {
  GEOMETRY()
  {
    if (!gGeoManager) {
      o2::mch::test::createStandaloneGeometry();
      o2::mch::geo::addAlignableVolumes(*gGeoManager);
    }
  };
}; // namespace boost::test_toolsBOOST_TEST_DONT_PRINT_LOG_VALUE(o2::mch::geo::MisAligner)structGEOMETRY

const std::vector<std::vector<std::string>> deSymNames{
  {"DE100", "DE103"},
  {"DE101", "DE102"},
  {"DE200", "DE203"},
  {"DE201", "DE202"},
  {"DE300", "DE303"},
  {"DE301", "DE302"},
  {"DE400", "DE403"},
  {"DE401", "DE402"},
  {"DE500", "DE501", "DE502", "DE503", "DE504", "DE514", "DE515", "DE516", "DE517"},
  {"DE505", "DE506", "DE507", "DE508", "DE509", "DE510", "DE511", "DE512", "DE513"},
  {"DE600", "DE601", "DE602", "DE603", "DE604", "DE614", "DE615", "DE616", "DE617"},
  {"DE605", "DE606", "DE607", "DE608", "DE609", "DE610", "DE611", "DE612", "DE613"},
  {"DE700", "DE701", "DE702", "DE703", "DE704", "DE705", "DE706", "DE720", "DE721", "DE722", "DE723", "DE724", "DE725"},
  {"DE707", "DE708", "DE709", "DE710", "DE711", "DE712", "DE713", "DE714", "DE715", "DE716", "DE717", "DE718", "DE719"},
  {"DE800", "DE801", "DE802", "DE803", "DE804", "DE805", "DE806", "DE820", "DE821", "DE822", "DE823", "DE824", "DE825"},
  {"DE807", "DE808", "DE809", "DE810", "DE811", "DE812", "DE813", "DE814", "DE815", "DE816", "DE817", "DE818", "DE819"},
  {"DE900", "DE901", "DE902", "DE903", "DE904", "DE905", "DE906", "DE920", "DE921", "DE922", "DE923", "DE924", "DE925"},
  {"DE907", "DE908", "DE909", "DE910", "DE911", "DE912", "DE913", "DE914", "DE915", "DE916", "DE917", "DE918", "DE919"},
  {"DE1000", "DE1001", "DE1002", "DE1003", "DE1004", "DE1005", "DE1006", "DE1020", "DE1021", "DE1022", "DE1023", "DE1024", "DE1025"},
  {"DE1007", "DE1008", "DE1009", "DE1010", "DE1011", "DE1012", "DE1013", "DE1014", "DE1015", "DE1016", "DE1017", "DE1018", "DE1019"}};

BOOST_FIXTURE_TEST_SUITE(geometrycreator, GEOMETRY)

BOOST_AUTO_TEST_CASE(ZeroMisAlignHalfChambers, *but::tolerance(0.00001))
{
  BOOST_REQUIRE(gGeoManager != nullptr);

  std::vector<o2::detectors::AlignParam> params;

  // The misaligner
  o2::mch::geo::MisAligner aGMA;

  aGMA.misAlign(params);

  int nvols = params.size();
  for (int i = 0; i < nvols; i++) {
    BOOST_TEST(params[i].getX() == 0.0);
    BOOST_TEST(params[i].getY() == 0.0);
    BOOST_TEST(params[i].getZ() == 0.0);
    BOOST_TEST(params[i].getPsi() == 0.0);
    BOOST_TEST(params[i].getTheta() == 0.0);
    BOOST_TEST(params[i].getPhi() == 0.0);
  }
}

BOOST_AUTO_TEST_CASE(MisAlignHalfChambers, *but::tolerance(0.00001))
{
  BOOST_REQUIRE(gGeoManager != nullptr);

  std::vector<o2::detectors::AlignParam> params;

  // The misaligner
  o2::mch::geo::MisAligner aGMA;

  // To generate module mislaignment (not mandatory)
  aGMA.setModuleCartMisAlig(0.1, 0.0, 0.2, 0.0, 0.3, 0.0);
  aGMA.setModuleAngMisAlig(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

  aGMA.misAlign(params);

  int nvols = params.size();
  for (int i = 0; i < nvols; i++) {
    if (params[i].getSymName() == fmt::format("MCH/HC{}", 0)) {
      BOOST_TEST(params[i].getX() == 0.1);
      BOOST_TEST(params[i].getY() == 0.2);
      BOOST_TEST(params[i].getZ() == 0.3);
      BOOST_TEST(params[i].getPsi() == 0.0);
      BOOST_TEST(params[i].getTheta() == 0.0);
      BOOST_TEST(params[i].getPhi() == 0.0);
      break;
    }
  }
}

BOOST_AUTO_TEST_CASE(MisAlignDetectionElements, *but::tolerance(0.00001))
{
  BOOST_REQUIRE(gGeoManager != nullptr);

  std::vector<o2::detectors::AlignParam> params;

  // The misaligner
  o2::mch::geo::MisAligner aGMA;

  // To generate detection element misalignment
  aGMA.setCartMisAlig(0.01, 0.0, 0.02, 0.0, 0.03, 0.0);
  aGMA.setAngMisAlig(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

  aGMA.misAlign(params);

  int nvols = params.size();
  for (int i = 0; i < nvols; i++) {
    std::cout << params[i].getSymName() << std::endl;
    if (params[i].getSymName() == fmt::format("MCH/HC{}/DE{}", 0, 100)) {
      BOOST_TEST(params[i].getX() == 0.01);
      BOOST_TEST(params[i].getY() == 0.02);
      BOOST_TEST(params[i].getZ() == 0.03);
      BOOST_TEST(params[i].getPsi() == 0.0);
      BOOST_TEST(params[i].getTheta() == 0.0);
      BOOST_TEST(params[i].getPhi() == 0.0);
      break;
    }
  }
  // for (int hc = 0; hc < 20; hc++) {
  //   for (int de = 0; de < deSymNames[hc].size(); de++) {
  //     BOOST_CHECK((gGeoManager->GetAlignableEntry((fmt::format("MCH/HC{}/{}", hc, deSymNames[hc][de].c_str())).c_str())));
  //   }
  // }
}

BOOST_AUTO_TEST_CASE(MisAlignHCDE, *but::tolerance(0.00001))
{
  BOOST_REQUIRE(gGeoManager != nullptr);

  std::vector<o2::detectors::AlignParam> params;

  // The misaligner
  o2::mch::geo::MisAligner aGMA;

  auto transformationB = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);
  auto tB = transformationB(100);
  o2::math_utils::Point3D<double> poB;
  tB.LocalToMaster(o2::math_utils::Point3D<double>{0, 0, 0}, poB);

  // To generate half chmaber misalignment
  aGMA.setModuleCartMisAlig(0.1, 0.0, 0.2, 0.0, 0.3, 0.0);
  aGMA.setModuleAngMisAlig(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  // To generate detection element misalignment
  aGMA.setCartMisAlig(0.01, 0.0, 0.02, 0.0, 0.03, 0.0);
  aGMA.setAngMisAlig(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

  aGMA.misAlign(params);

  auto transformation = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);
  auto t = transformation(100);
  o2::math_utils::Point3D<double> po;
  t.LocalToMaster(o2::math_utils::Point3D<double>{0, 0, 0}, po);

  BOOST_TEST(po.X() - poB.X() == 0.11);
  tB = transformationB(101);
  t = transformation(101);
  tB.LocalToMaster(o2::math_utils::Point3D<double>{0, 0, 0}, poB);
  t.LocalToMaster(o2::math_utils::Point3D<double>{0, 0, 0}, po);
  BOOST_TEST(po.X() - poB.X() == 0.09);
  std::cout << po.X() << " " << poB.X() << std::endl;
  int nvols = params.size();
  // for (int i = 0; i < nvols; i++) {
  //   std::cout << params[i].getSymName() << std::endl;
  //   if (params[i].getSymName() == fmt::format("MCH/HC{}/DE{}", 0, 100)) {
  //     BOOST_TEST(params[i].getX() == 0.11);
  //     // BOOST_TEST(params[i].getY() == 0.02);
  //     // BOOST_TEST(params[i].getZ() == 0.03);
  //     // BOOST_TEST(params[i].getPsi() == 0.0);
  //     // BOOST_TEST(params[i].getTheta() == 0.0);
  //     // BOOST_TEST(params[i].getPhi() == 0.0);
  //     break;
  //   }
  // }
  // for (int hc = 0; hc < 20; hc++) {
  //   for (int de = 0; de < deSymNames[hc].size(); de++) {
  //     BOOST_CHECK((gGeoManager->GetAlignableEntry((fmt::format("MCH/HC{}/{}", hc, deSymNames[hc][de].c_str())).c_str())));
  //   }
  // }
}

BOOST_AUTO_TEST_SUITE_END()
