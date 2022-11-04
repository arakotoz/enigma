// Modified copy of from O2/macro/CreateGRPMagFieldObject.C

#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <ctime>
#include <chrono>
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/CcdbApi.h"
#include "CommonTypes/Units.h"

#endif

using timePoint = long;
using CcdbApi = o2::ccdb::CcdbApi;
using GRPMagField = o2::parameters::GRPMagField;
using current = o2::units::Current_t;

// Simple macro to exemplify how to fill a GRPMagField object (GRP object containing the information on the magnetic field)

// alienv setenv O2Physics/latest -c root -l ~/cernbox/alice/enigma/macros/createGRPMagFieldObject.C++

void createGRPMagFieldObject(current l3 = 0,     // Ampere
                             current dipole = 0, // Ampere
                             bool isUniform = true,
                             std::string ccdbPath = "file:///tmp/CCDBSnapshot")
{

  GRPMagField grp;
  grp.setL3Current(l3);
  grp.setDipoleCurrent(dipole);
  grp.setFieldUniformity(isUniform);

  CcdbApi api;
  api.init(ccdbPath);
  std::map<std::string, std::string> metadata;
  metadata["responsible"] = "DCS";

  auto year = 2021, month = 10, day = 25, hour = 1, minute = 1, second = 1;
  std::tm tm{}; // get_time does not set all fields hence {}
  tm.tm_year = year - 1900;
  tm.tm_mon = month - 1;
  tm.tm_mday = day;
  tm.tm_hour = hour;
  tm.tm_min = minute;
  tm.tm_sec = second;
  std::time_t tt = std::mktime(&tm);
  std::time_t tStart = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::from_time_t(tt).time_since_epoch()).count();
  timePoint start = tStart;

  std::time_t tEnd = (tStart + 60 * 60 * 24 * 365) * 1000; // 1 year validity, just an example
  timePoint end = tEnd;
  // timePoint end = -1;
  //  long ts = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  if (end < 0) {
    end = (start + 60 * 60 * 10) * 1000; // start + 10h, in ms
  }

  api.storeAsTFileAny(&grp, "GLO/Config/GRPMagField", metadata, start * 1000, end); // making it 1-year valid to be sure we have something
}
