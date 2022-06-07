# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

#o2_add_library(
#   PUBLIC_LINK_LIBRARIES O2::DetectorsBase ROOT::RGL)


#o2_add_test(
#  geometry-misaligner
#  COMPONENT_NAME mch
#  SOURCES testGeometryMisAligner.cxx
#  PUBLIC_LINK_LIBRARIES O2::MCHGeometryTest
#  LABELS muon mch)


o2_add_test_root_macro(
  DisplayTracks.C
  PUBLIC_LINK_LIBRARIES O2::CommonUtils O2::CommonConstants O2::DetectorsBase ROOT::RGL
  LABELS muon mch)