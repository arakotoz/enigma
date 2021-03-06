// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// $Id$
//
//-----------------------------------------------------------------------------
/// \class MisAligner
///
/// This performs the misalignment on an existing muon arm geometry
/// based on the standard definition of the detector elements in
/// $ALICE_ROOT/MUON/data
///
/// --> User has to specify the magnitude of the alignments, in the Cartesian
/// co-ordiantes (which are used to apply translation misalignments) and in the
/// spherical co-ordinates (which are used to apply angular displacements)
///
/// --> If the constructor is used with no arguments, user has to set
/// misalignment ranges by hand using the methods :
/// SetApplyMisAlig, SetMaxCartMisAlig, SetMaxAngMisAlig, SetXYAngMisAligFactor
/// (last method takes account of the fact that the misalingment is greatest in
/// the XY plane, since the detection elements are fixed to a support structure
/// in this plane. Misalignments in the XZ and YZ plane will be very small
/// compared to those in the XY plane, which are small already - of the order
/// of microns)
///
/// Note : If the detection elements are allowed to be misaligned in all
/// directions, this has consequences for the alignment algorithm
/// (AliMUONAlignment), which needs to know the number of free parameters.
/// Eric only allowed 3 :  x,y,theta_xy, but in principle z and the other
/// two angles are alignable as well.
///
/// \author Bruce Becker, Javier Castillo
//-----------------------------------------------------------------------------

#include "MCHGeometryMisAligner/MisAligner.h"
// #include "AliMUONGeometryTransformer.h"
// #include "AliMUONGeometryModuleTransformer.h"
// #include "AliMUONGeometryDetElement.h"
// #include "AliMUONGeometryBuilder.h"
// #include "AliMpExMap.h"
// #include "AliMpExMapIterator.h"

// #include "AliMathBase.h"
// #include "AliLog.h"

#include <TClonesArray.h>
#include <TGeoMatrix.h>
#include <TMatrixDSym.h>
#include <TMath.h>
#include <TRandom.h>
#include <Riostream.h>

#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include "DetectorsCommonDataFormats/AlignParam.h"

#include "CCDB/CcdbApi.h"

#include "Framework/Logger.h"

ClassImp(o2::mch::geo::MisAligner);

namespace o2::mch::geo
{

bool MisAligner::isMatrixConvertedToAngles(const double* rot, double& psi, double& theta, double& phi) const
{
  /// Calculates the Euler angles in "x y z" notation
  /// using the rotation matrix
  /// Returns false in case the rotation angles can not be
  /// extracted from the matrix
  //
  if (std::abs(rot[0]) < 1e-7 || std::abs(rot[8]) < 1e-7) {
    LOG(ERROR) << "Failed to extract roll-pitch-yall angles!";
    return false;
  }
  psi = std::atan2(-rot[5], rot[8]);
  theta = std::asin(rot[2]);
  phi = std::atan2(-rot[1], rot[0]);
  return true;
}
//______________________________________________________________________________
MisAligner::MisAligner(double cartXMisAligM, double cartXMisAligW, double cartYMisAligM, double cartYMisAligW, double angMisAligM, double angMisAligW)
  : TObject(), mUseUni(kFALSE), mUseGaus(kTRUE), mXYAngMisAligFactor(0.0), mZCartMisAligFactor(0.0)
{
  /// Standard constructor
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 2; j++) {
      mDetElemMisAlig[i][j] = 0.0;
      mModuleMisAlig[i][j] = 0.0;
    }
  }
  mDetElemMisAlig[0][0] = cartXMisAligM;
  mDetElemMisAlig[0][1] = cartXMisAligW;
  mDetElemMisAlig[1][0] = cartYMisAligM;
  mDetElemMisAlig[1][1] = cartYMisAligW;
  mDetElemMisAlig[5][0] = angMisAligM;
  mDetElemMisAlig[5][1] = angMisAligW;
}

//______________________________________________________________________________
MisAligner::MisAligner(double cartMisAligM, double cartMisAligW, double angMisAligM, double angMisAligW)
  : TObject(), mUseUni(kFALSE), mUseGaus(kTRUE), mXYAngMisAligFactor(0.0), mZCartMisAligFactor(0.0)
{
  /// Standard constructor
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 2; j++) {
      mDetElemMisAlig[i][j] = 0.0;
      mModuleMisAlig[i][j] = 0.0;
    }
  }
  mDetElemMisAlig[0][0] = cartMisAligM;
  mDetElemMisAlig[0][1] = cartMisAligW;
  mDetElemMisAlig[1][0] = cartMisAligM;
  mDetElemMisAlig[1][1] = cartMisAligW;
  mDetElemMisAlig[5][0] = angMisAligM;
  mDetElemMisAlig[5][1] = angMisAligW;
}

//______________________________________________________________________________
MisAligner::MisAligner(double cartMisAlig, double angMisAlig)
  : TObject(),
    mUseUni(kTRUE),
    mUseGaus(kFALSE),
    mXYAngMisAligFactor(0.0),
    mZCartMisAligFactor(0.0)
{
  /// Standard constructor
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 2; j++) {
      mDetElemMisAlig[i][j] = 0.0;
      mModuleMisAlig[i][j] = 0.0;
    }
  }
  mDetElemMisAlig[0][1] = cartMisAlig;
  mDetElemMisAlig[1][1] = cartMisAlig;
  mDetElemMisAlig[5][1] = angMisAlig;
}

//_____________________________________________________________________________
MisAligner::MisAligner()
  : TObject(),
    mUseUni(kFALSE),
    mUseGaus(kTRUE),
    mXYAngMisAligFactor(0.0),
    mZCartMisAligFactor(0.0)
{
  /// Default constructor
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 2; j++) {
      mDetElemMisAlig[i][j] = 0.0;
      mModuleMisAlig[i][j] = 0.0;
    }
  }
}

/*
//______________________________________________________________________________
MisAligner::~MisAligner()
{
  /// Destructor
}
 */

//_________________________________________________________________________
void MisAligner::setXYAngMisAligFactor(double factor)
{
  /// Set XY angular misalign factor

  if (TMath::Abs(factor) > 1.0 && factor > 0.) {
    mXYAngMisAligFactor = factor;
    mDetElemMisAlig[3][0] = mDetElemMisAlig[5][0] * factor; // These lines were
    mDetElemMisAlig[3][1] = mDetElemMisAlig[5][1] * factor; // added to keep
    mDetElemMisAlig[4][0] = mDetElemMisAlig[5][0] * factor; // backward
    mDetElemMisAlig[4][1] = mDetElemMisAlig[5][1] * factor; // compatibility
  } else {
    LOG(ERROR) << "Invalid XY angular misalign factor, " << factor;
  }
}

//_________________________________________________________________________
void MisAligner::setZCartMisAligFactor(double factor)
{
  /// Set XY angular misalign factor
  if (TMath::Abs(factor) < 1.0 && factor > 0.) {
    mZCartMisAligFactor = factor;
    mDetElemMisAlig[2][0] = mDetElemMisAlig[0][0];          // These lines were added to
    mDetElemMisAlig[2][1] = mDetElemMisAlig[0][1] * factor; // keep backward compatibility
  } else {
    LOG(ERROR) << Form("Invalid Z cartesian misalign factor, %f", factor);
  }
}

//_________________________________________________________________________
void MisAligner::getUniMisAlign(double cartMisAlig[3], double angMisAlig[3], const double lParMisAlig[6][2]) const
{
  /// Misalign using uniform distribution
  /**
   misalign the centre of the local transformation
   rotation axes :
   fAngMisAlig[1,2,3] = [x,y,z]
   Assume that misalignment about the x and y axes (misalignment of z plane)
   is much smaller, since the entire detection plane has to be moved (the
   detection elements are on a support structure), while rotation of the x-y
   plane is more free.
   */
  cartMisAlig[0] = gRandom->Uniform(-lParMisAlig[0][1] + lParMisAlig[0][0], lParMisAlig[0][0] + lParMisAlig[0][1]);
  cartMisAlig[1] = gRandom->Uniform(-lParMisAlig[1][1] + lParMisAlig[1][0], lParMisAlig[1][0] + lParMisAlig[1][1]);
  cartMisAlig[2] = gRandom->Uniform(-lParMisAlig[2][1] + lParMisAlig[2][0], lParMisAlig[2][0] + lParMisAlig[2][1]);

  angMisAlig[0] = gRandom->Uniform(-lParMisAlig[3][1] + lParMisAlig[3][0], lParMisAlig[3][0] + lParMisAlig[3][1]);
  angMisAlig[1] = gRandom->Uniform(-lParMisAlig[4][1] + lParMisAlig[4][0], lParMisAlig[4][0] + lParMisAlig[4][1]);
  angMisAlig[2] = gRandom->Uniform(-lParMisAlig[5][1] + lParMisAlig[5][0], lParMisAlig[5][0] + lParMisAlig[5][1]); // degrees
}

//_________________________________________________________________________
void MisAligner::getGausMisAlign(double cartMisAlig[3], double angMisAlig[3], const double lParMisAlig[6][2]) const
{
  /// Misalign using gaussian distribution
  /**
   misalign the centre of the local transformation
   rotation axes :
   fAngMisAlig[1,2,3] = [x,y,z]
   Assume that misalignment about the x and y axes (misalignment of z plane)
   is much smaller, since the entire detection plane has to be moved (the
   detection elements are on a support structure), while rotation of the x-y
   plane is more free.
   */
  cartMisAlig[0] = gRandom->Gaus(lParMisAlig[0][0], lParMisAlig[0][1]); //, 3. * lParMisAlig[0][1]);
  cartMisAlig[1] = gRandom->Gaus(lParMisAlig[1][0], lParMisAlig[1][1]); //, 3. * lParMisAlig[1][1]);
  cartMisAlig[2] = gRandom->Gaus(lParMisAlig[2][0], lParMisAlig[2][1]); //, 3. * lParMisAlig[2][1]);

  angMisAlig[0] = gRandom->Gaus(lParMisAlig[3][0], lParMisAlig[3][1]); //, 3. * lParMisAlig[3][1]);
  angMisAlig[1] = gRandom->Gaus(lParMisAlig[4][0], lParMisAlig[4][1]); //, 3. * lParMisAlig[4][1]);
  angMisAlig[2] = gRandom->Gaus(lParMisAlig[5][0], lParMisAlig[5][1]); //, 3. * lParMisAlig[5][1]); // degrees
}

//_________________________________________________________________________
TGeoCombiTrans MisAligner::misAlignDetElem() const
{
  /// Misalign given transformation and return the misaligned transformation.
  /// Use misalignment parameters for detection elements.
  /// Note that applied misalignments are small deltas with respect to the detection
  /// element own ideal local reference frame. Thus deltaTransf represents
  /// the transformation to go from the misaligned d.e. local coordinates to the
  /// ideal d.e. local coordinates.
  /// Also note that this -is not- what is in the ALICE alignment framework known
  /// as local nor global (see MisAligner::MisAlign)

  double cartMisAlig[3] = {0, 0, 0};
  double angMisAlig[3] = {0, 0, 0};

  if (mUseUni) {
    getUniMisAlign(cartMisAlig, angMisAlig, mDetElemMisAlig);
  } else {
    if (!mUseGaus) {
      LOG(WARN) << Form("Neither uniform nor gausian distribution is set! Will use gausian...");
    }
    getGausMisAlign(cartMisAlig, angMisAlig, mDetElemMisAlig);
  }

  TGeoTranslation deltaTrans(cartMisAlig[0], cartMisAlig[1], cartMisAlig[2]);
  TGeoRotation deltaRot;
  deltaRot.RotateX(angMisAlig[0]);
  deltaRot.RotateY(angMisAlig[1]);
  deltaRot.RotateZ(angMisAlig[2]);

  TGeoCombiTrans deltaTransf(deltaTrans, deltaRot);
  // TGeoHMatrix newTransfMat = transform * deltaTransf;

  LOG(INFO) << Form("Rotated DE by %f about Z axis.", angMisAlig[2]);

  return TGeoCombiTrans(deltaTransf);
}

//_________________________________________________________________________
TGeoCombiTrans MisAligner::misAlignModule() const
{
  /// Misalign given transformation and return the misaligned transformation.
  /// Use misalignment parameters for modules.
  /// Note that applied misalignments are small deltas with respect to the module
  /// own ideal local reference frame. Thus deltaTransf represents
  /// the transformation to go from the misaligned module local coordinates to the
  /// ideal module local coordinates.
  /// Also note that this -is not- what is in the ALICE alignment framework known
  /// as local nor global (see MisAligner::MisAlign)

  double cartMisAlig[3] = {0, 0, 0};
  double angMisAlig[3] = {0, 0, 0};

  if (mUseUni) {
    getUniMisAlign(cartMisAlig, angMisAlig, mModuleMisAlig);
  } else {
    if (!mUseGaus) {
      LOG(WARN) << Form("Neither uniform nor gausian distribution is set! Will use gausian...");
    }
    getGausMisAlign(cartMisAlig, angMisAlig, mModuleMisAlig);
  }

  TGeoTranslation deltaTrans(cartMisAlig[0], cartMisAlig[1], cartMisAlig[2]);
  TGeoRotation deltaRot;
  deltaRot.RotateX(angMisAlig[0]);
  deltaRot.RotateY(angMisAlig[1]);
  deltaRot.RotateZ(angMisAlig[2]);

  TGeoCombiTrans deltaTransf(deltaTrans, deltaRot);
  // TGeoHMatrix newTransfMat = transform * deltaTransf;

  LOG(INFO) << Form("Rotated Module by %f about Z axis.", angMisAlig[2]);
  // return deltaTransf;
  return TGeoCombiTrans(deltaTransf);
}

//______________________________________________________________________
// void MisAligner::MisAlign(Bool_t verbose) const
void MisAligner::misAlign(std::vector<o2::detectors::AlignParam>& params, Bool_t verbose) const
{
  /// Takes the internal geometry module transformers, copies them to
  /// new geometry module transformers.
  /// Calculates  module misalignment parameters and applies these
  /// to the new module transformer.
  /// Calculates the module misalignment delta transformation in the
  /// Alice Alignment Framework newTransf = delta * oldTransf.
  /// Add a module misalignment to the new geometry transformer.
  /// Gets the Detection Elements from the module transformer.
  /// Calculates misalignment parameters and applies these
  /// to the local transformation of the Detection Element.
  /// Obtains the new global transformation by multiplying the new
  /// module transformer transformation with the new local transformation.
  /// Applies the new global transform to a new detection element.
  /// Adds the new detection element to a new module transformer.
  /// Calculates the d.e. misalignment delta transformation in the
  /// Alice Alignment Framework (newGlobalTransf = delta * oldGlobalTransf).
  /// Add a d.e. misalignment to the new geometry transformer.
  /// Adds the new module transformer to a new geometry transformer.
  /// Returns the new geometry transformer.

  o2::detectors::DetID detMCH("MCH");

  std::vector<std::vector<int>> DEofHC{{100, 103},
                                       {101, 102},
                                       {200, 203},
                                       {201, 202},
                                       {300, 303},
                                       {301, 302},
                                       {400, 403},
                                       {401, 402},
                                       {500, 501, 502, 503, 504, 514, 515, 516, 517},
                                       {505, 506, 507, 508, 509, 510, 511, 512, 513},
                                       {600, 601, 602, 603, 604, 614, 615, 616, 617},
                                       {605, 606, 607, 608, 609, 610, 611, 612, 613},
                                       {700, 701, 702, 703, 704, 705, 706, 720, 721, 722, 723, 724, 725},
                                       {707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719},
                                       {800, 801, 802, 803, 804, 805, 806, 820, 821, 822, 823, 824, 825},
                                       {807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819},
                                       {900, 901, 902, 903, 904, 905, 906, 920, 921, 922, 923, 924, 925},
                                       {907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919},
                                       {1000, 1001, 1002, 1003, 1004, 1005, 1006, 1020, 1021, 1022, 1023, 1024, 1025},
                                       {1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019}};

  // std::vector<o2::detectors::AlignParam> params;

  o2::detectors::AlignParam lAP;
  for (int hc = 0; hc < 20; hc++) { // module transformers
    // const AliMUONGeometryModuleTransformer* kModuleTransformer =
    //   transformer->GetModuleTransformer(iMt, true);

    // AliMUONGeometryModuleTransformer* newModuleTransformer =
    //   new AliMUONGeometryModuleTransformer(iMt);
    // newGeometryTransformer->AddModuleTransformer(newModuleTransformer);

    // TGeoCombiTrans moduleTransform =
    //   TGeoCombiTrans(*kModuleTransformer->GetTransformation());
    // New module transformation
    LOG(INFO) << "Will MisAlignModule " << hc;
    TGeoCombiTrans localDeltaTransform = misAlignModule();

    // localDeltaTransform.Print();
    std::string sname = fmt::format("MCH/HC{}", hc);
    LOG(INFO) << "symName is " << sname.c_str();
    // newModuleTransformer->SetTransformation(newModuleTransform);
    lAP.setSymName(sname.c_str());
    LOG(INFO) << "local delta params";
    double lPsi, lTheta, lPhi = 0.;
    if (!isMatrixConvertedToAngles(localDeltaTransform.GetRotationMatrix(), lPsi, lTheta, lPhi)) {
      LOG(ERROR) << "Problem extracting angles!";
    }
    LOG(INFO) << fmt::format("{} : {} | X: {:+f} Y: {:+f} Z: {:+f} | pitch: {:+f} roll: {:+f} yaw: {:+f}\n", lAP.getSymName(), lAP.getAlignableID(), localDeltaTransform.GetTranslation()[0],
                             localDeltaTransform.GetTranslation()[1], localDeltaTransform.GetTranslation()[2], lPsi, lTheta, lPhi);
    if (!lAP.setLocalParams(localDeltaTransform)) {
      LOG(ERROR) << "Could not set local params for " << sname.c_str();
    }
    LOG(INFO) << "global delta params";
    LOG(INFO) << fmt::format("{} : {} | X: {:+f} Y: {:+f} Z: {:+f} | pitch: {:+f} roll: {:+f} yaw: {:+f}\n", lAP.getSymName(), lAP.getAlignableID(), lAP.getX(),
                             lAP.getY(), lAP.getZ(), lAP.getPsi(), lAP.getTheta(), lAP.getPhi());
    // lAP.print();
    lAP.applyToGeometry();
    params.emplace_back(lAP);
    for (int de = 0; de < DEofHC[hc].size(); de++) {
      LOG(INFO) << "  Will MisAlignDetElem " << DEofHC[hc][de];
      localDeltaTransform = misAlignDetElem();

      sname = fmt::format("MCH/HC{}/DE{}", hc, DEofHC[hc][de]);
      LOG(INFO) << "  symName is " << sname.c_str();
      lAP.setSymName(sname.c_str());
      LOG(INFO) << "  local delta params";
      if (!isMatrixConvertedToAngles(localDeltaTransform.GetRotationMatrix(), lPsi, lTheta, lPhi)) {
        LOG(ERROR) << "Problem extracting angles!";
      }
      LOG(INFO) << fmt::format("{} : {} | X: {:+f} Y: {:+f} Z: {:+f} | pitch: {:+f} roll: {:+f} yaw: {:+f}\n", lAP.getSymName(), lAP.getAlignableID(), localDeltaTransform.GetTranslation()[0],
                               localDeltaTransform.GetTranslation()[1], localDeltaTransform.GetTranslation()[2], lPsi, lTheta, lPhi);
      if (!lAP.setLocalParams(localDeltaTransform)) {
        LOG(ERROR) << "  Could not set local params for " << sname.c_str();
      }
      LOG(INFO) << "  global delta params";
      LOG(INFO) << fmt::format("  {} : {} | X: {:+f} Y: {:+f} Z: {:+f} | pitch: {:+f} roll: {:+f} yaw: {:+f}\n", lAP.getSymName(), lAP.getAlignableID(), lAP.getX(),
                               lAP.getY(), lAP.getZ(), lAP.getPsi(), lAP.getTheta(), lAP.getPhi());
      params.emplace_back(lAP);
      lAP.applyToGeometry();
    }
    /*
    const std::string& ccdbHost = "http://localhost:8080";
    long tmin = 0;
    long tmax = -1;
    const std::string& objectPath = "";
    const std::string& fileName = "MCHMisAlignment.root";

    if (!ccdbHost.empty()) {
        std::string path = objectPath.empty() ? o2::base::NameConf::getAlignmentPath(detMCH) : objectPath;
        LOGP(info, "Storing alignment object on {}/{}", ccdbHost, path);
        o2::ccdb::CcdbApi api;
        map<string, string> metadata; // can be empty
        api.init(ccdbHost.c_str());   // or http://localhost:8080 for a local installation
        // store abitrary user object in strongly typed manner
        api.storeAsTFileAny(&params, path, metadata, tmin, tmax);
    }

    if (!fileName.empty()) {
        LOGP(info, "Storing ITS alignment in local file {}", fileName);
        TFile algFile(fileName.c_str(), "recreate");
        algFile.WriteObjectAny(&params, "std::vector<o2::detectors::AlignParam>", "alignment");
        algFile.Close();
    }
     */
    // Get delta transformation:
    // Tdelta = Tnew * Told.inverse
    // TGeoHMatrix deltaModuleTransform =
    //   AliMUONGeometryBuilder::Multiply(
    //     newModuleTransform,
    //     kModuleTransformer->GetTransformation()->Inverse());

    // // Create module mis alignment matrix
    // newGeometryTransformer
    //   ->AddMisAlignModule(kModuleTransformer->GetModuleId(), deltaModuleTransform);

    // AliMpExMap* detElements = kModuleTransformer->GetDetElementStore();

    // if (verbose)
    //   LOG(INFO) << Form("%i DEs in old GeometryStore  %i", detElements->GetSize(), iMt);

    // TIter next(detElements->CreateIterator());
    // AliMUONGeometryDetElement* detElement;

    // while ((detElement = static_cast<AliMUONGeometryDetElement*>(next()))) {
    //   /// make a new detection element
    //   AliMUONGeometryDetElement* newDetElement =
    //     new AliMUONGeometryDetElement(detElement->GetId(),
    //                                   detElement->GetVolumePath());

    //   // local transformation of this detection element.
    //   TGeoCombiTrans localTransform = TGeoCombiTrans(*detElement->GetLocalTransformation());
    //   TGeoCombiTrans newLocalTransform = MisAlignDetElem(localTransform);
    //   newDetElement->SetLocalTransformation(newLocalTransform);

    //   // global transformation
    //   TGeoHMatrix newGlobalTransform =
    //     AliMUONGeometryBuilder::Multiply(newModuleTransform,
    //                                      newLocalTransform);
    //   newDetElement->SetGlobalTransformation(newGlobalTransform);

    //   // add this det element to module
    //   newModuleTransformer->GetDetElementStore()->Add(newDetElement->GetId(),
    //                                                   newDetElement);

    //   // In the Alice Alignment Framework misalignment objects store
    //   // global delta transformation
    //   // Get detection "intermediate" global transformation
    //   TGeoHMatrix newOldGlobalTransform = newModuleTransform * localTransform;
    //   // Get detection element global delta transformation:
    //   // Tdelta = Tnew * Told.inverse
    //   TGeoHMatrix deltaGlobalTransform = AliMUONGeometryBuilder::Multiply(
    //     newGlobalTransform,
    //     newOldGlobalTransform.Inverse());

    //   // Create mis alignment matrix
    //   newGeometryTransformer
    //     ->AddMisAlignDetElement(detElement->GetId(), deltaGlobalTransform);
    // }

    if (verbose) {
      LOG(INFO) << Form("MisAligned half chamber %i", hc);
    }
    // newGeometryTransformer->AddModuleTransformer(newModuleTransformer);
  }
  // return newGeometryTransformer;
}

void MisAligner::setAlignmentResolution(const TClonesArray* misAlignArray, int rChId, double rChResX, double rChResY, double rDeResX, double rDeResY)
{
}

// void MisAligner::SetAlignmentResolution(const TClonesArray* misAlignArray, int rChId, double rChResX, double rChResY, double rDeResX, double rDeResY)
// {

//   int chIdMin = (rChId < 0) ? 0 : rChId;
//   int chIdMax = (rChId < 0) ? 9 : rChId;
//   double chResX = (rChResX < 0) ? mModuleMisAlig[0][1] : rChResX;
//   double chResY = (rChResY < 0) ? mModuleMisAlig[1][1] : rChResY;
//   double deResX = (rDeResX < 0) ? mDetElemMisAlig[0][1] : rDeResX;
//   double deResY = (rDeResY < 0) ? mDetElemMisAlig[1][1] : rDeResY;

//   TMatrixDSym mChCorrMatrix(6);
//   mChCorrMatrix[0][0] = chResX * chResX;
//   mChCorrMatrix[1][1] = chResY * chResY;
//   //  mChCorrMatrix.Print();

//   TMatrixDSym mDECorrMatrix(6);
//   mDECorrMatrix[0][0] = deResX * deResX;
//   mDECorrMatrix[1][1] = deResY * deResY;
//   //  mDECorrMatrix.Print();

//   detectors::AlignParam* alignMat = 0x0;

//   for (int chId = chIdMin; chId <= chIdMax; chId++) {
//     TString chName1;
//     TString chName2;

//     chName1 = Form("HC%d", chId * 2);
//     chName2 = Form("HC%d", chId * 2 + 1);

//     for (int i = 0; i < misAlignArray->GetEntries(); i++) {
//       alignMat = (detectors::AlignParam*)misAlignArray->At(i);
//       TString volName(alignMat->getSymName());
//       if ((volName.Contains(chName1) &&
//            ((volName.Last('/') == volName.Index(chName1) + chName1.Length()) ||
//             (volName.Length() == volName.Index(chName1) + chName1.Length()))) ||
//           (volName.Contains(chName2) &&
//            ((volName.Last('/') == volName.Index(chName2) + chName2.Length()) ||
//             (volName.Length() == volName.Index(chName2) + chName2.Length())))) {
//         volName.Remove(0, volName.Last('/') + 1);
//         if (volName.Contains("HC")) {
//           //    alignMat->Print("NULL");
//           alignMat->SetCorrMatrix(mChCorrMatrix);
//         } else if (volName.Contains("DE")) {
//           //    alignMat->Print("NULL");
//           alignMat->SetCorrMatrix(mDECorrMatrix);
//         }
//       }
//     }
//   }
// }
// } // namespace mch
} // namespace o2::mch::geo
