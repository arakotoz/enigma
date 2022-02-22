#include <TMath.h>
#include "Math/SMatrix.h"

// shortcuts to type of matrices, for track parameters
using SMatrix44Sym = ROOT::Math::SMatrix<float, 4, 4, ROOT::Math::MatRepSym<float, 4>>;
using SMatrix44Std = ROOT::Math::SMatrix<float, 4>;
using SMatrix4 = ROOT::Math::SVector<float, 4>;
using SVector2 = ROOT::Math::SVector<float, 2>;
using SMatrix22 = ROOT::Math::SMatrix<float, 2>;
using SMatrix24 = ROOT::Math::SMatrix<float, 2, 4>;
using SMatrix42 = ROOT::Math::SMatrix<float, 4, 2>;

// shortcuts to type of matrices, for rotation of the system of coordinates (MCS)
using SMatrix33Std = ROOT::Math::SMatrix<double, 3>;
using SVector3 = ROOT::Math::SVector<double, 3>;

// number of detectors
constexpr int mNDet = 10;

// detector planes position along z
constexpr float mZDet[mNDet] = { 50., 52., 60., 62., 70., 72., 80., 82.,
				 90., 92. };

// default measurement error in x, y, z
float mSigmaX[mNDet] = { 0.05, 0.05, 0.05, 0.05, 0.05, 
			 0.05, 0.05, 0.05, 0.05, 0.05};
float mSigmaY[mNDet] = { 0.05, 0.05, 0.05, 0.05, 0.05, 
			 0.05, 0.05, 0.05, 0.05, 0.05};
float mSigmaZ[mNDet] = { 0.05, 0.05, 0.05, 0.05, 0.05, 
			 0.05, 0.05, 0.05, 0.05, 0.05};

// Multiple Coulomb Scattering sigma_theta
float mSigmaThetaMCS = 0.1 * TMath::DegToRad();  // RMS in a plane

// origin of the tracks
constexpr float mOrigin[3] = { 0., 0., 0. };

// limit in the theta angle of the tracks
constexpr float mThetaMax = 10. * TMath::DegToRad();

// the track
struct Track {

  // clusters
  float mX[mNDet];
  float mY[mNDet];
  float mZ[mNDet];
  // hits
  float mXh[mNDet];
  float mYh[mNDet];
  float mZh[mNDet];

  // track: z, phi, theta
  float mTrkZ;
  float mTrkPhi;
  float mTrkTheta;

  // 
  float mTrackChi2 = 0.;

  // track parameters vector: x, y, phi, tanl
  SMatrix4 mParameters{};
  // track covariance matrix: <x,x>, <x,y>, <x,phi>, ...
  SMatrix44Sym mCovariances{};
  
  // set track parameters
  void setX(float x) { mParameters(0) = x; }
  void setY(float y) { mParameters(1) = y; }
  void setPhi(float phi) { mParameters(2) = phi; }
  void setTanl(float tanl) { mParameters(3) = tanl; }

  // get track parameters
  float getX() const { return mParameters(0); }
  float getY() const { return mParameters(1); }
  float getPhi() const { return mParameters(2); }
  float getTanl() const { return mParameters(3); }

  // get track parameter errors
  float getSigma2X() const { return mCovariances(0, 0); }
  float getSigma2Y() const { return mCovariances(1, 1); }
  float getSigma2Phi() const { return mCovariances(2, 2); }
  float getSigma2Tanl() const { return mCovariances(3, 3); }
  
  void setZ(float z) { mTrkZ = z; }
  float getZ() const { return mTrkZ; }
  
  void setGenPhi(float phi) { mTrkPhi = phi; }
  float getGenPhi() const { return mTrkPhi; }
  
  void setGenTheta(float theta) { mTrkTheta = theta; }
  float getGenTheta() const { return mTrkTheta; }
  
  const SMatrix44Sym& getCovariances() const { return mCovariances; }
  void setCovariances(const SMatrix44Sym& covariances) { mCovariances = covariances; }

  float getTrackChi2() const { return mTrackChi2 / (float)(2 * mNDet - 4); }
  void setTrackChi2(float chi2) { mTrackChi2 = chi2; }

};

//_____________________________________________________________________________
void init(Track& trk) {

  // use the first two points
  auto dx = trk.mX[1] - trk.mX[0];
  auto dy = trk.mY[1] - trk.mY[0];
  auto dz = trk.mZ[1] - trk.mZ[0];
  auto dr = TMath::Sqrt(dx * dx + dy * dy);

  auto x0 = trk.mX[0];
  auto y0 = trk.mY[0];
  auto phi0 = TMath::ATan2(dy, dx);  // [-PI,+PI]
  phi0 += (phi0 < 0.) ? TMath::TwoPi() : 0.;  // bring phi0 to [0,2PI]
  auto tanl0 = dz / dr;

  trk.setX(x0);
  trk.setY(y0);
  trk.setZ(trk.mZ[0]);
  trk.setPhi(phi0);
  trk.setTanl(tanl0);

  auto sigma2DX = mSigmaX[0] * mSigmaX[0] + mSigmaX[1] * mSigmaX[1];
  auto sigma2DY = mSigmaY[0] * mSigmaY[0] + mSigmaY[1] * mSigmaY[1];
  auto sigma2DZ = mSigmaZ[0] * mSigmaZ[0] + mSigmaZ[1] * mSigmaZ[1];

  auto dx2 = dx * dx;
  auto dy2 = dy * dy;
  auto dz2 = dz * dz;
  auto dr2 = dx2 + dy2;
  auto sigma2Phi = (sigma2DX * dy2 + sigma2DY * dx2) / (dr2 * dr2);
  auto sigma2Tanl = (sigma2DZ * dr2 * dr2 + dz2 * (sigma2DX * dx2 + sigma2DY * dy2)) / (dr2 * dr2 * dr2);

  SMatrix44Sym lastParamCov;
  lastParamCov(0, 0) = mSigmaX[0] * mSigmaX[0];
  lastParamCov(1, 1) = mSigmaY[0] * mSigmaY[0];
  lastParamCov(2, 2) = sigma2Phi;
  lastParamCov(3, 3) = sigma2Tanl;

  trk.setCovariances(lastParamCov);
  trk.setTrackChi2(0.);
  
}

//_____________________________________________________________________________
void fit(Track& trk) {
  
  for (int iDet = 1; iDet < mNDet; iDet++) {

    // PROPAGATE
    
    auto dz = trk.mZ[iDet] - trk.getZ();
    auto phi0 = trk.getPhi();
    auto tanl0 = trk.getTanl();
    auto invtanl0 = 1./tanl0;

    trk.mParameters(0) += dz * TMath::Cos(phi0) * invtanl0;
    trk.mParameters(1) += dz * TMath::Sin(phi0) * invtanl0;
    trk.setZ(trk.mZ[iDet]);
    //printf("%f  %f  %f  %f \n", trk.getX(), trk.getY(), trk.mX[iDet], trk.mY[iDet]); 

    // calculate Jacobian
    SMatrix44Std jacob = ROOT::Math::SMatrixIdentity();
    jacob(0, 2) = -dz * TMath::Sin(phi0) * invtanl0;
    jacob(0, 3) = -dz * TMath::Cos(phi0) * invtanl0 * invtanl0;
    jacob(1, 2) = +dz * TMath::Cos(phi0) * invtanl0;
    jacob(1, 3) = -dz * TMath::Sin(phi0) * invtanl0 * invtanl0;

    // extrapolate track parameter covariances to "zEnd"
    trk.setCovariances(ROOT::Math::Similarity(jacob, trk.mCovariances));

    // UPDATE
    
    SMatrix24 H_k;
    H_k(0, 0) = 1.0;
    H_k(1, 1) = 1.0;
    
    SMatrix22 V_k;
    V_k(0, 0) = mSigmaX[iDet] * mSigmaX[iDet];
    V_k(1, 1) = mSigmaY[iDet] * mSigmaY[iDet];
    
    SVector2 m_k(trk.mX[iDet], trk.mY[iDet]), r_k_kminus1;    

    // covariance of residuals
    SMatrix22 invResCov = (V_k + ROOT::Math::Similarity(H_k, trk.mCovariances));
    invResCov.Invert();
    
    // Kalman gain matrix
    SMatrix42 K_k = trk.mCovariances * ROOT::Math::Transpose(H_k) * invResCov;
    
    // update parameters
    r_k_kminus1 = m_k - H_k * trk.mParameters; // Residuals of prediction
    trk.mParameters = trk.mParameters + K_k * r_k_kminus1;

    // update covariances
    SMatrix44Std I = ROOT::Math::SMatrixIdentity();
    SMatrix44Std updCov = (I - K_k * H_k) * trk.mCovariances;
    // std ---> sym matrix
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
	trk.mCovariances(i, j) = updCov(i, j);
      }
    }

    auto addChi2Track = ROOT::Math::Similarity(r_k_kminus1, invResCov);
    trk.mTrackChi2 += addChi2Track;

  }
}
