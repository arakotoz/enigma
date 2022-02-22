#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <vector>

#include <TRandom.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TH1F.h>

#include "KF_01.h"

#endif

// track parameters: x, y, slopex, slopey
void KF_01() {

  TH1F* hPullsX = new TH1F("hPullsX", "hPullsX", 500, -10., +10.);
  TH1F* hPullsY = new TH1F("hPullsY", "hPullsY", 500, -10., +10.);
  TH1F* hPullsSlopeX = new TH1F("hPullsSlopeX", "hPullsSlopeX", 500, -10., +10.);
  TH1F* hPullsSlopeY = new TH1F("hPullsSlopeY", "hPullsSlopeY", 500, -10., +10.);
  
  TH1F* hDifX = new TH1F("hDifX", "hDifX", 500, -0.2, +0.2);
  TH1F* hDifY = new TH1F("hDifY", "hDifY", 500, -0.2, +0.2);

  TH1F* hTrkChi2 = new TH1F("hTrkChi2", "hTrkChi2", 500, 0., 10.);
  
  gRandom->SetSeed(0);
  
  // activate the Multiple Coulomb Scattering
  bool withMCS = false;
  mSigmaThetaMCS = 0.01 * TMath::DegToRad();
  
  // sigma of the hit measurement in x, y, z
  for (int iDet = 0; iDet < mNDet; iDet++) {
    mSigmaX[iDet] = 0.015;
    mSigmaY[iDet] = 0.015;
    mSigmaZ[iDet] = 0.000;
  }
  
  // number of tracks to generate
  const int nTracks = 1000000;
  std::vector<Track> tracks;
  tracks.reserve(nTracks);

  // graph to draw the hits of one track
  TGraph2D* gTrk = new TGraph2D();
  gTrk->SetMarkerStyle(20);    
  TCanvas *cgTrk = new TCanvas("cgTrk","One track",0,0,600,400);

  // generate tracks
  Track trk;
  float t, vx, vy, vz, xi, yi, zi, xh, yh, zh, x, y, z, phi, theta;
  float phi_mcs, theta_mcs;
  SMatrix33Std A;
  SVector3 vp, v;
  float cosphi, costheta, sinphi, sintheta;
  float cosphi_mcs, costheta_mcs, sinphi_mcs, sintheta_mcs;
  for (int iTrack = 0; iTrack < nTracks; iTrack++) {
    
    // initial direction at track origin
    phi = gRandom->Rndm() * TMath::TwoPi();
    theta = gRandom->Rndm() * mThetaMax;
    trk.setGenPhi(phi);
    trk.setGenTheta(theta);
    vx = TMath::Sin(theta) * TMath::Cos(phi);
    vy = TMath::Sin(theta) * TMath::Sin(phi);
    vz = TMath::Cos(theta);
    //printf("phi %f theta %f \n", phi, theta);

    // transport track through the detector planes
    xi = mOrigin[0];
    yi = mOrigin[1];
    zi = mOrigin[2];
    for (int iDet = 0; iDet < mNDet; iDet++) {

      // propagate the track to the current detector plane
      t = (mZDet[iDet] - zi) / vz;
      xh = xi + t * vx;
      yh = yi + t * vy;
      zh = mZDet[iDet];

      // set the "hit" coordinates
      trk.mXh[iDet] = xh;
      trk.mYh[iDet] = yh;
      trk.mZh[iDet] = zh;

      // apply measurement error
      x = gRandom->Gaus(xh, mSigmaX[iDet]);
      y = gRandom->Gaus(yh, mSigmaY[iDet]);
      z = zh;

      // control histograms
      hDifX->Fill(x - xh);
      hDifY->Fill(y - yh);
      //printf("Det %02d  x %7.4f  y %7.4f  z %7.4f \n", iDet, x, y, zh);

      // set measured (x, y)
      trk.mX[iDet] = x;
      trk.mY[iDet] = y;
      trk.mZ[iDet] = z;

      // plot control track
      if (iTrack == 0) {
	gTrk->SetPoint(iDet, zh, x, y);
      }

      if (withMCS) {
	// apply MCS
	// transform the SoC with "z" along the track direction
	cosphi = TMath::Cos(phi);
	sinphi = TMath::Sin(phi);
	costheta = TMath::Cos(theta);
	sintheta = TMath::Sin(theta);
	A(0, 0) = cosphi * costheta;
	A(0, 1) = sinphi * costheta;
	A(0, 2) = - sintheta;
	A(1, 0) = - sinphi;
	A(1, 1) = cosphi;
	A(1, 2) = 0.;
	A(2, 0) = cosphi * sintheta;
	A(2, 1) = sinphi * sintheta;
	A(2, 2) = costheta;
	// for the inverse transformation
	A.Invert();

	// generate Gaussian MCS
	phi_mcs = gRandom->Rndm() * TMath::TwoPi();
	theta_mcs = abs(gRandom->Gaus(0., mSigmaThetaMCS));      
	cosphi_mcs = TMath::Cos(phi_mcs);
	sinphi_mcs = TMath::Sin(phi_mcs);
	costheta_mcs = TMath::Cos(theta_mcs);
	sintheta_mcs = TMath::Sin(theta_mcs);
	vp(0) = sintheta_mcs * cosphi_mcs;
	vp(1) = sintheta_mcs * sinphi_mcs;
	vp(2) = costheta_mcs;
	v = A * vp;
	vx = v(0);
	vy = v(1);
	vz = v(2);
	xi = xh;
	yi = yh;
	zi = zh;
      }
      
    }
    if (iTrack == 0) {
      gTrk->Draw("P");
    }
    tracks.push_back(trk);
  }

  // loop over generated tracks and fit
  int ntrk = 0;
  float sigmaX, sigmaY, sigmaSlopeX, sigmaSlopeY;
  float genSlopeX, genSlopeY;
  for (auto& trk : tracks) {
    genSlopeX = TMath::Tan(trk.getGenTheta()) * TMath::Cos(trk.getGenPhi());
    genSlopeY = TMath::Tan(trk.getGenTheta()) * TMath::Sin(trk.getGenPhi());
    init(trk);
    //auto& cov = trk.getCovariances();
    //std::cout << cov << std::endl;
    //printf("trk: %d  %f  %f  %f  %f \n", ntrk, trk.getX(), trk.getY(), trk.getSlopeX(), trk.getTheta());
    //auto& cov = trk.getCovariances();
    //std::cout << cov << std::endl;
    fit(trk);
    sigmaX = TMath::Sqrt(trk.getSigma2X());
    sigmaY = TMath::Sqrt(trk.getSigma2Y());
    sigmaSlopeX = TMath::Sqrt(trk.getSigma2SlopeX());
    sigmaSlopeY = TMath::Sqrt(trk.getSigma2SlopeY());
    //printf("%d   %f  %f  %f  --->  %f  %f  %f  ,  %f  %f \n", ntrk, trk.mX[mNDet - 1], trk.mY[mNDet - 1], trk.mZ[mNDet - 1], trk.getX(), trk.getY(), trk.getZ(), sigmaX, sigmaY);
    //printf("slopes:   %f (%f)   %f (%f) \n", trk.getSlopeX(), genSlopeX, trk.getSlopeY(), genSlopeY);
    // calculate pulls
    auto dx = trk.getX() - trk.mXh[mNDet - 1];
    auto dy = trk.getY() - trk.mYh[mNDet - 1];
    auto dslopex = trk.getSlopeX() - genSlopeX;
    auto dslopey = trk.getSlopeY() - genSlopeY;
    hPullsX->Fill(dx / sigmaX);
    hPullsY->Fill(dy / sigmaY);
    hPullsSlopeX->Fill(dslopex / sigmaSlopeX);
    hPullsSlopeY->Fill(dslopey / sigmaSlopeY);
    hTrkChi2->Fill(trk.getTrackChi2());
    ntrk++;
  }
  printf("Found %d tracks.\n", ntrk);

  TCanvas *c1 = new TCanvas("c1","Pulls",0,0,800,600);
  c1->Divide(2,2);
  c1->cd(1);
  hPullsX->Fit("gaus");
  c1->cd(2);
  hPullsY->Fit("gaus");
  c1->cd(3);
  hPullsSlopeX->Fit("gaus");
  c1->cd(4);
  hPullsSlopeY->Fit("gaus");

  TCanvas *c2 = new TCanvas("c2","Errors",0,0,800,600);
  c2->Divide(2,2);
  c2->cd(1);
  hDifX->Fit("gaus");
  c2->cd(2);
  hDifY->Fit("gaus");

  TCanvas *c3 = new TCanvas("c3","TrkChi2",0,0,800,600);
  c3->Divide(2,2);
  c3->cd(1);
  gPad->SetGridx(true);
  hTrkChi2->Draw();
  
}

