#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <Rtypes.h>
#include <ROOT/RDataFrame.hxx>

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"

#include "recoInfo.h"

#endif

// root -l
// .L ~/cernbox/alice/enigma/macros/extractClusterInfo.C++
// extractClusterInfo()

// fileStop = 44 for reco-with-mille/old-ctf/pass1
// fileStop = 4315 for prealigned/old-ctf
// fileStop = 173 for reco-with-mille/old-ctf/pass2
// fileStop = 819 for reco-with-mille/new-ctf/new-pass1
// fileStop = 819 for prealigned/new-ctf

void extractClusterInfo(const Bool_t doVerbosePrint = true,
                        const Int_t printPeriod = 10,
                        const Int_t fileStop = 2,
                        const bool preferAlignedFile = true,
                        const bool useMilleAlignment = true,
                        const bool useNewCTFs = true,
                        const bool keepClustersInTracksOnly = true)
{
  ROOT::EnableImplicitMT(0);
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  // geometry

  const bool applyMisalignment = false;
  o2::base::GeometryManager::loadGeometry("", applyMisalignment, preferAlignedFile);
  o2::mft::GeometryTGeo* geom = o2::mft::GeometryTGeo::Instance();
  geom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G));

  // dictionary

  std::string dictFileName = "MFTdictionary.bin";
  if (useNewCTFs) {
    dictFileName = "o2_mft_dictionary.root";
  }
  o2::itsmft::TopologyDictionary* dict = nullptr;
  try {
    dict = o2::itsmft::TopologyDictionary::loadFrom(dictFileName);
  } catch (std::exception e) {
    std::cout << "Error " << e.what() << std::endl;
    return;
  }
  dict->readFromFile(dictFileName);

  // mapping

  o2::itsmft::ChipMappingMFT chipMappingMFT;

  // cluster and track chains

  const int runN = 505713;
  std::string basePath = "/Users/andry/cernbox/alice/mft/pilotbeam";
  std::string alignStatus = "";
  if (preferAlignedFile || applyMisalignment) {
    if (useMilleAlignment) {
      if (useNewCTFs) {
        alignStatus = "reco-with-mille/new-ctf/new-pass1";
      } else {
        alignStatus = "reco-with-mille/old-ctf/pass1";
      }
    } else {
      if (useNewCTFs) {
        alignStatus = "prealigned/new-ctf";
      } else {
        alignStatus = "prealigned/old-ctf";
      }
    }
  } else {
    if (useMilleAlignment) {
      alignStatus = "reco-with-mille/old-ctf/pass2";
    } else {
      alignStatus = "idealgeo/old-ctf";
    }
  }
  std::stringstream generalPathSs;
  generalPathSs << basePath << "/" << runN << "/" << alignStatus;
  std::string generalPath = generalPathSs.str();

  const Int_t fileStart = 1;
  TChain mftclusterChain("o2sim");
  TChain mfttrackChain("o2sim");

  for (Int_t ii = fileStart; ii <= fileStop; ii++) {
    std::stringstream ss;
    if (ii < 100) {
      ss << generalPath << "/"
         << std::setw(3) << std::setfill('0') << ii;
    } else {
      ss << generalPath << "/" << ii;
    }
    std::string filePath = ss.str();
    mftclusterChain.Add(Form("%s/mftclusters.root", filePath.c_str()));
    mfttrackChain.Add(Form("%s/mfttracks.root", filePath.c_str()));
  }
  std::cout << "Number of files per chain = " << fileStop << std::endl;

  std::vector<o2::itsmft::CompClusterExt> compClusters, *compClustersP = &compClusters;
  std::vector<unsigned char> clusterPatterns, *clusterPatternsP = &clusterPatterns;
  std::vector<o2::itsmft::ROFRecord> clustersRof, *clustersRofP = &clustersRof;
  std::vector<Hit> mftHits;
  mftclusterChain.SetBranchAddress("MFTClusterComp", &compClustersP);
  mftclusterChain.SetBranchAddress("MFTClusterPatt", &clusterPatternsP);
  mftclusterChain.SetBranchAddress("MFTClustersROF", &clustersRofP);
  Int_t nEntriesClusterChain = mftclusterChain.GetEntries();
  std::cout << "Number of cluster chain entries = " << nEntriesClusterChain << std::endl;

  std::vector<o2::mft::TrackMFT> mftTracks, *mftTracksP = &mftTracks;
  std::vector<int> trackClusterRefs, *trackClusterRefsP = &trackClusterRefs;
  mfttrackChain.SetBranchAddress("MFTTrack", &mftTracksP);
  mfttrackChain.SetBranchAddress("MFTTrackClusIdx", &trackClusterRefsP);
  Int_t nEntriesTrackChain = mfttrackChain.GetEntries();
  std::cout << "Number of track chain entries = " << nEntriesTrackChain << std::endl;

  assert(nEntriesClusterChain == nEntriesTrackChain);
  const Int_t nRof = nEntriesClusterChain;

  // output tree
  // with general use compression level 505
  // cf https://root.cern.ch/doc/master/structROOT_1_1RCompressionSetting_1_1EDefaults.html#a47faae5d3e4bb7b1941775f764730596aa27e7f29058cc84d676f20aea9b86c30

  TFile hfile(Form("%s/outtree.root", generalPath.c_str()), "recreate", "", 505);

  HitStruct hitInfo;

  TTree* tree = new TTree("recoInfo", "the reco info tree");
  tree->Branch("rofIdx", &hitInfo.rofIdx, "rofIdx/i");
  tree->Branch("bc", &hitInfo.bc, "bc/i");
  tree->Branch("orbit", &hitInfo.orbit, "orbit/i");
  tree->Branch("sensor", &hitInfo.sensor, "sensor/s");
  tree->Branch("layer", &hitInfo.layer, "layer/s");
  tree->Branch("disk", &hitInfo.disk, "disk/s");
  tree->Branch("half", &hitInfo.half, "half/s");
  tree->Branch("trackIdx", &hitInfo.trackIdx, "trackIdx/I");
  tree->Branch("measuredGlobalX", &hitInfo.measuredGlobalX, "measuredGlobalX/D");
  tree->Branch("measuredGlobalY", &hitInfo.measuredGlobalY, "measuredGlobalY/D");
  tree->Branch("measuredGlobalZ", &hitInfo.measuredGlobalZ, "measuredGlobalZ/D");
  tree->Branch("measuredLocalX", &hitInfo.measuredLocalX, "measuredLocalX/D");
  tree->Branch("measuredLocalY", &hitInfo.measuredLocalY, "measuredLocalY/D");
  tree->Branch("measuredLocalZ", &hitInfo.measuredLocalZ, "measuredLocalZ/D");
  tree->Branch("measuredSigmaX2", &hitInfo.measuredSigmaX2, "measuredSigmaX2/D");
  tree->Branch("measuredSigmaY2", &hitInfo.measuredSigmaY2, "measuredSigmaY2/D");
  tree->Branch("measuredSigmaZ2", &hitInfo.measuredSigmaZ2, "measuredSigmaZ2/D");
  tree->Branch("recoGlobalX", &hitInfo.recoGlobalX, "recoGlobalX/D");
  tree->Branch("recoGlobalY", &hitInfo.recoGlobalY, "recoGlobalY/D");
  tree->Branch("recoGlobalZ", &hitInfo.recoGlobalZ, "recoGlobalZ/D");
  tree->Branch("recoLocalX", &hitInfo.recoLocalX, "recoLocalX/D");
  tree->Branch("recoLocalY", &hitInfo.recoLocalY, "recoLocalY/D");
  tree->Branch("recoLocalZ", &hitInfo.recoLocalZ, "recoLocalZ/D");
  tree->Branch("residualX", &hitInfo.residualX, "residualX/D");
  tree->Branch("residualY", &hitInfo.residualY, "residualY/D");
  tree->Branch("residualZ", &hitInfo.residualZ, "residualZ/D");
  tree->Branch("residualLocalX", &hitInfo.residualLocalX, "residualLocalX/D");
  tree->Branch("residualLocalY", &hitInfo.residualLocalY, "residualLocalY/D");
  tree->Branch("residualLocalZ", &hitInfo.residualLocalZ, "residualLocalZ/D");

  // loop on both chains

  Int_t trackIdx = 0;
  Int_t nclsTotal = 0;
  Int_t nclsInTracks = 0;
  Int_t nTrackCA = 0;
  Int_t nRofNoTrack = 0;
  Int_t current_strobe = 0;

  for (Int_t irof = 0; irof < nRof; irof++) {

    mftclusterChain.GetEntry(irof);
    mfttrackChain.GetEntry(irof);

    // Cache compact clusters

    mftHits.clear();
    mftHits.reserve(compClusters.size());
    std::vector<unsigned char>::iterator pattIt = clusterPatterns.begin();

    // loop on clusters

    for (auto& c : compClusters) {
      mftHits.emplace_back(c, pattIt, geom, chipMappingMFT, dict, irof);
    }

    if (current_strobe >= clustersRof.size()) {
      current_strobe = 0;
    }
    assert(current_strobe < clustersRof.size());
    const auto& rofRec = clustersRof[current_strobe].getBCData();
    for (auto icls = 0; icls < compClusters.size(); icls++) {
      mftHits[icls].setOrbit(rofRec.orbit);
      mftHits[icls].setBc(rofRec.bc);
    }

    // loop on tracks

    if (mftTracks.size() == 0) {
      nRofNoTrack++;
    }

    for (auto& track : mftTracks) {
      auto ncls = track.getNumberOfPoints();
      auto offset = track.getExternalClusterIndexOffset();
      for (auto icls = 0; icls < ncls; icls++) {
        auto clsEntry = trackClusterRefs[offset + icls];
        assert(clsEntry < mftHits.size());
        track.propagateParamToZlinear(mftHits[clsEntry].clusterGlobalZ());
        mftHits[clsEntry].setTrackPosition(
          track.getX(),
          track.getY(),
          track.getZ(),
          Hit::isGlobal);
        mftHits[clsEntry].globalToLocal(
          compClusters[clsEntry].getChipID(), geom);
        mftHits[clsEntry].setTrackIdx(trackIdx);
      }
      if (track.isCA()) {
        nTrackCA++;
      }
      trackIdx++;
      nclsInTracks += ncls;
    }

    if ((irof % (nRof / printPeriod) == 0) && doVerbosePrint) {
      std::cout << "\n###### ROF " << irof
                << " found " << mftHits.size()
                << " MFT clusters, "
                << mftTracks.size()
                << " MFT tracks"
                << std::endl;
      Int_t index = 0;
      std::cout << "\t mftHits[" << index << "]\t: ";
      mftHits[index].print();
      index = mftHits.size() - 1;
      std::cout << "\t mftHits[" << index << "]\t: ";
      mftHits[index].print();
    }

    // fill output tree

    for (auto currentHit : mftHits) {
      hitInfo = currentHit.getHitStruct();
      if (keepClustersInTracksOnly) {
        if (currentHit.isInTrack()) {
          tree->Fill();
        }
      } else {
        tree->Fill();
      }
    }

    nclsTotal += mftHits.size();
    current_strobe++;
  }

  // save output tree

  tree->Write();

  std::chrono::steady_clock::time_point stop_time = std::chrono::steady_clock::now();

  std::cout << "============= SUMMARY ============= " << std::endl;
  std::cout << "Total nb clusters : \t\t" << nclsTotal << std::endl;
  std::cout << "Total nb of tracks : \t\t" << trackIdx + 1 << std::endl;
  std::cout << "Total nb clusters in tracks: \t" << nclsInTracks << std::endl;
  std::cout << "Total nb of CA tracks : \t" << nTrackCA << std::endl;
  std::cout << "Total nb of no track ROFs : \t" << nRofNoTrack << std::endl;
  std::cout << "----------------------------------- " << endl;
  std::cout << "Execution time: \t\t"
            << std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time).count()
            << " seconds" << endl;
}
