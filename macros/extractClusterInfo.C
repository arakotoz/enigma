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
#include<TString.h>
#include <TTree.h>
#include <Rtypes.h>

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

void extractClusterInfo(const Bool_t doVerbosePrint = true, 
                        const Int_t printPeriod = 10, 
                        const Int_t fileStop = 4315,
                        const bool preferAlignedFile = true,
                        const bool keepClustersInTracksOnly = true
                        )
{
    std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

    // geometry

    const bool applyMisalignment = false;
    o2::base::GeometryManager::loadGeometry("", applyMisalignment, preferAlignedFile);
    o2::mft::GeometryTGeo *geom = o2::mft::GeometryTGeo::Instance();
    geom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G));

    //dictionary

    std::string dictFileName = "MFTdictionary.bin";
    std::ifstream dictFile(dictFileName.c_str());
    if (!dictFile.good()) {
        std::cout << "Error: MFT dictionnary file " << dictFileName << " not found!" << std::endl;
        return;
    }
    o2::itsmft::TopologyDictionary dict;
    dict.readBinaryFile(dictFileName);

    // mapping

    o2::itsmft::ChipMappingMFT chipMappingMFT;

    // cluster and track chains

    std::string generalPath = "/Users/andry/cernbox/alice/mft/pilotbeam/505713/";
    std::string alignStatus = "";
    if (preferAlignedFile || applyMisalignment) {
        alignStatus = "prealigned";
    } else {
        alignStatus = "idealgeo";
    }
    generalPath += alignStatus;


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

	TTree* tree = new TTree("recoInfo","the reco info tree");
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
	tree->Branch("measuredSigmaX2", &hitInfo.measuredSigmaX2, "measuredSigmaX2/D");
	tree->Branch("measuredSigmaY2", &hitInfo.measuredSigmaY2, "measuredSigmaY2/D");
	tree->Branch("measuredSigmaZ2", &hitInfo.measuredSigmaZ2, "measuredSigmaZ2/D");
	tree->Branch("recoGlobalX", &hitInfo.recoGlobalX, "recoGlobalX/D");
	tree->Branch("recoGlobalY", &hitInfo.recoGlobalY, "recoGlobalY/D");
	tree->Branch("recoGlobalZ", &hitInfo.recoGlobalZ, "recoGlobalZ/D");
	tree->Branch("residualX", &hitInfo.residualX, "residualX/D");
	tree->Branch("residualY", &hitInfo.residualY, "residualY/D");
	tree->Branch("residualZ", &hitInfo.residualZ, "residualZ/D");

    // loop on both chains

    Int_t trackIdx = 0;
    Int_t nclsTotal = 0;
    Int_t nclsInTracks = 0;
    Int_t nTrackCA = 0;
    Int_t nRofNoTrack = 0;
    Int_t current_strobe = 0;

    for (Int_t irof = 0; irof < nRof; irof++ ) {

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

        if ( current_strobe >= clustersRof.size() ) {
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
                assert(offset+icls < compClusters.size());
                track.propagateParamToZlinear(mftHits[offset+icls].clusterGlobalZ());
                mftHits[offset+icls].setTrackGlobalPosition(
                    track.getX(),
                    track.getY(),
                    track.getZ()
                );
                mftHits[offset+icls].setTrackIdx(trackIdx);
            }
            if ( track.isCA() ) {
                nTrackCA++;
            }
            trackIdx++;
            nclsInTracks += ncls;
        }

        if ((irof % (nRof/printPeriod) == 0) && doVerbosePrint) {
            std::cout << "\n###### ROF " << irof 
                      << " found " << mftHits.size() 
                      << " MFT clusters, " 
                      << mftTracks.size() 
                      << " MFT tracks"
                      << std::endl;
            Int_t index = 0;
            std::cout << "\t mftHits[" << index << "]\t: " ; 
            mftHits[index].print();
            index = mftHits.size()-1;
            std::cout << "\t mftHits[" << index << "]\t: " ; 
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

    std::chrono::steady_clock::time_point stop_time =  std::chrono::steady_clock::now();

    std::cout << "============= SUMMARY ============= " << std::endl;
    std::cout << "Total nb clusters : \t\t" <<  nclsTotal << std::endl;
    std::cout << "Total nb of tracks : \t\t" << trackIdx+1 << std::endl;
    std::cout << "Total nb clusters in tracks: \t" <<  nclsInTracks << std::endl;
    std::cout << "Total nb of CA tracks : \t" << nTrackCA << std::endl;
    std::cout << "Total nb of no track ROFs : \t" << nRofNoTrack << std::endl;
    std::cout << "----------------------------------- " << endl;
    std::cout << "Execution time: \t\t" 
              << std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time).count()
              << " seconds" << endl;

}

