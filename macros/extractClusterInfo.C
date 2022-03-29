#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <TChain.h>
#include <Rtypes.h>

#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"

#include "recoInfo.h"

#endif

void extractClusterInfo()
{
    // geometry

    const bool applyMisalignment = true;
    const bool preferAlignedFile = true; // o2sim_geometry-aligned.root
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

    // cluster chain

    TChain mftclusterChain("o2sim");
    mftclusterChain.Add("/Users/andry/cernbox/alice/mft/pilotbeam/505713/a_raw_0110_tf_053-outdir/mftclusters.root");

    std::vector<o2::itsmft::CompClusterExt> compClusters, *compClustersP = &compClusters;
    std::vector<Hit> mftHits;
    mftclusterChain.SetBranchAddress("MFTClusterComp", &compClustersP);
    Int_t nEntriesClusterChain = mftclusterChain.GetEntries();
    std::cout << "Number of cluster entries = " << nEntriesClusterChain << std::endl;

    // track chain

    TChain mfttrackChain("o2sim");
    mfttrackChain.Add("/Users/andry/cernbox/alice/mft/pilotbeam/505713/a_raw_0110_tf_053-outdir/mfttracks.root");

    std::vector<o2::mft::TrackMFT> mftTracks, *mftTracksP = &mftTracks;
    std::vector<int> trackClusterRefs, *trackClusterRefsP = &trackClusterRefs;
    mfttrackChain.SetBranchAddress("MFTTrack", &mftTracksP);
    mfttrackChain.SetBranchAddress("MFTTrackClusIdx", &trackClusterRefsP);
    Int_t nEntriesTrackChain = mfttrackChain.GetEntries();
    std::cout << "Number of track entries = " << nEntriesTrackChain << std::endl;

    assert(nEntriesClusterChain == nEntriesTrackChain);
    Int_t nEntries = nEntriesClusterChain;

    // loop on both chains

    Int_t trackIdx = 0;
    Int_t nclsTotal = 0;
    Int_t nclsInTracks = 0;
    for (Int_t ii = 0; ii < nEntries; ii++ ) {

        mftclusterChain.GetEntry(ii);
        mfttrackChain.GetEntry(ii);

        // Cache compact clusters

        mftHits.clear();
        mftHits.reserve(compClusters.size());

        // loop on clusters

        for (auto& c : compClusters) {
            mftHits.emplace_back(c, geom, chipMappingMFT, dict);
        }

        // loop on tracks

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
            trackIdx++;
            nclsInTracks += ncls;
        }

        if (ii % 10 == 0) {
            std::cout << "###### Entry " << ii 
                      << " found " << mftHits.size() 
                      << " MFT clusters" << std::endl;
            Int_t index = 0;
            std::cout << "---> mftHits[" << index << "]" << std::endl; 
            mftHits[index].print();
            index = mftHits.size()-1;
            std::cout << "---> mftHits[" << index << "]" << std::endl; 
            mftHits[index].print();    
        }
        nclsTotal += mftHits.size();
    }

    std::cout << "============= SUMMARY ============= " << std::endl;
    std::cout << "Total nb clusters : \t" <<  nclsTotal << std::endl;
    std::cout << "Total nb of tracks : \t" << trackIdx+1 << std::endl;
    std::cout << "Total nb clusters in tracks: \t" <<  nclsInTracks << std::endl;
}

