#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <TChain.h>
#include <Rtypes.h>

#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
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

    // chain

    TChain mftclusterChain("o2sim");
    mftclusterChain.Add("/Users/andry/cernbox/alice/mft/pilotbeam/505713/a_raw_0110_tf_053-outdir/mftclusters.root");

    std::vector<o2::itsmft::CompClusterExt> compClusters, *compClustersP = &compClusters;
    std::vector<Hit> mftHits;
    mftclusterChain.SetBranchAddress("MFTClusterComp", &compClustersP);
    Int_t nEntries = mftclusterChain.GetEntries();
    std::cout << "Number of entries = " << nEntries << std::endl;

    for (Int_t ii = 0; ii < nEntries; ii++ ) {
        mftclusterChain.GetEntry(ii);
        for (auto& c : compClusters) {
            mftHits.emplace_back(c, geom, chipMappingMFT, dict);
        }
    }
    std::cout << "Found " << mftHits.size() << " MFT clusters" << std::endl;
    Int_t index = 0;
    std::cout << "---> mftHits[" << index << "]" << std::endl; 
    mftHits[index].print();
    index = mftHits.size()-1;
    std::cout << "---> mftHits[" << index << "]" << std::endl; 
    mftHits[index].print();
}

