#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <Rtypes.h>
#include <TChain.h>

#include <MathUtils/Cartesian.h>

#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "ReconstructionDataFormats/BaseCluster.h"
#include "MFTTracking/IOUtils.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"

using MFTCluster = o2::BaseCluster<double>;


void convertCompactClusters(std::vector<o2::itsmft::CompClusterExt> clusters,
                            o2::mft::GeometryTGeo* geom,
                            std::vector<MFTCluster>& output,
                            o2::itsmft::TopologyDictionary& dict);


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
    std::vector<MFTCluster> mftClusters;
    mftclusterChain.SetBranchAddress("MFTClusterComp", &compClustersP);
    Int_t nEntries = mftclusterChain.GetEntries();
    std::cout << "Number of entries = " << nEntries << std::endl;

    for (Int_t ii = 0; ii < nEntries; ii++ ) {
        mftclusterChain.GetEntry(ii);
        convertCompactClusters(compClusters, geom, mftClusters, dict);
    }
    std::cout << "Found " << mftClusters.size() << " MFT clusters" << std::endl;
    std::cout << "sensor id (X, Y, Z)/(SigmaY2, SigmaYZ, SigmaZ2) cnt bits" << std::endl;
    Int_t index = 0;
    std::cout << "---> mftClusters[" << index << "]" << std::endl; 
    std::cout << mftClusters[index] << std::endl;
    index = mftClusters.size()-1;
    std::cout << "---> mftClusters[" << index << "]" << std::endl; 
    std::cout << mftClusters[index] << std::endl;
}

//_________________________________________________________
/// convert compact clusters to 3D spacepoints into std::vector<o2::BaseCluster<float>>
void convertCompactClusters(std::vector<o2::itsmft::CompClusterExt> clusters,
                            o2::mft::GeometryTGeo* geom,
                            std::vector<MFTCluster>& output,
                            o2::itsmft::TopologyDictionary& dict)
{
    for (auto& c : clusters) {
        auto chipID = c.getChipID();
        auto pattID = c.getPatternID();
        o2::math_utils::Point3D<double> locXYZ;
        double sigmaX2 = o2::mft::ioutils::DefClusError2Row, sigmaY2 = o2::mft::ioutils::DefClusError2Col;
        if (pattID != o2::itsmft::CompCluster::InvalidPatternID) {
            sigmaX2 = dict.getErr2X(pattID); // ALPIDE local Y coordinate => MFT global X coordinate (ALPIDE rows)
            sigmaY2 = dict.getErr2Z(pattID); // ALPIDE local Z coordinate => MFT global Y coordinate (ALPIDE columns)
            locXYZ = dict.getClusterCoordinates(c);
        }
        auto gloXYZ = geom->getMatrixL2G(chipID) * locXYZ; // Transformation to the local --> global
        auto& cl3d = output.emplace_back(c.getSensorID(), gloXYZ);
        double cook = 1.0; // WARNING!! COOKED
        cl3d.setErrors(cook * sigmaX2, cook * sigmaY2, 0);
    }
}
