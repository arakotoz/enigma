#include <bitset>
#include <iomanip>
#include <iostream>
#include <string>

#include <TObject.h>
#include <Rtypes.h>

#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "MathUtils/Cartesian.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"
#include "MFTTracking/IOUtils.h"

struct ClusterStruct {
    UShort_t sensor; // sensor id
    UShort_t layer; // layer id
    UShort_t disk; // disk id
    UShort_t half; // half id
    Bool_t isInTrack; // is attached to a track ?
    Double_t measuredGlobalX; // cluster x, global frame (cm)
    Double_t measuredGlobalY; // cluster y, global frame (cm)
    Double_t measuredGlobalZ; // cluster z, global frame (cm)
    Double_t measuredSigmaX2; // cluster, variance x (cm2)
    Double_t measuredSigmaY2; // cluster, variance y (cm2)
    Double_t measuredSigmaZ2; // cluster, variance z (cm2)
};

struct AlignHitStruct {
    UShort_t sensor; // sensor id
    UShort_t layer; // layer id        
    UShort_t disk; // disk id
    UShort_t half; // half id
    Double_t residualX; // (track x - cluster x), global frame (cm)
    Double_t residualY; // (track y - cluster y), global frame (cm)
    Double_t residualZ; // (track z - cluster z), global frame (cm)
};

class Hit
{
    public:
        Hit() = default;
        ~Hit() = default;

        // constructors

        Hit(const Int_t sensor,
            o2::itsmft::ChipMappingMFT mapping,
            const o2::math_utils::Point3D<double>& clusterGlobalPosition);
        Hit(const Int_t sensor,
            o2::itsmft::ChipMappingMFT mapping,
            const double clusterGlobalX,
            const double clusterGlobalY,
            const double clusterGlobalZ);

        // short-named getters

        double clusterGlobalX() const { return mGlobalMeasuredPosition.X(); }
        double clusterGlobalY() const { return mGlobalMeasuredPosition.Y(); }
        double clusterGlobalZ() const { return mGlobalMeasuredPosition.Z(); }
        double clusterSigmaX2() const { return mMeasuredSigmaX2; }
        double clusterSigmaY2() const { return mMeasuredSigmaY2; }
        double clusterSigmaZ2() const { return mMeasuredSigmaZ2; }
        UShort_t disk() const { return mDisk; }
        UShort_t half() const { return mHalf; }
        Bool_t isInTrack() const { return mIsInTrack; }
        UShort_t layer() const { return mLayer; }
        double residualX() const; // (cm)
        double residualY() const; // (cm)
        double residualZ() const; // (cm)
        UShort_t sensor() const { return mSensor; }
        double trackGlobalX() const { return mGlobalRecoPosition.X(); }
        double trackGlobalY() const { return mGlobalRecoPosition.Y(); }
        double trackGlobalZ() const { return mGlobalRecoPosition.Z(); }

        // getters

        AlignHitStruct getAlignHitStruct() const;
        o2::math_utils::Point3D<double> getClusterGlobalXYZ() const
        { 
            return mGlobalMeasuredPosition; 
        }
        o2::math_utils::Point3D<double>& getClusterGlobalXYZ()
        { 
            return mGlobalMeasuredPosition; 
        }
        ClusterStruct getClusterStruct() const;
        o2::math_utils::Point3D<double> getTrackGlobalXYZ() const
        { 
            return mGlobalRecoPosition; 
        }
        o2::math_utils::Point3D<double>& getTrackGlobalXYZ()
        { 
            return mGlobalRecoPosition; 
        }

        // print

        void print() const;

        // setters

        void setClusterGlobalPosition(
            const o2::math_utils::Point3D<double>& position)
        {
            mGlobalMeasuredPosition = position;
        }
        void setClusterGlobalPosition(const double x,
            const double y, const double z)
        {
            setClusterGlobalX(x);
            setClusterGlobalY(y);
            setClusterGlobalZ(z);
        }
        void setClusterGlobalX(const double x)
        {
            mGlobalMeasuredPosition.SetX(x);
        }
        void setClusterGlobalY(const double y)
        {
            mGlobalMeasuredPosition.SetY(y);
        }
        void setClusterGlobalZ(const double z)
        {
            mGlobalMeasuredPosition.SetZ(z);
        }
        void setIsInTrack(){ mIsInTrack = kTRUE; }
        void setMeasuredErrors(double sx2, double sy2, double sz2)
        {
            setMeasuredSigmaX2(sx2);
            setMeasuredSigmaY2(sy2);
            setMeasuredSigmaZ2(sz2);
        }
        void setMeasuredSigmaX2(double v) { mMeasuredSigmaX2 = v; }
        void setMeasuredSigmaY2(double v) { mMeasuredSigmaY2 = v; }
        void setMeasuredSigmaZ2(double v) { mMeasuredSigmaZ2 = v; }
        void setSensor(Int_t sensor, o2::itsmft::ChipMappingMFT mapping);
        void setTrackGlobalPosition(
            const o2::math_utils::Point3D<double>& position)
        {
            mGlobalRecoPosition = position;
        }
        void setTrackGlobalPosition(const double x,
            const double y, const double z)
        {
            setTrackGlobalX(x);
            setTrackGlobalY(y);
            setTrackGlobalZ(z);
        }
        void setTrackGlobalX(const double x)
        {
            mGlobalRecoPosition.SetX(x);
        }
        void setTrackGlobalY(const double y)
        {
            mGlobalRecoPosition.SetY(y);
        }
        void setTrackGlobalZ(const double z)
        {
            mGlobalRecoPosition.SetZ(z);
        }

    protected:
        // sensor id
        UShort_t mSensor = 0;
        // layer id
        UShort_t mLayer = 0;
        // disk id
        UShort_t mDisk = 0;
        // half id
        UShort_t mHalf = 0;
        // does this hit belong to a track ?
        Bool_t mIsInTrack = kFALSE;
        // Cartesian position (cm, in Global frame) of the reconstructed track
        // analytically propagated to the z position of the cluster
        o2::math_utils::Point3D<double> mGlobalRecoPosition; 
        // Cartesian position (cm, in Global frame) of the cluster
        o2::math_utils::Point3D<double> mGlobalMeasuredPosition;
        // Errors on the cluster cartesian position (cm2)
        double mMeasuredSigmaX2 = 0.;
        double mMeasuredSigmaY2 = 0.;
        double mMeasuredSigmaZ2 = 0.;

        ClassDefNV(Hit, 1);
};

//__________________________________________________________________________
Hit::Hit(const Int_t sensor,
         o2::itsmft::ChipMappingMFT mapping,
         const o2::math_utils::Point3D<double>& clusterGlobalPosition)
{
    setSensor(sensor, mapping);
    setClusterGlobalPosition(clusterGlobalPosition);
    setTrackGlobalPosition(0., 0., 0.);
}

//__________________________________________________________________________
Hit::Hit(const Int_t sensor,
         o2::itsmft::ChipMappingMFT mapping,
         const double clusterGlobalX,
         const double clusterGlobalY,
         const double clusterGlobalZ)
{
    setSensor(sensor, mapping);
    setClusterGlobalPosition(clusterGlobalX, clusterGlobalY, clusterGlobalZ);
    setTrackGlobalPosition(0., 0., 0.);
}

//__________________________________________________________________________
double Hit::residualX() const
{
    /// return (track x - cluster x) at this layer, x in global frame (cm)
    double residual = 0;
    if (isInTrack()) {
        residual = mGlobalRecoPosition.X() - mGlobalMeasuredPosition.X();
    } else {
        std::cout << "Hit::residualX() = 0 - "
                  << "WARNING, cluster NOT associated to a track !!!" 
                  << std::endl;
    }
    return residual;
}

//__________________________________________________________________________
double Hit::residualY() const
{
    /// return (track y - cluster y) at this layer, y in global frame (cm)
    double residual = 0;
    if (isInTrack()) {
        residual = mGlobalRecoPosition.Y() - mGlobalMeasuredPosition.Y();
    } else {
        std::cout << "Hit::residualY() = 0 - "
                  << "WARNING, cluster NOT associated to a track !!!" 
                  << std::endl;
    }
    return residual;
}

//__________________________________________________________________________
double Hit::residualZ() const
{
    /// return (track z - cluster z) at this layer, z in global frame (cm)
    double residual = 0;
    if (isInTrack()) {
        residual = mGlobalRecoPosition.Z() - mGlobalMeasuredPosition.Z();
    } else {
        std::cout << "Hit::residualZ() = 0 - "
                  << "WARNING, cluster NOT associated to a track !!!" 
                  << std::endl;
    }
    return residual;
}

//__________________________________________________________________________
AlignHitStruct Hit::getAlignHitStruct() const
{
    AlignHitStruct hit{};
    hit.sensor = sensor();
    hit.layer = layer();
    hit.disk = disk();
    hit.half = half();
    hit.residualX = residualX();
    hit.residualY = residualY();
    hit.residualZ = residualZ();
    return hit;
}

//__________________________________________________________________________
ClusterStruct Hit::getClusterStruct() const
{
    ClusterStruct cluster{};
    cluster.sensor = sensor();
    cluster.layer = layer();
    cluster.disk = disk();
    cluster.half = half();
    cluster.isInTrack = isInTrack();
    cluster.measuredGlobalX = clusterGlobalX();
    cluster.measuredGlobalY = clusterGlobalY();
    cluster.measuredGlobalZ = clusterGlobalZ();
    cluster.measuredSigmaX2 = clusterSigmaX2();
    cluster.measuredSigmaY2 = clusterSigmaY2();
    cluster.measuredSigmaY2 = clusterSigmaZ2();
    return cluster;
}

//__________________________________________________________________________
void Hit::print() const
{
    std::cout << "Hit " << std::noshowpos
              << "s " << std::setw(4) << sensor() << " "
              << "l " << std::setw(2) << layer() << " "
              << "d " << std::setw(2) << disk() << " "
              << "h " << std::setw(1) << half() << " "
              << "cluster (" 
              <<  std::showpos << std::scientific << std::setprecision(3)
              << clusterGlobalX() << ", "
              << clusterGlobalY() << ", "
              << clusterGlobalZ() << ") "
              << std::noshowpos << std::setprecision(3)
              << "err2 ("
              << clusterSigmaX2() << ", "
              << clusterSigmaY2() << ", "
              << clusterSigmaZ2() << ") "
              << "track ("
              << std::showpos << std::setprecision(3)
              << trackGlobalX() << ", "
              << trackGlobalY() << ", "
              << trackGlobalZ() << ") "
              << std::noshowpos << std::setprecision(6)
              << std::endl;
}

//__________________________________________________________________________
void Hit::setSensor(Int_t sensor, o2::itsmft::ChipMappingMFT mapping)
{
    if (sensor < mapping.getNChips()){
        o2::itsmft::MFTChipMappingData chipMapping = (mapping.getChipMappingData())[sensor];
        mSensor = chipMapping.globalChipSWID;
        mLayer = (UShort_t)chipMapping.layer;
        mDisk = (UShort_t)chipMapping.disk;
        mHalf = (UShort_t)chipMapping.half;
    } else {
        std::cout << "Hit::setSensor() - "
                  << "ERROR, sensor id " << sensor 
                  << " >= " << mapping.getNChips() << std::endl;
    }
}


//_________________________________________________________
/// convert compact clusters to 3D spacepoints into std::vector<o2::BaseCluster<float>>
void convertCompactClusters(std::vector<o2::itsmft::CompClusterExt> clusters,
                            o2::mft::GeometryTGeo* geom,
                            o2::itsmft::ChipMappingMFT chipMappingMFT,
                            std::vector<Hit>& output,
                            o2::itsmft::TopologyDictionary& dict)
{
    // inspired from Detectors/ITSMFT/MFT/tracking/src/IOUtils.cxx
    for (auto& c : clusters) {
        auto chipID = c.getChipID();
        auto pattID = c.getPatternID();
        o2::math_utils::Point3D<double> locXYZ;
        double sigmaX2 = o2::mft::ioutils::DefClusError2Row, sigmaY2 = o2::mft::ioutils::DefClusError2Col; //Dummy COG errors (about half pixel size)
        if (pattID != o2::itsmft::CompCluster::InvalidPatternID) {
            sigmaX2 = dict.getErr2X(pattID); // ALPIDE local Y coordinate => MFT global X coordinate (ALPIDE rows)
            sigmaY2 = dict.getErr2Z(pattID); // ALPIDE local Z coordinate => MFT global Y coordinate (ALPIDE columns)
            locXYZ = dict.getClusterCoordinates(c);
        }
        auto gloXYZ = geom->getMatrixL2G(chipID) * locXYZ; // Transformation to the local --> global
        auto& cl3d = output.emplace_back(c.getSensorID(), chipMappingMFT, gloXYZ);
        double cook = 1.0; // WARNING!! COOKED
        cl3d.setMeasuredErrors(cook * sigmaX2, cook * sigmaY2, 0.);
    }
}
