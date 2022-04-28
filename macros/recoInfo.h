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

struct HitStruct {
    UInt_t rofIdx; // readout frame index
    UShort_t sensor; // sensor id
    UShort_t layer; // layer id
    UShort_t disk; // disk id
    UShort_t half; // half id
    Int_t trackIdx; // if attached to a track, track index > -1
    Double_t measuredGlobalX; // cluster x, global frame (cm)
    Double_t measuredGlobalY; // cluster y, global frame (cm)
    Double_t measuredGlobalZ; // cluster z, global frame (cm)
    Double_t measuredSigmaX2; // cluster, variance x (cm2)
    Double_t measuredSigmaY2; // cluster, variance y (cm2)
    Double_t measuredSigmaZ2; // cluster, variance z (cm2)
    Double_t recoGlobalX; // track x, global frame (cm)
    Double_t recoGlobalY; // track y, global frame (cm)
    Double_t recoGlobalZ; // track z, global frame (cm)
    Double_t residualX; // (track x - cluster x), global frame (cm)
    Double_t residualY; // (track y - cluster y), global frame (cm)
    Double_t residualZ; // (track z - cluster z), global frame (cm)
};

class Hit
{
    public:
        Hit();
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
        Hit(o2::itsmft::CompClusterExt cluster,
            o2::mft::GeometryTGeo* geom,
            o2::itsmft::ChipMappingMFT chipMappingMFT,
            o2::itsmft::TopologyDictionary& dict,
            const UInt_t rof = 0);

        // convert a compact cluster to 3D spacepoint stored into a Hit

        void convertCompactCluster(o2::itsmft::CompClusterExt c,
                                   o2::mft::GeometryTGeo* geom,
                                   o2::itsmft::ChipMappingMFT chipMappingMFT,
                                   o2::itsmft::TopologyDictionary& dict);

        // short-named getters

        double clusterGlobalX() const { return mGlobalMeasuredPosition.X(); }
        double clusterGlobalY() const { return mGlobalMeasuredPosition.Y(); }
        double clusterGlobalZ() const { return mGlobalMeasuredPosition.Z(); }
        double clusterSigmaX2() const { return mMeasuredSigmaX2; }
        double clusterSigmaY2() const { return mMeasuredSigmaY2; }
        double clusterSigmaZ2() const { return mMeasuredSigmaZ2; }
        UShort_t disk() const { return mDisk; }
        UShort_t half() const { return mHalf; }
        Bool_t isInTrack() const { return (mTrackIdx > -1) ? kTRUE : kFALSE; }
        UShort_t layer() const { return mLayer; }
        double residualX() const; // (cm)
        double residualY() const; // (cm)
        double residualZ() const; // (cm)
        UInt_t rofIdx() const { return mRofIdx; }
        UShort_t sensor() const { return mSensor; }
        double trackGlobalX() const { return mGlobalRecoPosition.X(); }
        double trackGlobalY() const { return mGlobalRecoPosition.Y(); }
        double trackGlobalZ() const { return mGlobalRecoPosition.Z(); }
        Int_t trackIdx() const { return mTrackIdx; }

        // getters

        o2::math_utils::Point3D<double> getClusterGlobalXYZ() const
        { 
            return mGlobalMeasuredPosition; 
        }
        o2::math_utils::Point3D<double>& getClusterGlobalXYZ()
        { 
            return mGlobalMeasuredPosition; 
        }
        HitStruct getHitStruct() const;
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
        void setMeasuredErrors(double sx2, double sy2, double sz2)
        {
            setMeasuredSigmaX2(sx2);
            setMeasuredSigmaY2(sy2);
            setMeasuredSigmaZ2(sz2);
        }
        void setMeasuredSigmaX2(double v) { mMeasuredSigmaX2 = v; }
        void setMeasuredSigmaY2(double v) { mMeasuredSigmaY2 = v; }
        void setMeasuredSigmaZ2(double v) { mMeasuredSigmaZ2 = v; }
        void setRofIdx(const UInt_t idx){ mRofIdx = idx; }
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
        void setTrackIdx(const Int_t idx){ mTrackIdx = idx; }

    protected:
        // readout frame index
        UInt_t mRofIdx = 0;
        // sensor id
        UShort_t mSensor = 0;
        // layer id
        UShort_t mLayer = 0;
        // disk id
        UShort_t mDisk = 0;
        // half id
        UShort_t mHalf = 0;
        // track index 
        Int_t mTrackIdx = -1;
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
Hit::Hit() 
  : mRofIdx(0),
    mSensor(0),
    mLayer(0),
    mDisk(0),
    mHalf(0),
    mTrackIdx(-1),
    mMeasuredSigmaX2(0.),
    mMeasuredSigmaY2(0.),
    mMeasuredSigmaZ2(0.)
{
    setClusterGlobalPosition(0., 0., 0.);
    setTrackGlobalPosition(0., 0., 0.);
    setMeasuredErrors(0., 0., 0.);
}

//__________________________________________________________________________
Hit::Hit(const Int_t sensor,
         o2::itsmft::ChipMappingMFT mapping,
         const o2::math_utils::Point3D<double>& clusterGlobalPosition)
{
    setSensor(sensor, mapping);
    setClusterGlobalPosition(clusterGlobalPosition);
    setTrackGlobalPosition(0., 0., 0.);
    setMeasuredErrors(0., 0., 0.);
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
    setMeasuredErrors(0., 0., 0.);
}

//__________________________________________________________________________
Hit::Hit(o2::itsmft::CompClusterExt cluster,
         o2::mft::GeometryTGeo* geom,
         o2::itsmft::ChipMappingMFT chipMappingMFT,
         o2::itsmft::TopologyDictionary& dict,
         const UInt_t rof)
{
    setTrackGlobalPosition(0., 0., 0.);
    setMeasuredErrors(0., 0., 0.);
    setRofIdx(rof);
    convertCompactCluster(cluster, geom, chipMappingMFT, dict);
}

//__________________________________________________________________________
double Hit::residualX() const
{
    /// return (track x - cluster x) at this layer, x in global frame (cm)
    double residual = 0;
    if (isInTrack()) {
        residual = mGlobalRecoPosition.X() - mGlobalMeasuredPosition.X();
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
    }
    return residual;
}

//__________________________________________________________________________
HitStruct Hit::getHitStruct() const
{
    HitStruct hit{};
    hit.rofIdx = rofIdx();
    hit.sensor = sensor();
    hit.layer = layer();
    hit.disk = disk();
    hit.half = half();
    hit.trackIdx = trackIdx();
    hit.measuredGlobalX = clusterGlobalX();
    hit.measuredGlobalY = clusterGlobalY();
    hit.measuredGlobalZ = clusterGlobalZ();
    hit.measuredSigmaX2 = clusterSigmaX2();
    hit.measuredSigmaY2 = clusterSigmaY2();
    hit.measuredSigmaZ2 = clusterSigmaZ2();
    hit.recoGlobalX = trackGlobalX(); 
    hit.recoGlobalY = trackGlobalY(); 
    hit.recoGlobalZ = trackGlobalZ(); 
    hit.residualX = residualX();
    hit.residualY = residualY(); 
    hit.residualZ = residualZ(); 
    return hit;
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
              << clusterSigmaZ2() << ") ";
    if ( isInTrack() ) {
        std::cout << "track " << trackIdx() << " ("
              << std::showpos << std::setprecision(3)
              << trackGlobalX() << ", "
              << trackGlobalY() << ", "
              << trackGlobalZ() << ") ";
    }
    std::cout << std::noshowpos << std::setprecision(6)
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
void Hit::convertCompactCluster(o2::itsmft::CompClusterExt c,
                                o2::mft::GeometryTGeo* geom,
                                o2::itsmft::ChipMappingMFT chipMappingMFT,
                                o2::itsmft::TopologyDictionary& dict)
{
    /// convert a compact cluster to 3D spacepoint stored into a Hit
    // lines from Detectors/ITSMFT/MFT/tracking/src/IOUtils.cxx
    auto chipID = c.getChipID();
    auto pattID = c.getPatternID();
    o2::math_utils::Point3D<double> locXYZ;
    //Dummy COG errors (about half pixel size)
    double sigmaX2 = o2::mft::ioutils::DefClusError2Row;
    double sigmaY2 = o2::mft::ioutils::DefClusError2Col; 
    if (pattID != o2::itsmft::CompCluster::InvalidPatternID) {
        // ALPIDE local Y coordinate => MFT global X coordinate (ALPIDE rows)
        sigmaX2 = dict.getErr2X(pattID); 
        // ALPIDE local Z coordinate => MFT global Y coordinate (ALPIDE columns)
        sigmaY2 = dict.getErr2Z(pattID); 
        locXYZ = dict.getClusterCoordinates(c);
    }
    // Transformation to the local --> global
    auto gloXYZ = geom->getMatrixL2G(chipID) * locXYZ;
    setSensor(c.getSensorID(), chipMappingMFT);
    setClusterGlobalPosition(gloXYZ);
    double cook = 1.0; // WARNING!! COOKED
    setMeasuredErrors(cook * sigmaX2, cook * sigmaY2, 0.);
}
