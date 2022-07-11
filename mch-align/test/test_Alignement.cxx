#include <cmath>
#include <string>
#include <tuple>
#include <vector>
#include <unordered_map>

#include <gsl/span>

// #include <TSystem.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TParameter.h>
#include <TDatabasePDG.h>
#include <Math/Vector4D.h>

#include "Framework/Logger.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "CommonUtils/NameConf.h"
#include "CommonConstants/LHCConstants.h"
#include "MCHGeometryTransformer/Transformations.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackExtrap.h"
#include "ReconstructionDataFormats/TrackMCHMID.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/Digit.h"

#include "SimulationDataFormat/MCCompLabel.h"

//  Test with alignment codes
#include "MCHAlign/Alignment.h"
#include "MCHTracking/Track.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/Cluster.h"
#include "MCHTracking/TrackExtrap.h"


//  Test with alignment codes
#include "MCHAlign/Alignment.h"



using namespace o2;


// const std::string_view prefix{"./simul_results/MCH_100events_fmu"}; // prefix of simulation directory
/*
const uint32_t nOrbitsPerTF = 128;
=======
const std::string_view prefix{"./simul_results/MCH_10events_fmu"}; // prefix of simulation directory
const uint32_t nOrbitsPerTF = 128;

>>>>>>> b38067db918b88ee3beed6b5b242096a4714d21b
// first orbit of the first TF(Trigger Flag) of each run
// for MCH-MID matching
const std::unordered_map<uint32_t, uint32_t> firstTForbit0perRun{
  {505207, 133875},
  {505217, 14225007},
  {505278, 1349340},
  {505285, 1488862},
  {505303, 2615411},
  {505397, 5093945},
  {505404, 19196217},
  {505405, 28537913},
  {505406, 41107641},
  {505413, 452530},
  {505440, 13320708},
  {505443, 26546564},
  {505446, 177711},
  {505548, 88037114},
  {505582, 295044346},
  {505600, 417241082},
  {505623, 10445984},
  {505629, 126979},
  {505637, 338969},
  {505645, 188222},
  {505658, 81044},
  {505669, 328291},
  {505673, 30988},
  {505713, 620506},
  {505720, 5359903}};
<<<<<<< HEAD
*/


struct TrackInfo {
  TrackInfo(const mch::TrackMCH& mch) : mchTrack(mch) {}

  const mch::TrackMCH& mchTrack;
  mch::TrackParam paramAtVertex{};
  double dca = 0.;
  double rAbs = 0.;
  std::vector<const mch::Digit*> digits{};
  double mchTime = -1.;
  double mchTimeRMS = 0.;
  double mchTimeSt12 = -1.;
  double mchTimeRMSSt12 = 0.;
  double mchTimeSt345 = -1.;
  double mchTimeRMSSt345 = 0.;
  std::vector<const mch::Digit*> digitsAtClusterPos{};
  double mchTimeAtClusterPos = -1.;
  double mchTimeRMSAtClusterPos = 0.;
  double mchTimeAtClusterPosSt12 = -1.;
  double mchTimeRMSAtClusterPosSt12 = 0.;
  double mchTimeAtClusterPosSt345 = -1.;
  double mchTimeRMSAtClusterPosSt345 = 0.;
  int midTime = -1;
};


o2::mch::geo::TransformationCreator transformation;
o2::mch::Alignment* test_align = new o2::mch::Alignment();


static const double muMass = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();
uint16_t minNSamplesSignal = 17;
double signalParam[4] = {80., 16., 12., 1.2};
uint16_t minNSamplesBackground = 14;
double backgroundParam[4] = {18., 24., -20., 7.};
int bcIntegrationRange = 6; // time window ([-range, range]) to integrate digits
int minNDigitsSignal = 10; // minimum number of digits passing the signal cuts to select signal events

constexpr double pi() { return 3.14159265358979323846; }
std::tuple<TFile*, TTreeReader*> LoadData(const char* fileName, const char* treeName);
void LoadDigits(TrackInfo& trackInfo, const std::vector<mch::Cluster>& clusters, const std::vector<mch::Digit>& digits,
                bool selectSignal, bool rejectBackground);
void computeMCHTime(const std::vector<const mch::Digit*>& digits, double& mean, double& rms, int deMin = 100, int deMax = 1025);
const dataformats::TrackMCHMID* FindMuon(uint32_t iMCHTrack, const std::vector<dataformats::TrackMCHMID>& muonTracks);
bool ExtrapToVertex(TrackInfo& trackInfo);
bool IsSelected(TrackInfo& trackInfo);
bool IsSignal(TrackInfo& trackInfo);
bool IsReconstructible(TrackInfo& trackInfo);
void CreateHistosAtVertex(std::vector<TH1*>& histos, const char* extension);
void FillHistosAtVertex(const TrackInfo& trackInfo, std::vector<TH1*>& histos);
void DrawHistosAtVertex(std::vector<TH1*> histos[2]);
void CreateTimeHistos(std::vector<TH1*>& histos, const char* extension);
void FillTimeHistos(const std::vector<const mch::Digit*>& digits, double mchTime, double mchTimeRMS, int midTime,
                    gsl::span<TH1*> histos, int deMin = 100, int deMax = 1025);
void DrawTimeHistos(gsl::span<TH1*> histos, const char* extension);
void CreateChargeHistos(std::vector<TH1*>& histos, const char* extension);
void FillChargeHistos(const std::vector<const mch::Digit*>& digits, gsl::span<TH1*> histos, int deMin = 100, int deMax = 1025);
void DrawChargeHistos(gsl::span<TH1*> histos, const char* extension);
void CreateCorrelationHistos(std::vector<TH1*>& histos);
void FillCorrelationHistos(const std::vector<const mch::Digit*>& digits, TH1* hist, double timeRef, int deMin = 100, int deMax = 1025);
void DrawCorrelationHistos(std::vector<TH1*>& histos);
double signalCut(double* x, double* p);
double backgroundCut(double* x, double* p);
void WriteHistos(TFile* f, const char* dirName, const std::vector<TH1*>& histos);
/*
void test_Alignement(int runNumber, std::string mchFileName, bool applyTrackSelection = false,
                   bool selectSignal = false, bool rejectBackground = false, std::string outFileName = "");
*/

mch::Track* MCHFormatConvert(mch::TrackMCH& track, const std::vector<mch::Cluster>& clusters);

//_________________________________________________________________________________________________
void test_Alignement(std::string prefix, std::string mchFileName, std::string recDataFileName, std::string recConsFileName, std::string outFileName = "Alignment_data.root", bool doAlign = true, Double_t weightRecord = 1, bool applyTrackSelection = false,
  bool selectSignal = false, bool rejectBackground = false)
{
  /// show the characteristics of the reconstructed tracks
  /// store the ouput histograms in outFileName if any
  // make sure the correct mapping is loaded
  auto& segmentation = mch::mapping::segmentation(300);
  if (segmentation.nofPads() != 27873) {
    LOG(error) << "wrong mapping implementation";
    LOG(error) << "do gSystem->Load(\"libO2MCHMappingImpl4\") before compiling this macro";
    exit(-1);
  }

  // load magnetic field (from o2sim_grp.root) and geometry (from o2sim_geometry.root)
  // and prepare track extrapolation to vertex (0,0,0)

  // getGRPFileName() if load with no prefix, it will search in the current directory
  // or you can specify a std::string_view object that indicate the files' directory

  cout << "Loading magnetic field and geometry from: " << prefix << endl;
  const auto grp = parameters::GRPObject::loadFrom(base::NameConf::getGRPFileName(prefix));
  base::Propagator::initFieldFromGRP(grp);
  mch::TrackExtrap::setField(); // need to set field value in TrackExtrap.cxx
  mch::TrackExtrap::useExtrapV2();
  base::GeometryManager::loadGeometry(prefix);
  transformation = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);

  // test_align = o2::mch::Alignment();

  // MCH-MID matching, not useful for now since MID is not yet included in this stage
  /*
  // find the first orbit of the first TF of this run
  uint32_t firstTForbit0(0);
  auto itOrbit0 = firstTForbit0perRun.find(runNumber);
  if (itOrbit0 != firstTForbit0perRun.end()) {
    firstTForbit0 = itOrbit0->second;
  } else {
    LOG(warning) << "first orbit not found for this run";
  }
  */

  // load tracks

  // input file: MCH -> mchtracks.root
  cout << "Reading tracks..." <<endl;
  std::string DataFile = Form("%s/%s", prefix.c_str(), mchFileName.c_str());
  auto [fMCH, mchReader] = LoadData(DataFile.c_str(), "o2sim");


  // tree entry do not correspond to tracks, there is only one single entry; all tracks and other data type
  // are stored in a list.
  TTreeReaderValue<std::vector<mch::ROFRecord>> mchROFs = {*mchReader, "trackrofs"};
  TTreeReaderValue<std::vector<mch::TrackMCH>> mchTracks = {*mchReader, "tracks"};
  TTreeReaderValue<std::vector<mch::Cluster>> mchClusters = {*mchReader, "trackclusters"};
  TTreeReaderValue<std::vector<MCCompLabel>> mchLabels = {*mchReader, "tracklabels"}; // for simulation test
  // TTreeReaderValue<std::vector<mch::Digit>> mchDigits = {*mchReader, "trackdigits"};


  // cout << mchReader->GetEntries() << " entries loaded" <<endl;

  /*
  // input file: MUON -> ?
  auto [fMUON, muonReader] = LoadData(muonFileName.c_str(), "o2sim");
  TTreeReaderValue<std::vector<dataformats::TrackMCHMID>> muonTracks = {*muonReader, "tracks"};
  int nTF = muonReader->GetEntries(false);
  if (mchReader->GetEntries(false) != nTF) {
    LOG(error) << mchFileName << " and " << muonFileName << " do not contain the same number of TF";
    exit(-1);
  }
  */

  // open output file if any
  std::string Path_file = Form("./Records/%s", outFileName.c_str());
  TFile* fOut = TFile::Open(Path_file.c_str(), "RECREATE");

  // initialize aligner
  test_align->init(recDataFileName, recConsFileName);

  // load initial global parameters and constraints for Millepede
  // Maybe there will be some initial parameters?
  /*
    // read initial params and constraints from CCDB?

    test_align->InitGlobalParameters(Double_t* par);
  */

  // config for detector constraints:
  /*

  */

  // create root file for saving results of alignment process(GlobalFit)
  // 3 branchs: params, errors, pulls

  int NGlobalPar = test_align->fNGlobal;
  double params[NGlobalPar];
  double errors[NGlobalPar];
  double pulls[NGlobalPar];

  TTree *TreeOut = new TTree("Alignment_records", "Alignment_records");
  TreeOut->Branch("NGlobalPar", &NGlobalPar, "NGlobalPar/I");
  TreeOut->Branch("params", params, "params[NGlobalPar]/D");
  TreeOut->Branch("errors", errors, "errors[NGlobalPar]/D");
  TreeOut->Branch("pulls", pulls, "pulls[NGlobalPar]/D");


  cout << "Start processing..." <<endl;
  // processing for each track
  while(mchReader->Next()){

    // loop for all tracks
    int nb_tracks = mchTracks->size();
    for(int id_track = 0; id_track < nb_tracks; ++ id_track){
      // Data format converision: O2 convention to MCH internal convention:
      // TrackMCH -> Track
      auto mchTrack = mchTracks->at(id_track);
      int nb_clusters = mchTrack.getNClusters();
      cout << "=======================================================================" <<endl;
      cout << "Start processing for track: " << id_track <<endl;
      cout << "=======================================================================" <<endl;

      cout << nb_clusters << " clusters attached for track " << id_track <<endl;


      // Check if current track is valid i.e. correctly assigned to identified particle
      auto mchLabel = mchLabels->at(id_track);
      cout << "Current track is valid? " << mchLabel.isValid() <<endl;
      cout << "Cureent track ID: " << mchLabel.getTrackID() << " for event: " << mchLabel.getEventID() <<endl;
      
      if(!mchLabel.isValid()){
        // skip non valid tracks
        continue;
      } // for simulation test


      // Check if current track is connected with other track
      /*
        cluster uid: the unique ID of the cluster from the chamber ID, detection element ID and cluster index
        a problem occured here while inspecting the clusters for each track
          skip this step
      */

      // Data format conversion
    
      // mch::Track* alignTrack = MCHFormatConvert(mchTrack, *mchClusters);

      /// Format conversion from TrackMCH to Track(MCH internal use)
      cout << "Start track format conversion..." <<endl; 
      mch::Track convertedTrack = mch::Track();
      // int CurrentChamber = alignTrack->getCurrentChamber();

      // Create a reference track param
      /*
      TMatrixD Param0{5, 1};
      Param0.SetMatrixArray(mchTrack.getParameters());
      Param0.Print();
      */
      auto Param0 = mchTrack.getParameters();
      double Z0 = mchTrack.getZ();
      mch::TrackParam extrapParam = mch::TrackParam(Z0, Param0);
      //mch::TrackParam* extrapParam = new mch::TrackParam();

      // Get clusters for current track
      int id_cluster_first = mchTrack.getFirstClusterIdx();
      int id_cluster_last = mchTrack.getLastClusterIdx();
      cout << "Current track's cluster starts at index: " << id_cluster_first << " stops at: " << id_cluster_last << endl;
      cout << "old param => x slope: " << Param0[1] << " y slope: " << Param0[3] <<endl;
      for(int id_cluster = id_cluster_first; id_cluster < id_cluster_last + 1 ; ++id_cluster){
        
        cout << "-----------------------------------------------------------------------" <<endl;

        // mch::TrackParam* extrapParam = new mch::TrackParam(Z0, Param0);

        /*
        extrapParam->setZ(Z0);
        extrapParam->setParameters(Param0);
        extrapParam->print();
        */

        mch::Cluster* cluster = &(mchClusters->at(id_cluster));
        const double Z_cluster = cluster->getZ();
        const int DEId_cluster = cluster->getDEId();
        const int CId_cluster = cluster->getChamberId();
        const int ind_cluster = cluster->getClusterIndex();
        cout << "Cluster: " << id_cluster << " at Z = " << Z_cluster <<endl;
        cout << "Cluster ID  " << "DET: " << DEId_cluster << " Chamber: " << CId_cluster << " Index: " << ind_cluster<<endl;
        
        // cout << "-----------------------------------------------------------------------" <<endl;
        // Perform extrapolation to get all track parameters
        /*
          will there be some extra bias added while processing extrapolation ?
        */
        if(mch::TrackExtrap::extrapToZ(extrapParam, Z_cluster)){
          LOG(info) << Form("Extrapolation succes to cluster: %i", id_cluster);
        }else{
          LOG(fatal) << Form("Extrapolation failed to cluster: %i", id_cluster);
        }

        cout << "new param => x slope: " << extrapParam.getNonBendingSlope() << " y slope: " << extrapParam.getBendingSlope() <<endl;
        cout << "-----------------------------------------------------------------------" <<endl;
        
        // Link current cluster to the extrapolated track param
        extrapParam.setClusterPtr(cluster);

        // Add current param to track
        convertedTrack.setCurrentParam(extrapParam, CId_cluster);
        convertedTrack.addParamAtCluster(extrapParam);
        /*
        auto itParam = convertedTrack.begin();
        for(;itParam!=convertedTrack.end();++itParam){
          const mch::Cluster* itCluster = itParam->getClusterPtr();
          LOG(info) << Form("cluster ID: %i    chamber ID: %i", itCluster->getDEId(), itCluster->getChamberId());
          // LOG(info) << Form("Z : %f", itParam->getZ());
          // cout << &itParam <<endl;
        }
        */
      }

      cout << convertedTrack.getNClusters() << " clusters' param are converted into MCH track format." <<endl;
      cout << "-----------------------------------------------------------------------" <<endl;
      cout << "Track format conversion done."<<endl;


      // cout << alignTrack->getNClusters() <<" clusters for current track"<< endl;
      // auto param_ptr = convertedTrack->begin();
      // const mch::Cluster* cluster = param_ptr->getClusterPtr();
      // cout << "first cluster detector ID: " << cluster->getDEId() <<endl;
      // no problem with alignTrack

      AliMillePedeRecord* mchRecord = test_align->ProcessTrack(convertedTrack, transformation, doAlign, weightRecord);
      
      cout << "=======================================================================" <<endl;
      cout << "Processing done for track: " << id_track <<endl;
      cout << "=======================================================================" <<endl;
      // test_align->ProcessTrack(mchRecord);

      cout << endl;
      cout << endl;
      cout << endl;

    }

    // Data storage finished(fMillepede.fTreeData) with all tracks:
    /*
      track records are stored in test_algin.fMillepede.fTreeData (with branch fRecord and same indexing with tracks)
      and can be read using ReadRecordData(Long_t recID)
    */
    
    cout << "Start alignment process..." <<endl;
    // Process alignment:
    cout << "=======================================================================" <<endl;
    cout << "Start global fitting..."<<endl;
    cout << "=======================================================================" <<endl;
    // Process global fit for each track:
    test_align->GlobalFit(params, errors, pulls);
    // Fill output tree:
    TreeOut->Fill();
    cout << "=======================================================================" <<endl;
    cout << "Global fitting done."<<endl;
    cout << "=======================================================================" <<endl;
    
    

  }
  
  // TreeOut->Print();
  // Close files and store all tracks' records
  test_align->terminate();
  LOG(info)<<"Alignment finished";
  
  if(TreeOut){
    if(fOut->IsWritable()){
      fOut->cd();
      TreeOut->Write();
      TreeOut->Delete();
      if(fOut){
        fOut->Close();
        fOut->Delete();
        LOG(info) << "Data storage done." <<endl;
      }
    }

  }else{
    LOG(fatal) << "Output tree is empty." <<endl;
  }

  LOG(info)<<"Test done!";
  /*
  auto t = transformation(1025);
  TMatrixD mat(3,4);
  t.GetTransformMatrix(mat);
  cout << mat(0,0) <<endl;
  */

}

//_________________________________________________________________________________________________
mch::Track* MCHFormatConvert(mch::TrackMCH& track, const std::vector<mch::Cluster>& clusters)
{

  /// Format conversion from TrackMCH to Track(MCH internal use)
  cout << "Start track format conversion..." <<endl; 
  static mch::Track* convertedTrack = new mch::Track();
  // int CurrentChamber = alignTrack->getCurrentChamber();

  // Create a reference track param
  auto Param0 = track.getParameters();
  double Z0 = track.getZ();
  static mch::TrackParam* extrapParam = new mch::TrackParam(Z0, Param0);


  // Get clusters for current track
  int id_cluster_first = track.getFirstClusterIdx();
  int id_cluster_last = track.getLastClusterIdx();
  cout << "Current track's cluster starts at index: " << id_cluster_first << " stops at: " << id_cluster_last << endl;
  cout << "old param => x slope: " << Param0[1] << " y slope: " << Param0[3] <<endl;
  for(int id_cluster = id_cluster_first; id_cluster < id_cluster_last + 1 ; ++id_cluster){
    
    cout << "-----------------------------------------------------------------------" <<endl;
    auto cluster = clusters.at(id_cluster);
    double Z_cluster = cluster.getZ();
    int DEId_cluster = cluster.getDEId();
    int CId_cluster = cluster.getChamberId();
    int ind_cluster = cluster.getClusterIndex();
    cout << "Cluster: " << id_cluster << " at Z = " << Z_cluster <<endl;
    cout << "Cluster ID  " << "DET: " << DEId_cluster << " Chamber: " << CId_cluster << " Index: " << ind_cluster<<endl;
    
    // cout << "-----------------------------------------------------------------------" <<endl;
    // Perform extrapolation to get all track parameters
    if(mch::TrackExtrap::extrapToZ(*extrapParam, Z_cluster)){
      cout << "Extrapolation succes to cluster: " << id_cluster <<endl;
    }else{
      cout << "Extrapolation failed to cluster: " << id_cluster <<endl;
    }

    cout << "new param => x slope: " << extrapParam->getNonBendingSlope() << " y slope: " << extrapParam->getBendingSlope() <<endl;
    cout << "-----------------------------------------------------------------------" <<endl;
    // Link current cluster to the extrapolated track param
    extrapParam->setClusterPtr(&cluster);

    // Add current param to track
    convertedTrack->setCurrentParam(*extrapParam, CId_cluster);
    convertedTrack->addParamAtCluster(*extrapParam);

  }

  cout << convertedTrack->getNClusters() << " clusters' param are converted into MCH track format." <<endl;
  cout << "-----------------------------------------------------------------------" <<endl;
  cout << "Track format conversion done."<<endl;
  return convertedTrack;
  

}


//_________________________________________________________________________________________________
std::tuple<TFile*, TTreeReader*> LoadData(const char* fileName, const char* treeName)
{
  /// open the input file and get the intput tree

  TFile* f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie()) {
    LOG(error) << "opening file " << fileName << " failed";
    exit(-1);
  }

  TTreeReader* r = new TTreeReader(treeName, f);
  if (r->IsZombie()) {
    LOG(error) << "tree " << treeName << " not found";
    exit(-1);
  }

  return std::make_tuple(f, r);
}

//_________________________________________________________________________________________________
void LoadDigits(TrackInfo& trackInfo, const std::vector<mch::Cluster>& clusters, const std::vector<mch::Digit>& digits,
                bool selectSignal, bool rejectBackground)
{
  /// fill the lists of digits associated to the track

  int nClusterOnTopOfNoDigit(0);

  for (int iCl = trackInfo.mchTrack.getFirstClusterIdx(); iCl <= trackInfo.mchTrack.getLastClusterIdx(); ++iCl) {

    const auto& cluster = clusters[iCl];

    // get the pads at the cluster position
    math_utils::Point3D<float> global{cluster.x, cluster.y, cluster.z};
    auto t = transformation(cluster.getDEId());
    auto local = t^(global);
    int padIDNB(-1), padIDB(-1);
    auto& segmentation = mch::mapping::segmentation(cluster.getDEId());
    bool padsFound = segmentation.findPadPairByPosition(local.x(), local.y(), padIDB, padIDNB);
    bool padFoundNB = padsFound || segmentation.isValid(padIDNB);
    bool padFoundB = padsFound || segmentation.isValid(padIDB);
    if (!padFoundNB && !padFoundB) {
      LOG(warning) << "cluster on top of no pad";
    }

    bool digitFound(false);
    for (uint32_t iDig = 0; iDig < cluster.nDigits; ++iDig) {
      const auto& digit = digits[cluster.firstDigit + iDig];
      if (selectSignal) {
        double nSample = digit.getNofSamples();
        if (digit.getNofSamples() < minNSamplesSignal || digit.getADC() < signalCut(&nSample, signalParam)) {
          continue;
        }
      }
      if (rejectBackground) {
        double nSample = digit.getNofSamples();
        if (digit.getNofSamples() < minNSamplesBackground || digit.getADC() < backgroundCut(&nSample, backgroundParam)) {
          continue;
        }
      }
      trackInfo.digits.push_back(&digit);
      if ((padFoundNB && digit.getPadID() == padIDNB) || (padFoundB && digit.getPadID() == padIDB)) {
        trackInfo.digitsAtClusterPos.push_back(&digit);
        digitFound = true;
      }
    }
    if (!digitFound) {
      ++nClusterOnTopOfNoDigit;
    }
  }

  if (nClusterOnTopOfNoDigit > 0 && trackInfo.midTime >= 0) {
    LOG(warning) << "matched track with " << nClusterOnTopOfNoDigit << "/"
                 << trackInfo.mchTrack.getNClusters() << " clusters on top of no digit";
  }
}

//_________________________________________________________________________________________________
void computeMCHTime(const std::vector<const mch::Digit*>& digits, double& mean, double& rms, int deMin, int deMax)
{
  /// compute the average time and time dispersion of MCH digits

  if (digits.empty()) {
    LOG(error) << "cannot compute mch time";
    return;
  }

  mean = 0.;
  double t2 = 0.;
  double n = 0.;
  for (const auto digit : digits) {
    if (digit->getDetID() >= deMin && digit->getDetID() <= deMax) {
      mean += digit->getTime();
      t2 += static_cast<double>(digit->getTime()) * digit->getTime();
      n += 1.;
    }
  }

  if (n == 0.) {
    LOG(error) << "cannot compute mch time";
    mean = -1.;
    return;
  }

  mean /= n;
  t2 /= n;
  rms = sqrt(t2 - mean * mean);
}

//_________________________________________________________________________________________________
const dataformats::TrackMCHMID* FindMuon(uint32_t iMCHTrack, const std::vector<dataformats::TrackMCHMID>& muonTracks)
{
  /// find the MCH-MID matched track corresponding to this MCH track
  for (const auto& muon : muonTracks) {
    if (muon.getMCHRef().getIndex() == iMCHTrack) {
      return &muon;
    }
  }
  return nullptr;
}

//_________________________________________________________________________________________________
bool ExtrapToVertex(TrackInfo& trackInfo)
{
  /// compute the track parameters at vertex, at DCA and at the end of the absorber
  /// return false if the propagation fails

  // extrapolate to vertex
  trackInfo.paramAtVertex.setZ(trackInfo.mchTrack.getZ());
  trackInfo.paramAtVertex.setParameters(trackInfo.mchTrack.getParameters());
  if (!mch::TrackExtrap::extrapToVertex(trackInfo.paramAtVertex, 0., 0., 0., 0., 0.)) {
    return false;
  }

  // extrapolate to DCA
  mch::TrackParam trackParamAtDCA(trackInfo.mchTrack.getZ(), trackInfo.mchTrack.getParameters());
  if (!mch::TrackExtrap::extrapToVertexWithoutBranson(trackParamAtDCA, 0.)) {
    return false;
  }
  double dcaX = trackParamAtDCA.getNonBendingCoor();
  double dcaY = trackParamAtDCA.getBendingCoor();
  trackInfo.dca = sqrt(dcaX * dcaX + dcaY * dcaY);

  // extrapolate to the end of the absorber
  mch::TrackParam trackParamAtRAbs(trackInfo.mchTrack.getZ(), trackInfo.mchTrack.getParameters());
  if (!mch::TrackExtrap::extrapToZ(trackParamAtRAbs, -505.)) {
    return false;
  }
  double xAbs = trackParamAtRAbs.getNonBendingCoor();
  double yAbs = trackParamAtRAbs.getBendingCoor();
  trackInfo.rAbs = sqrt(xAbs * xAbs + yAbs * yAbs);

  return true;
}

//_________________________________________________________________________________________________
bool IsSelected(TrackInfo& trackInfo)
{
  /// apply standard track selections + pDCA

  static const double sigmaPDCA23 = 80.;
  static const double sigmaPDCA310 = 54.;
  static const double nSigmaPDCA = 6.;
  static const double relPRes = 0.0004;
  static const double slopeRes = 0.0005;

  double thetaAbs = TMath::ATan(trackInfo.rAbs / 505.) * TMath::RadToDeg();
  if (thetaAbs < 2. || thetaAbs > 10.) {
    return false;
  }

  double p = trackInfo.paramAtVertex.p();
  double eta = 0.5 * log((p + trackInfo.paramAtVertex.pz()) / (p - trackInfo.paramAtVertex.pz()));
  if (eta < -4. || eta > -2.5) {
    return false;
  }

  double pDCA = trackInfo.mchTrack.getP() * trackInfo.dca;
  double sigmaPDCA = (thetaAbs < 3) ? sigmaPDCA23 : sigmaPDCA310;
  double nrp = nSigmaPDCA * relPRes * p;
  double pResEffect = sigmaPDCA / (1. - nrp / (1. + nrp));
  double slopeResEffect = 535. * slopeRes * p;
  double sigmaPDCAWithRes = TMath::Sqrt(pResEffect * pResEffect + slopeResEffect * slopeResEffect);
  if (pDCA > nSigmaPDCA * sigmaPDCAWithRes) {
    return false;
  }

  return true;
}

//_________________________________________________________________________________________________
bool IsSignal(TrackInfo& trackInfo)
{
  /// check if the track has still enough digits in time to pass the signal selection

  int nDigits(0);

  for (const auto digit : trackInfo.digits) {
    if (digit->getTime() >= trackInfo.mchTime - bcIntegrationRange && digit->getTime() <= trackInfo.mchTime + bcIntegrationRange) {
      ++nDigits;
    }
  }

  return nDigits > minNDigitsSignal;
}

//_________________________________________________________________________________________________
bool IsReconstructible(TrackInfo& trackInfo)
{
  /// check if the track has still enough digits to be reconstructible

  bool hasDigits[10] = {false, false, false, false, false, false, false, false, false, false};
  for (const auto digit : trackInfo.digits) {
    hasDigits[digit->getDetID() / 100 - 1] = true;
  }

  int nFiredChambersSt45 = 0;
  for (int i = 6; i < 10; ++i) {
    if (hasDigits[i]) {
      ++nFiredChambersSt45;
    }
  }

  return (hasDigits[0] || hasDigits[1]) && (hasDigits[2] || hasDigits[3]) && (hasDigits[4] || hasDigits[5]) && nFiredChambersSt45 >= 3;
}

//_________________________________________________________________________________________________
void CreateHistosAtVertex(std::vector<TH1*>& histos, const char* extension)
{
  /// create single muon histograms at vertex

  histos.emplace_back(new TH1F(Form("pT%s", extension), "pT;p_{T} (GeV/c)", 300, 0., 30.));
  histos.emplace_back(new TH1F(Form("eta%s", extension), "eta;eta", 200, -4.5, -2.));
  histos.emplace_back(new TH1F(Form("phi%s", extension), "phi;phi", 360, 0., 360.));
  histos.emplace_back(new TH1F(Form("dca%s", extension), "DCA;DCA (cm)", 500, 0., 500.));
  histos.emplace_back(new TH1F(Form("pDCA23%s", extension), "pDCA for #theta_{abs} < 3#circ;pDCA (GeV.cm/c)", 2500, 0., 5000.));
  histos.emplace_back(new TH1F(Form("pDCA310%s", extension), "pDCA for #theta_{abs} > 3#circ;pDCA (GeV.cm/c)", 2500, 0., 5000.));
  histos.emplace_back(new TH1F(Form("rAbs%s", extension), "rAbs;R_{abs} (cm)", 1000, 0., 100.));
  histos.emplace_back(new TH1F(Form("nClusters%s", extension), "number of clusters per track;n_{clusters}", 20, 0., 20.));
  histos.emplace_back(new TH1F(Form("chi2%s", extension), "normalized #chi^{2};#chi^{2} / ndf", 500, 0., 50.));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillHistosAtVertex(const TrackInfo& trackInfo, std::vector<TH1*>& histos)
{
  /// fill single muon histograms at vertex

  double thetaAbs = TMath::ATan(trackInfo.rAbs / 505.) * TMath::RadToDeg();
  double pDCA = trackInfo.mchTrack.getP() * trackInfo.dca;
  double pT = sqrt(trackInfo.paramAtVertex.px() * trackInfo.paramAtVertex.px() +
                   trackInfo.paramAtVertex.py() * trackInfo.paramAtVertex.py());
  double p = trackInfo.paramAtVertex.p();
  double eta = 0.5 * log((p + trackInfo.paramAtVertex.pz()) / (p - trackInfo.paramAtVertex.pz()));
  double phi = 180. + atan2(-trackInfo.paramAtVertex.px(), -trackInfo.paramAtVertex.py()) / pi() * 180.;

  histos[0]->Fill(pT);
  histos[1]->Fill(eta);
  histos[2]->Fill(phi);
  histos[3]->Fill(trackInfo.dca);
  if (thetaAbs < 3) {
    histos[4]->Fill(pDCA);
  } else {
    histos[5]->Fill(pDCA);
  }
  histos[6]->Fill(trackInfo.rAbs);
  histos[7]->Fill(trackInfo.mchTrack.getNClusters());
  histos[8]->Fill(trackInfo.mchTrack.getChi2OverNDF());
}

//_________________________________________________________________________________________________
void DrawHistosAtVertex(std::vector<TH1*> histos[2])
{
  /// draw histograms at vertex

  // find the optimal number of pads
  int nPadsx(1), nPadsy(1);
  while ((int)histos[0].size() > nPadsx * nPadsy) {
    if (nPadsx == nPadsy) {
      ++nPadsx;
    } else {
      ++nPadsy;
    }
  }

  // draw histograms
  TCanvas* cHist = new TCanvas("histos", "histos", 10, 10, TMath::Max(nPadsx * 300, 1200), TMath::Max(nPadsy * 300, 900));
  cHist->Divide(nPadsx, nPadsy);
  for (int i = 0; i < (int)histos[0].size(); ++i) {
    cHist->cd((i / nPadsx) * nPadsx + i % nPadsx + 1);
    gPad->SetLogy();
    histos[0][i]->SetStats(false);
    histos[0][i]->SetLineColor(4);
    histos[0][i]->Draw();
    histos[1][i]->SetLineColor(2);
    histos[1][i]->Draw("same");
  }

  // add a legend
  TLegend* lHist = new TLegend(0.5, 0.65, 0.9, 0.8);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos[0][0], Form("%g mch tracks", histos[0][0]->GetEntries()), "l");
  lHist->AddEntry(histos[1][0], Form("%g muon tracks", histos[1][0]->GetEntries()), "l");
  cHist->cd(1);
  lHist->Draw("same");
}

//_________________________________________________________________________________________________
void CreateTimeHistos(std::vector<TH1*>& histos, const char* extension)
{
  /// create track time histograms

  histos.emplace_back(new TH1F(Form("timeResVsMCH%s", extension), "#Deltat vs <MCH time>;#Deltat (BC)", 8001, -2000.25, 2000.25));
  histos.emplace_back(new TH1F(Form("timeResVsMCH%sMatch", extension), "#Deltat vs <MCH time>;#Deltat (BC)", 8001, -2000.25, 2000.25));
  histos.emplace_back(new TH1F(Form("timeResVsMID%s", extension), "#Deltat vs MID time (matched tracks);#Deltat (BC)", 4001, -2000.5, 2000.5));
  histos.emplace_back(new TH1F(Form("timeRMS%s", extension), "MCH time dispersion;#sigmat (BC)", 4001, -0.25, 2000.25));
  histos.emplace_back(new TH1F(Form("timeRMS%sMatch", extension), "MCH time dispersion;#sigmat (BC)", 4001, -0.25, 2000.25));
  histos.emplace_back(new TH1F(Form("timeDiffMCHMID%s", extension), "<MCH time> - MID time (matched tracks);#Deltat (BC)", 8001, -2000.25, 2000.25));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillTimeHistos(const std::vector<const mch::Digit*>& digits, double mchTime, double mchTimeRMS, int midTime,
                    gsl::span<TH1*> histos, int deMin, int deMax)
{
  /// fill track time histograms

  for (const auto digit : digits) {
    if (digit->getDetID() >= deMin && digit->getDetID() <= deMax) {
      histos[0]->Fill(digit->getTime() - mchTime);
    }
  }
  histos[3]->Fill(mchTimeRMS);

  if (midTime >= 0) {
    for (const auto digit : digits) {
      if (digit->getDetID() >= deMin && digit->getDetID() <= deMax) {
        histos[1]->Fill(digit->getTime() - mchTime);
        histos[2]->Fill(digit->getTime() - midTime);
      }
    }
    histos[4]->Fill(mchTimeRMS);
    histos[5]->Fill(mchTime - midTime);
  }
}

//_________________________________________________________________________________________________
void DrawTimeHistos(gsl::span<TH1*> histos, const char* extension)
{
  /// draw track time histograms

  TCanvas* cHist = new TCanvas(Form("cTime%s", extension), Form("cTime%s", extension), 10, 10, 800, 800);
  cHist->Divide(2, 2);
  cHist->cd(1);
  gPad->SetLogy();
  histos[0]->SetStats(false);
  histos[0]->SetLineColor(4);
  histos[0]->Draw();
  histos[1]->SetLineColor(2);
  histos[1]->Draw("same");
  cHist->cd(2);
  gPad->SetLogy();
  histos[2]->SetStats(false);
  histos[2]->SetLineColor(2);
  histos[2]->Draw();
  cHist->cd(3);
  gPad->SetLogy();
  histos[3]->SetStats(false);
  histos[3]->SetLineColor(4);
  histos[3]->Draw();
  histos[4]->SetLineColor(2);
  histos[4]->Draw("same");
  cHist->cd(4);
  gPad->SetLogy();
  histos[5]->SetStats(false);
  histos[5]->SetLineColor(2);
  histos[5]->Draw();

  TLegend* lHist = new TLegend(0.1, 0.8, 0.5, 0.9);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos[0], "all tracks", "l");
  lHist->AddEntry(histos[1], "matched tracks", "l");
  cHist->cd(1);
  lHist->Draw("same");
}

//_________________________________________________________________________________________________
void CreateChargeHistos(std::vector<TH1*>& histos, const char* extension)
{
  /// create track charge histograms

  histos.emplace_back(new TH1F(Form("ADC%s", extension), "ADC;ADC", 10001, -0.5, 10000.5));
  histos.emplace_back(new TH1F(Form("Samples%s", extension), "N samples;N samples", 1024, -0.5, 1023.5));
  histos.emplace_back(new TH2F(Form("ADCvsSample%s", extension), "ADC vs N samples (all tracks);N samples;ADC", 1024, -0.5, 1023.5, 10001, -0.5, 100009.5));

  histos.emplace_back(new TH1F(Form("ADC%sMatch", extension), "ADC;ADC", 10001, -0.5, 10000.5));
  histos.emplace_back(new TH1F(Form("Samples%sMatch", extension), "N samples;N samples", 1024, -0.5, 1023.5));
  histos.emplace_back(new TH2F(Form("ADCvsSample%sMatch", extension), "ADC vs N samples (matched tracks);N samples;ADC", 1024, -0.5, 1023.5, 10001, -0.5, 100009.5));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillChargeHistos(const std::vector<const mch::Digit*>& digits, gsl::span<TH1*> histos, int deMin, int deMax)
{
  /// fill track charge histograms
  for (const auto digit : digits) {
    if (digit->getDetID() >= deMin && digit->getDetID() <= deMax) {
      histos[0]->Fill(digit->getADC());
      histos[1]->Fill(digit->getNofSamples());
      histos[2]->Fill(digit->getNofSamples(), digit->getADC());
    }
  }
}

//_________________________________________________________________________________________________
void DrawChargeHistos(gsl::span<TH1*> histos, const char* extension)
{
  /// draw track charge histograms

  TCanvas* cHist = new TCanvas(Form("cCharge%s", extension), Form("cCharge%s", extension), 10, 10, 800, 800);
  cHist->Divide(2, 2);
  cHist->cd(1);
  gPad->SetLogy();
  histos[0]->SetStats(false);
  histos[0]->SetLineColor(4);
  histos[0]->Draw();
  histos[3]->SetLineColor(2);
  histos[3]->Draw("same");
  cHist->cd(2);
  gPad->SetLogy();
  histos[1]->SetStats(false);
  histos[1]->SetLineColor(4);
  histos[1]->Draw();
  histos[4]->SetLineColor(2);
  histos[4]->Draw("same");
  cHist->cd(3);
  gPad->SetLogz();
  histos[2]->SetStats(false);
  histos[2]->Draw("colz");
  cHist->cd(4);
  gPad->SetLogz();
  histos[5]->SetStats(false);
  histos[5]->Draw("colz");

  TLegend* lHist = new TLegend(0.5, 0.8, 0.9, 0.9);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos[0], "all tracks", "l");
  lHist->AddEntry(histos[3], "matched tracks", "l");
  cHist->cd(1);
  lHist->Draw("same");

  static TF1* fSignal = new TF1("fSignal", signalCut, 0, 1023, 4);
  fSignal->SetParameters(signalParam);
  fSignal->SetLineColor(2);
  static TF1* fBackground = new TF1("fBackground", backgroundCut, 0, 1023, 4);
  fBackground->SetParameters(backgroundParam);
  fBackground->SetLineColor(4);
  cHist->cd(3);
  fSignal->Draw("same");
  fBackground->Draw("same");
  cHist->cd(4);
  fSignal->Draw("same");
  fBackground->Draw("same");
}

//_________________________________________________________________________________________________
void CreateCorrelationHistos(std::vector<TH1*>& histos)
{
  /// create correlation histograms between number of digits and total charge

  histos.emplace_back(new TH2F("ChargevsNDigits", "Charge vs N digits;N digits;ADC", 100, 0, 100, 10000, 0, 100000));
  histos.emplace_back(new TH2F("ChargevsNDigitsSt1", "Charge vs N digits (St1);N digits;ADC", 100, 0, 100, 10000, 0, 100000));
  histos.emplace_back(new TH2F("ChargevsNDigitsSt2", "Charge vs N digits (St2);N digits;ADC", 100, 0, 100, 10000, 0, 100000));
  histos.emplace_back(new TH2F("ChargevsNDigitsSt345", "Charge vs N digits (St345);N digits;ADC", 100, 0, 100, 10000, 0, 100000));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillCorrelationHistos(const std::vector<const mch::Digit*>& digits, TH1* hist, double timeRef, int deMin, int deMax)
{
  /// fill correlation histograms between number of digits and total charge

  uint32_t charge(0);
  int nDigits(0);

  for (const auto digit : digits) {
    if (digit->getDetID() < deMin || digit->getDetID() > deMax) {
      continue;
    }
    if (digit->getTime() < timeRef - bcIntegrationRange || digit->getTime() > timeRef + bcIntegrationRange) {
      continue;
    }
    charge += digit->getADC();
    ++nDigits;
  }

  hist->Fill(nDigits, charge);
}

//_________________________________________________________________________________________________
void DrawCorrelationHistos(std::vector<TH1*>& histos)
{
  /// draw correlation histograms between number of digits and total charge

  TCanvas* cCorr = new TCanvas("cCorr", "cCorr", 10, 10, 800, 800);
  cCorr->Divide(2, 2);
  for (int i = 0; i < 4; ++i) {
    cCorr->cd(i + 1);
    gPad->SetLogz();
    histos[i]->Draw("boxcolz");
  }
}

//_________________________________________________________________________________________________
double signalCut(double* x, double* p)
{
  /// function used to select the signal
  double x0 = pow(p[0] / p[2], 1. / p[3]) + p[1];
  if (x[0] < x0) {
    return p[0];
  } else {
    return p[2] * pow(x[0] - p[1], p[3]);
  }
}

//_________________________________________________________________________________________________
double backgroundCut(double* x, double* p)
{
  /// function used to select the signal
  double x0 = (p[3] * p[2] - p[1] * p[0]) / (p[3] - p[1]);
  if (x[0] < x0) {
    return p[1] * (x[0] - p[0]);
  } else {
    return p[3] * (x[0] - p[2]);
  }
}

//_________________________________________________________________________________________________
void WriteHistos(TFile* f, const char* dirName, const std::vector<TH1*>& histos)
{
  /// write histograms in the subdirectory dirName

  f->mkdir(dirName, dirName, true);
  f->cd(dirName);

  for (auto h : histos) {
    h->Write();
  }
}



///////////////////////////////////////////////////////////////////////////////////////////////////
//                                           Update                                              //
///////////////////////////////////////////////////////////////////////////////////////////////////
/*
1. Issue while checking if the track is connected:
    two clusters identified within the same chamber ?
    the cluster index ordering violated in the above case.

2. Problem while creating MCH type cluster: same namespace used for common data format


*/



