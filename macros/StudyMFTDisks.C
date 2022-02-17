using MFTCluster = o2::BaseCluster<float>;
using MFTTrack = o2::mft::TrackMFT;
using o2::MCTrack;
using o2::itsmft::CompClusterExt;
o2::itsmft::ChipMappingMFT mftChipMapper;
o2::itsmft::TopologyDictionary dict;




std::vector<TH2D*> histoClsXYinDiskNb(5);
std::vector<TH2D*> histoClsXYRedundantinDiskNb(5);

std::vector<TH2D*> histoClsPhiEtainDiskNb(5);
std::vector<TH2D*> histoClsPhiEtaRedundantinDiskNb(5);

std::vector<TH2D*> histoClsXYinLayerNb(10);
std::vector<TH2D*> histoClsPhiEtainLayerNb(10);


std::vector<TH2D*> histoKinPhiEta_ncls(5);

TH2D *histoKinPhiEta=0;

TTree* mftTrackTree=0;

void BookHistos();

o2::math_utils::Point3D<float> convertToGlobalCoordinates(CompClusterExt cluster, o2::mft::GeometryTGeo *gman);

//StudyMFTDisks("outputfile_studyMFTDisks.root", "mftclusters.root", "mfttracks.root")

void StudyMFTDisks(const Char_t *ofname = "outputfile_studyMFTDisks.root", const Char_t *clusterFileName = "mftclusters.root", const Char_t *recoFileName = "mfttracks.root")
{
  BookHistos();



  // Geometry and matrix transformations
  std::string inputGeom = "o2sim_geometry.root";
  //o2::base::GeometryManager::loadGeometry(inputGeom);
  o2::base::GeometryManager::loadGeometry("", true);
  //see here: https://github.com/AliceO2Group/AliceO2/tree/dev/Detectors/Base
  o2::mft::GeometryTGeo *gman = o2::mft::GeometryTGeo::Instance();
  gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));

  //mftTracks
  TFile* trkFileIn = new TFile(recoFileName);
  mftTrackTree = (TTree*)trkFileIn->Get("o2sim");
  std::vector<o2::mft::TrackMFT> trackMFTVec, *trackMFTVecP = &trackMFTVec;
  mftTrackTree->SetBranchAddress("MFTTrack", &trackMFTVecP);

  std::vector<int> trackExtClsVec, *trackExtClsVecP = &trackExtClsVec;
  mftTrackTree->SetBranchAddress("MFTTrackClusIdx", &trackExtClsVecP);
  int nEntries=mftTrackTree->GetEntries();

  // Cluster pattern dictionary
  std::string dictfile = "MFTdictionary.bin";

  std::ifstream file(dictfile.c_str());
  if (file.good()) {
    printf("Running with dictionary: %s \n", dictfile.c_str());
    dict.readBinaryFile(dictfile);
  } else {
    printf("Can not run without dictionary !\n");
    return;
  }

  // Cluster file: initializing cluster tree
  TFile clusterFile(clusterFileName);
  TTree* clsTree = (TTree*)clusterFile.Get("o2sim");
  std::vector<CompClusterExt> clsVec, *clsVecP = &clsVec;
  clsTree->SetBranchAddress("MFTClusterComp", &clsVecP);


  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clsLabels = nullptr; // This goes to global variables
  if (clsTree->GetBranch("MFTClusterMCTruth"))
  { // This goes to LoadClusters
    clsTree->SetBranchAddress("MFTClusterMCTruth", &clsLabels);
  }
  else
  {
    printf("No Monte-Carlo information in this file\n");
  }




  for (int iEntry = 0; iEntry < nEntries ; iEntry++)
  {
    clsTree -> GetEntry(iEntry);
    mftTrackTree->GetEntry(iEntry);

    std::vector<MFTTrack> mMFTTracks;
    std::vector<int> mtrackExtClsIDs;




    mMFTTracks.swap(trackMFTVec);
    mtrackExtClsIDs.swap(trackExtClsVec);
    int srcIDTrue, trkIDTrue, evnIDTrue;
    //printf("nbTrack %lu\n", mMFTTracks.size());
    int iTrack = 0;
    std::array<std::array<int,2>,5> clsEntriesForRedundancy;
    for (auto& mftTrack : mMFTTracks)
    {//loop over the MFT tracks

      for(auto idisk = 0 ; idisk < 5 ; idisk++)
      {
        clsEntriesForRedundancy[idisk]={-1,-1};
      }
      auto ncls = mftTrack.getNumberOfPoints();
      auto offset = mftTrack.getExternalClusterIndexOffset();

      histoKinPhiEta->Fill(mftTrack.getPhi(),mftTrack.getEta());
      if (ncls>=5)
      {
        histoKinPhiEta_ncls[0]->Fill(mftTrack.getPhi(),mftTrack.getEta());
      }
      if (ncls>=6)
      {
        histoKinPhiEta_ncls[1]->Fill(mftTrack.getPhi(),mftTrack.getEta());
      }
      if (ncls>=7)
      {
        histoKinPhiEta_ncls[2]->Fill(mftTrack.getPhi(),mftTrack.getEta());
      }
      if (ncls>=8)
      {
        histoKinPhiEta_ncls[3]->Fill(mftTrack.getPhi(),mftTrack.getEta());
      }
      if (ncls>=9)
      {
        histoKinPhiEta_ncls[4]->Fill(mftTrack.getPhi(),mftTrack.getEta());
      }
      //printf("For now all is fine %d\n", iTrack);

      for (int icls = 0; icls < ncls; ++icls)//cluster loop 1
      {

        auto clsEntry = mtrackExtClsIDs[offset + icls];
        //printf("clsEntry = %d\n", clsEntry);

        auto cluster = clsVec[clsEntry];
        //printf("Entering cluster loop %d\n", icls);
        auto clsLayer =  mftChipMapper.chip2Layer(cluster.getChipID());

        int clsMFTdiskID = clsLayer/2; //entier pour root


        if (clsEntriesForRedundancy[clsMFTdiskID][0]!=-1)
        {
          clsEntriesForRedundancy[clsMFTdiskID][1]=clsEntry;
        }
        else
        {
          clsEntriesForRedundancy[clsMFTdiskID][0]=clsEntry;
        }
        //printf("Entering cluster loop %d\n", icls);
        //Get the global position of the cluster

        // Transformation to the local --> global
        auto gloC = convertToGlobalCoordinates(cluster, gman);


        TVector3 v;
        v.SetXYZ(gloC.x(),gloC.y(),gloC.z());

        double etaCls = v.Eta();
        double phiCls = v.Phi();

        histoClsPhiEtainDiskNb[clsMFTdiskID]->Fill(phiCls,etaCls);
        histoClsPhiEtainLayerNb[clsLayer]->Fill(phiCls,etaCls);

        histoClsXYinDiskNb[clsMFTdiskID]->Fill(gloC.x(),gloC.y());
        histoClsXYinLayerNb[clsLayer]->Fill(gloC.x(),gloC.y());
        //-----------------------Form("histoClsXYinDiskNb_%d", i)
      }//end of cluster loop 1


      for(auto idisk = 0 ; idisk < 5 ; idisk++)
      {
        if ((clsEntriesForRedundancy[idisk][0]!=-1) && (clsEntriesForRedundancy[idisk][1]!=-1))
        {
          o2::itsmft::CompClusterExt cluster1 = clsVec[clsEntriesForRedundancy[idisk][0]];
          o2::itsmft::CompClusterExt cluster2 = clsVec[clsEntriesForRedundancy[idisk][1]];

          auto gloC1 = convertToGlobalCoordinates(cluster1, gman);
          auto gloC2 = convertToGlobalCoordinates(cluster2, gman);

          TVector3 v1;
          TVector3 v2;
          v1.SetXYZ(gloC1.x(),gloC1.y(),gloC1.z());
          v2.SetXYZ(gloC2.x(),gloC2.y(),gloC2.z());

          histoClsPhiEtaRedundantinDiskNb[idisk]->Fill(v1.Phi(),v1.Eta());
          histoClsPhiEtaRedundantinDiskNb[idisk]->Fill(v2.Phi(),v2.Eta());

          histoClsXYRedundantinDiskNb[idisk]->Fill(gloC1.x(),gloC1.y());
          histoClsXYRedundantinDiskNb[idisk]->Fill(gloC2.x(),gloC2.y());
        }
      }


      iTrack++;
    }//end of loop over mftTracks

  }//end of loop on entries




  TFile of(ofname, "RECREATE");//output file
  //Write everything in one output file
  of.cd();

  for(auto i = 0 ; i < 5 ; i++)
  {
    histoClsXYinDiskNb[i]->Write();
    histoClsXYRedundantinDiskNb[i]->Write();
    histoClsPhiEtainDiskNb[i]->Write();
    histoClsPhiEtaRedundantinDiskNb[i]->Write();

    histoKinPhiEta_ncls[i]->Write();
  }

  for(auto i = 0 ; i < 10 ; i++)
  {
    histoClsXYinLayerNb[i]->Write();
    histoClsPhiEtainLayerNb[i]->Write();
  }

  histoKinPhiEta->Write();

  of.Close();

}




o2::math_utils::Point3D<float> convertToGlobalCoordinates(CompClusterExt cluster, o2::mft::GeometryTGeo *gman)
{
  //Get the global position of the cluster
  auto chipID = cluster.getChipID();
  auto pattID = cluster.getPatternID();
  o2::math_utils::Point3D<float> locC;

  if (pattID != o2::itsmft::CompCluster::InvalidPatternID)
  {
    locC = dict.getClusterCoordinates(cluster);
  }

  // Transformation to the local --> global
  auto gloC = gman->getMatrixL2G(chipID) * locC;

  return gloC;
}

void BookHistos()
{
  //Histograms for each MFT disk (disk numbers from 0 to 4)
  for(auto i = 0 ; i < 5 ; i++)
  {
    histoClsXYinDiskNb[i] = new TH2D(Form("histoClsXYinDiskNb_%d", i), Form("Y versus X of clusters in the disk %d", i), 400, -20, 20, 400, -20, 20);
    histoClsXYinDiskNb[i]->SetXTitle("x (cm) position of clusters");
    histoClsXYinDiskNb[i]->SetYTitle("y (cm) position of clusters");
    histoClsXYinDiskNb[i]->Sumw2();

    histoClsXYRedundantinDiskNb[i] = new TH2D(Form("histoClsXYRedundantinDiskNb_%d", i), Form("Y versus X of redundant clusters in the disk %d", i), 400, -20, 20, 400, -20, 20);
    histoClsXYRedundantinDiskNb[i]->SetXTitle("x (cm) position of redundant clusters");
    histoClsXYRedundantinDiskNb[i]->SetYTitle("y (cm) position of redundant clusters");
    histoClsXYRedundantinDiskNb[i]->Sumw2();

    histoClsPhiEtainDiskNb[i] = new TH2D(Form("histoClsPhiEtainDiskNb_%d", i), Form("#eta versus #varphi of clusters in the disk %d", i), 600, -TMath::Pi(), TMath::Pi(), 35, -4.5, -1.);
    histoClsPhiEtainDiskNb[i]->SetXTitle("#varphi of clusters");
    histoClsPhiEtainDiskNb[i]->SetYTitle("#eta of clusters");
    histoClsPhiEtainDiskNb[i]->Sumw2();

    histoClsPhiEtaRedundantinDiskNb[i] = new TH2D(Form("histoClsPhiEtaRedundantinDiskNb_%d", i), Form("#eta versus #varphi of redundant clusters in the disk %d", i), 600, -TMath::Pi(), TMath::Pi(), 35, -4.5, -1.);
    histoClsPhiEtaRedundantinDiskNb[i]->SetXTitle("#varphi of clusters");
    histoClsPhiEtaRedundantinDiskNb[i]->SetYTitle("#eta of clusters");
    histoClsPhiEtaRedundantinDiskNb[i]->Sumw2();
  }

  //Histograms for each MFT layer (layer numbers from 0 to 9)
  for(auto i = 0 ; i < 10 ; i++)
  {
    histoClsXYinLayerNb[i] = new TH2D(Form("histoClsXYinLayerNb_%d", i), Form("Y versus X of clusters in the layer %d", i), 400, -20, 20, 400, -20, 20);
    histoClsXYinLayerNb[i]->SetXTitle("x (cm) position of clusters");
    histoClsXYinLayerNb[i]->SetYTitle("y (cm) position of clusters");
    histoClsXYinLayerNb[i]->Sumw2();

    histoClsPhiEtainLayerNb[i] = new TH2D(Form("histoClsPhiEtainLayerNb_%d", i), Form("#eta versus #varphi of clusters in the layer %d", i), 600, -TMath::Pi(), TMath::Pi(), 35, -4.5, -1.);
    histoClsPhiEtainLayerNb[i]->SetXTitle("#varphi of clusters");
    histoClsPhiEtainLayerNb[i]->SetYTitle("#eta of clusters");
    histoClsPhiEtainLayerNb[i]->Sumw2();

  }

  histoKinPhiEta = new TH2D("histoKinPhiEta", "#eta_{kin} versus #varphi_{kin} of MFT tracks", 600, -TMath::Pi(), TMath::Pi(), 35, -4.5, -1.);
  histoKinPhiEta->SetXTitle("#varphi_{kin}");
  histoKinPhiEta->SetYTitle("#eta_{kin}");
  histoKinPhiEta->Sumw2();

  for(auto i = 5 ; i < 10 ; i++)
  {
    histoKinPhiEta_ncls[i-5] = new TH2D(Form("histoKinPhiEta_ncls_%d", i), Form("#eta_{kin} versus #varphi_{kin} of MFT tracks with ncls #\geq %d", i), 600, -TMath::Pi(), TMath::Pi(), 35, -4.5, -1.);
    histoKinPhiEta_ncls[i-5]->SetXTitle("#varphi_{kin}");
    histoKinPhiEta_ncls[i-5]->SetYTitle("#eta_{kin}");
    histoKinPhiEta_ncls[i-5]->Sumw2();
  }

}
