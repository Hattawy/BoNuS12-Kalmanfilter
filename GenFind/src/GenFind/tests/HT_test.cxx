void HT_test(
    int i_event = 20,
    int Ntracks = 3)
{


// test the efficiency of separating hits by track. when reconstructing, also fill a "truth" vector
// which holds which track each hit is supposed to go to. Compare against the result, check performance
//
  
  gSystem->Load("libGenFind.so");
  using namespace genfind;
  using namespace ROOT::Math;

  TFile * f = new TFile("simple_example_out.root","READ");
  TTree * t = (TTree*)gROOT->FindObject("EVENT");
  if(!t) { std::cout << " Tree not found " << std::endl; return;}
  int nentries=t->GetEntries();

  if(i_event+Ntracks>nentries){
    std::cout << "Only " << nentries << " in file, this parameter set fails" << std::endl;
    return;
  }

  std::vector<DD4hep::Simulation::Geant4Tracker::Hit*> * g4hits = nullptr;
  std::vector<DD4hep::Simulation::Geant4Tracker::Hit*> * g4hits2 = nullptr;
  t->SetBranchAddress("SiVertexBarrelHits", &g4hits);
  t->SetBranchAddress("SiTrackerBarrelHits", &g4hits2);

  HoughTransform* ht = new HoughTransform();

  TH2F * hxy   = new TH2F("hxy",   "hxy",   100, -200, 200, 100, -200, 200);
  TH2F * hxz   = new TH2F("hxz",   "hxz",   100, -1000, 1000, 100, -1000, 1000);
  TH2F * hyz   = new TH2F("hyz",   "hyz",   100, -1000, 1000, 100, -1000, 1000);
  TH2F * hxy2  = new TH2F("hxy2",  "hxy2",  100, -1010 ,1010, 100, -1010, 1010);
  TH2F * huv   = new TH2F("huv",   "huv",   100, -120 , 120, 100, -120, 120);



  //make composite tracks with real datapoints
  std::vector<XYZTVector> hits ;
  std::vector<int>  tracktruth ;
  double degree = TMath::Pi()/180.0;
  for(int i_track = 0; i_track < Ntracks; i_track++){ // made tracks by combining events

    t->GetEntry(i_event);
    i_event++; // this causes an error, because it increments before you do anything to the track
    
    for( auto thit : (*g4hits) ) {
      hits.push_back( {thit->position.X(), thit->position.Y(), thit->position.Z(), 0.0} );
      tracktruth.push_back(i_track);
      std::cout << thit->position.X() << " , " <<  thit->position.Y() << std::endl;
      hxy->Fill(thit->position.X(), thit->position.Y());
      hxz->Fill(thit->position.X(), thit->position.Z());
      hyz->Fill(thit->position.Y(), thit->position.Z());
      hxy2->Fill(thit->position.X(), thit->position.Y());
    }

    for( auto thit : (*g4hits2) ) {
      hits.push_back( {thit->position.X(), thit->position.Y(), thit->position.Z(), 0.0} );
      tracktruth.push_back(i_track);
      std::cout << thit->position.X() << " , " <<  thit->position.Y() << std::endl;
      hxy->Fill(thit->position.X(), thit->position.Y());
      hxz->Fill(thit->position.X(), thit->position.Z());
      hyz->Fill(thit->position.Y(), thit->position.Z());
      hxy2->Fill(thit->position.X(), thit->position.Y());
    }
  }
  std::cout << "Created 'event' by combining " << Ntracks << " real events from file\n";
  
  
  //make track with artificial linear points.
  /*for(int i_hit = 0; i_hit < 100; i_hit++){  // make test hits in straight line
    double x_incr = 1;
    double y_incr = 1;
    double z_incr = 1;
    hits.push_back( {x_incr*i_hit+10 , y_incr*i_hit, z_incr*i_hit , 0.0} );
    std::cout << x_incr*i_hit << " , " <<  y_incr*i_hit << std::endl;
    hxy->Fill(x_incr*i_hit+10, y_incr*i_hit);
    hxz->Fill(x_incr*i_hit+10, z_incr*i_hit);
    hyz->Fill(y_incr*i_hit, z_incr*i_hit);
    hxy2->Fill(x_incr*i_hit+10, y_incr*i_hit);

  }

  double x_extra = 0.0;
  double y_extra = 50.0;
  double z_extra = 0.0;
  
  hits.push_back( {x_extra, y_extra, z_extra, 0.0} );
  std::cout << x_extra << " , " <<  y_extra << std::endl;
  hxy->Fill(x_extra, y_extra);
  hxz->Fill(x_extra, z_extra);
  hyz->Fill(y_extra, z_extra);
  hxy2->Fill(x_extra, y_extra);*/


  //  auto mg = ht->FillHoughTransforms(hits);


  // use track truth and the tracknumber in the findtracks result to test performance.
  // First, are all the hits properly separated,
  // Second, are all the hits accounted for
  
  std::vector<ROOT::Math::XYZTVector> conf_hits = ht->UseConformalMap(hits);

  for (auto conf_hit : conf_hits ){
    huv->Fill(conf_hit.X(),conf_hit.Y()) ; 
  
  }

  auto mg = ht->FillHoughTransforms(conf_hits);
  std::vector<std::vector<std::tuple<int,double,double>>> FoundTracks = ht->FindTracks(conf_hits);
  //std::vector<std::vector<std::tuple<int,double,double>>> FoundTracks = ht->FindTracks(hits);
  std::cout << "FindTracks found " << FoundTracks.size() << " tracks " << endl;
  std::vector<TH1F*> FoundTrackHists;
  int tracknumber=0;
  for( std::vector<std::tuple<int,double,double>> atrack : FoundTracks ) {
    tracknumber++;
    TH1F* singletrack = (TH1F*)hxy2->Clone();
    singletrack->Reset();
    singletrack->SetLineColor(2+FoundTrackHists.size());

    int ref_truth_track = tracktruth[std::get<0>(atrack[0])];
    int ref_diff = tracknumber - ref_truth_track;
    //std::cout << "ref diff is " << ref_diff << std::endl;
    double grouped_test = 0;
    int hit_diff=0;
    int numhit=0;
    for( std::tuple<int,double,double> thit : atrack ) {
         hit_diff = hit_diff + tracknumber - tracktruth[std::get<0>(thit)]-ref_diff;
         numhit++; 
         //std::cout << "track number " << tracknumber << " hit number " << std::get<0>(thit) << std::endl;
         //std::cout << "track truth " << tracktruth[std::get<0>(thit)]; 
         //std::cout << " x value " << std::get<1>(thit) << " y value " << std::get<2>(thit) << std::endl;
         //singletrack->Fill( std::get<1>(thit), std::get<2>(thit) );
         singletrack->Fill( hits[std::get<0>(thit)].X(), hits[std::get<0>(thit)].Y() );
    }
    if (numhit>0) {
      grouped_test = hit_diff/numhit;
    } else {
      grouped_test = -1;
    }
    std::cout << "grouped test " << grouped_test << std::endl;
    FoundTrackHists.push_back(singletrack);
  }

// now, check for missed hits, hits that belonged to a track but were left out instead.
//
//  doubled_hits = 0;  // more complicated to check this, don't right now.
  int missed_hits = 0;
  int total_hits =0;
  for( int hit_index=0; hit_index < tracktruth.size(); hit_index++ ) {
    total_hits++;
    //std::cout << " hit index " << hit_index << std::endl;
    bool hit_listed = false;
    for( std::vector<std::tuple<int,double,double>> atrack : FoundTracks ) {
      for ( std::tuple<int,double,double> thit : atrack ) {
        if (hit_index == std::get<0>(thit)) {
          hit_listed = true;
        }
      }   
    }
    if (!hit_listed)  {
      missed_hits++;
    }
  } 
  std::cout << "missed hits " << missed_hits << " total hits " << total_hits << std::endl;

  TCanvas * c = new TCanvas();

  c->Divide(2,2);
  c->cd(1);
  mg->Draw("a");
  std::cout << "multigraph list is " << mg->GetListOfGraphs()->GetSize() << endl;
  std::cout << mg->GetListOfGraphs() << endl;
  std::cout << mg->GetListOfFunctions() << endl;
  c->cd(2);
  mg->GetYaxis()->SetRangeUser(-50,50);
  mg->Draw("a");
  c->cd(3);
  huv->Draw();
  c->cd(4);
  ht->fIntersect->Draw();

  c = new TCanvas();
  c->Divide(2,2);
  c->cd(1);
  hxy2->Draw("box");
  c->cd(2);
  if(FoundTrackHists.size()>0) {FoundTrackHists[0]->Draw("box");}
  c->cd(3);
  if(FoundTrackHists.size()>1) {FoundTrackHists[1]->Draw("box");}
  c->cd(4);
  if(FoundTrackHists.size()>2) {FoundTrackHists[2]->Draw("box");}


}
