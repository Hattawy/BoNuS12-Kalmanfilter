void hough_transform3(
    int i_event = 20,
    int Ntracks = 1)
{
  gSystem->Load("libGenFind.so");
  using namespace genfind;
  using namespace ROOT::Math;
  double degree = TMath::Pi()/180.0;

  TFile * f = new TFile("simple_example_out2.root","READ");
  TTree * t = (TTree*)gROOT->FindObject("EVENT");
  if(!t) { std::cout << " Tree not found " << std::endl; return;}
  int nentries=t->GetEntries();

  if(i_event+Ntracks>nentries){
    std::cout << "Only " << nentries << " in file, this parameter set fails" << std::endl;
    return;
  }

  // Setup branches
  std::vector<DD4hep::Simulation::Geant4Tracker::Hit*> * g4hits = nullptr;
  std::vector<DD4hep::Simulation::Geant4Tracker::Hit*> * g4hits2 = nullptr;
  t->SetBranchAddress("SiVertexBarrelHits", &g4hits);
  t->SetBranchAddress("SiTrackerBarrelHits", &g4hits2);

  // Hough Transform class that has all methods. (Would be best to move some member functions to just functions)
  HoughTransform* ht = new HoughTransform();

  std::vector<XYZTVector> hits ;

  // Get the track hits
  for(int i_track = 0; i_track < Ntracks; i_track++){

    t->GetEntry(i_event);
    i_event++;
    
    // SiTrackerBarrelHits
    for( auto thit : (*g4hits2) ) {
      hits.push_back( {thit->position.X(), thit->position.Y(), thit->position.Z(), 0.0} );
      std::cout << thit->position.X() << " , " <<  thit->position.Y() << std::endl;
    }
    // SiVertexBarrelHits
    for( auto thit : (*g4hits) ) {
      hits.push_back( {thit->position.X(), thit->position.Y(), thit->position.Z(), 0.0} );
      std::cout << thit->position.X() << " , " <<  thit->position.Y() << std::endl;
    }

  }


  // ----------------------------------------------

  // A point at the origin
  XYZTVector zero_ref = {0,0,0,0};

  // ----------------------------------------------
  // This returns the hits in conformal space relative to zero_ref=(x0,y0)
  // i.e.: u = (x-x0)/((x-x0)^2+(y-y0)^2), v = -(y-y0)/((x-x0)^2+(y-y0)^2)
  auto chits = ht->GetConformalCoordinates(zero_ref, hits);
  std::vector<XYZTVector> all_chits ;
  all_chits.insert(all_chits.end(), chits.begin(), chits.end());

  // ----------------------------------
  TGraph* gr_UV = new TGraph();
  TGraph* gr_XY = new TGraph();

  // Fill graph with position space XY coordinates
  int i_point = 0;
  for(auto ahit : hits) {
    gr_XY->SetPoint(i_point, ahit.X(), ahit.Y() );
    i_point++;
  }

  // Fill graph with conformal space XY coordinates
  i_point = 0;
  for(auto ahit : all_chits) {
    //std::cout << ahit.X() << " , " <<  ahit.Y() << std::endl;
    gr_UV->SetPoint(i_point, ahit.X(), ahit.Y() );
    i_point++;
  }

  // Fill the histograms for line finding
  auto mg = ht->FillHoughTransforms(all_chits);

  //std::cout << " n peaks : " << ht->FindPeaks(all_chits) << std::endl;

  //auto pretrack_hits = ht->GetPreTrackHits(master_hits,4.0);

  //std::vector<TH1F*> hists;
  //std::vector<TH1F*> alltrack;
  //std::vector<TH1F*> vertextrack;
  //for( auto atrack : pretrack_hits ) {
  //  TH1F* singletrack = (TH1F*)htrack1->Clone();
  //  singletrack->Reset();
  //  singletrack->SetLineColor(2+hists.size());
  //  TH1F* singlevertextrack = (TH1F*)htrack2->Clone();
  //  singlevertextrack->Reset();
  //  singlevertextrack->SetLineColor(2+hists.size());
  //  TH1F* fPhi2 = (TH1F*)ht->fPhi->Clone();
  //  fPhi2->Reset();
  //  fPhi2->SetLineColor(2+hists.size());
  //  for( auto thit : atrack ) {
  //    singletrack->Fill( std::get<0>(thit).X(), std::get<0>(thit).Y());
  //    singlevertextrack->Fill( std::get<0>(thit).X(), std::get<0>(thit).Y());
  //    hpeaky->Fill( std::get<0>(thit).X(), std::get<0>(thit).Y());
  //    hpeaky2->Fill(std::get<0>(thit).X(), std::get<0>(thit).Y());
  //    fPhi2->Fill(TMath::ATan(std::get<1>(thit).Y()/std::get<1>(thit).X())/degree);
  //  }
  //  hists.push_back(fPhi2);
  //  alltrack.push_back(singletrack);
  //  vertextrack.push_back(singlevertextrack);
  //}
//  auto trackone = pretrack_hits.at(0);
//  for( auto thit : trackone ) {
//    htrack1->Fill( std::get<0>(thit).X(), std::get<0>(thit).Y());
//  }
//  auto tracktwo = pretrack_hits.at(1);
//  for( auto thit : tracktwo ) {
//    htrack2->Fill( std::get<0>(thit).X(), std::get<0>(thit).Y());
//  }

  //for(auto t: ht->fHTRoots) {
  //  htheta->Fill(std::get<2>(t));
  //}

  TCanvas * c = new TCanvas();

  gr_XY->SetMarkerStyle(20);
  gr_UV->SetMarkerStyle(20);
  gr_XY->SetMarkerSize(0.5);
  gr_UV->SetMarkerSize(0.5);

  c->Divide(2,2);
  c->cd(1);
  gr_XY->Draw("ap");
  //for(auto ahist : alltrack){
  //  ahist->Draw("box,same");
  //}
  //htrack1->SetLineColor(2);
  //htrack1->Draw("box,same");
  //htrack2->SetLineColor(3);
  //htrack2->Draw("box,same");
  //hpeaky->Draw("box,same");
  c->cd(2);
  gr_UV->Draw("ap");
  //c->cd(3);
  //hxy2->Draw("box");
  //for(auto ahist : vertextrack){
  //  ahist->Draw("box,same");
  //}
  //hpeaky2->SetLineColor(2);
  //hpeaky2->Draw("box,same");
  c->cd(4);
  //hphi->Draw();
  //hrphi->Draw("lego2");

  //c = new TCanvas();
  //c->Divide(2,2);
  //c->cd(1);
  //hxy->Draw("box");
  //c->cd(2);
  //hxz->Draw("box");
  //c->cd(3);
  //hyz->Draw("box");
  //

  //c = new TCanvas();
  //c->Divide(2,2);
  //c->cd(1);
  //ht->fPhi->Draw();
  //for(auto ahist : hists){
  //ahist->Draw("same");
  //}
  ////ht->fRhoTheta0->Draw("colz");
  //c->cd(2);
  //htheta->Draw();

  //c->cd(3);
  //mg->Draw("a");
  //mg->GetYaxis()->SetRangeUser(-0.10,0.10);
  //mg->Draw("a");

  c = new TCanvas();
  c->Divide(2,2);
  c->cd(1);
  ht->fRhoTheta0->Draw("colz");
  c->cd(2);
  ht->fRhoTheta1->Draw("colz");
  c->cd(3);
  ht->fRhoTheta2->Draw("colz");

}
