void hough_transform(){

  using namespace genfind;
  using namespace ROOT::Math;

  HoughTransform ht;

  TH2F * hxy = new TH2F("hxy", "hxy", 50,-1,1,50,0,1);
  TH2F * huv = new TH2F("huv", "huv", 50,-10,10,50,-1,1);
  TH1F * hphi = new TH1F("huv", "huv", 50,-4,4);

  std::vector<XYZTVector> hits ;

  for(int i_track = 0; i_track < 10; i_track++){
    double R  = gRandom->Uniform();
    double x0 = gRandom->Uniform(-0.5,0.5);
    double dx = gRandom->Uniform(-0.05, 0.05);

    for(int i = 0; i < 20; i++){
      double x = x0 + dx*double(i);
      double y = TMath::Sqrt(R*R-x*x);
      hits.push_back( {x,y,0.0,0.0});
      hxy->Fill(x,y);
    }

    auto chits = ht.GetConformalCoordinates(hits);

    for(auto ahit : chits) {
      huv->Fill( ahit.X(), ahit.Y() );
      hphi->Fill(ahit.Phi());
    }
  }

  TCanvas * c = new TCanvas();

  c->Divide(2,2);
  c->cd(1);
  hxy->Draw("colz");
  c->cd(2);
  huv->Draw("colz");
  c->cd(3);
  hphi->Draw("colz");

}
