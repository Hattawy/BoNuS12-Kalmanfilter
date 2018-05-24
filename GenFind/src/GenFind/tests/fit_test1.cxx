#include "ConstField.h"
#include "Exception.h"
#include "FieldManager.h"
#include "KalmanFitterRefTrack.h"
#include "DAF.h"
#include "StateOnPlane.h"
#include "Track.h"
#include "TrackPoint.h"

#include <MaterialEffects.h>
#include <RKTrackRep.h>
#include <TGeoMaterialInterface.h>

#include <EventDisplay.h>

#include <HelixTrackModel.h>
#include <MeasurementCreator.h>
#include <WireMeasurement.h>
#include "PlanarMeasurement.h"

#include <TDatabasePDG.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TRandom.h>
#include <TVector3.h>
#include <vector>

#include "TDatabasePDG.h"
#include <TMath.h>

void fit_test1(
    int i_event = 26,
    int Ntracks = 1)
{

  using namespace genfind;
  using namespace ROOT::Math;

  TFile * f = new TFile("simple_example_out.root","READ");
  TTree * t = (TTree*)gROOT->FindObject("EVENT");
  if(!t) { std::cout << " Tree not found " << std::endl; return;}

  std::vector<DD4hep::Simulation::Geant4Particle*> * mc_particles = nullptr;
  std::vector<DD4hep::Simulation::Geant4Tracker::Hit*> * g4hits = nullptr;
  std::vector<DD4hep::Simulation::Geant4Tracker::Hit*> * g4hits2 = nullptr;
  t->SetBranchAddress("SiVertexBarrelHits", &g4hits);
  t->SetBranchAddress("SiTrackerBarrelHits", &g4hits2);
  t->SetBranchAddress("MCParticles", &mc_particles);

  // Load the geometry
  DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();
  lcdd.fromCompact("JLEIC.xml");
  const double position[3]     = {0,0,0}; // position to calculate magnetic field at (the origin in this case)
  double bField[3]             = {0,0,0};
  lcdd.field().magneticField(position,bField);
  double Bz                    = bField[2]/dd4hep::tesla;
  std::cout << " Magnetic Field Bz             = " << Bz << std::endl;

  //auto& id_decoder = DD4hep::DDRec::IDDecoder::getInstance();
  auto vol_man = lcdd.volumeManager();

  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0., Bz*10.0)); // gentfit uses kilo-Gauss
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();
  display->setOptions("X");

  // Add the geometry to eve display
  TGeoNode*       node1 = gGeoManager->GetTopNode();
  TEveGeoTopNode* its   = new TEveGeoTopNode(gGeoManager, node1);
  gEve->AddGlobalElement(its);

  //------------------------------------------------
  //
  double degree = TMath::Pi()/180.0;
  const int pdg = 11; // particle pdg code

  //genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();
  genfit::AbsKalmanFitter* fitter = new genfit::DAF();
  
  // helix track model
  //const double charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/(3.);
  //genfit::HelixTrackModel* helix = new genfit::HelixTrackModel(pos, mom, charge);

  // trackrep
  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

  // smeared start state
  genfit::MeasuredStateOnPlane stateSmeared(rep);

  const int detId(0); // detector ID
  int planeId(0); // detector plane ID
  int hitId(0); // hit ID

  double detectorResolution(0.0001); // resolution of planar detectors
  TMatrixDSym hitCov(2);
  hitCov.UnitMatrix();
  hitCov *= detectorResolution*detectorResolution;
  TVectorD hitCoords(2);

  for(int i_track = 0; i_track < Ntracks; i_track++) {

    t->GetEntry(i_event);
    i_event++;

    if( (*g4hits).size() == 0 ) {
      continue;
    }
    auto first_hit_pos  = (*((*g4hits)[0])).position;

    // approximate covariance
    TMatrixDSym covM(6);
    double res1 = 0.01;
    for (int i = 0; i < 3; ++i)
      covM(i,i) = res1*res1;
    for (int i = 3; i < 6; ++i)
      covM(i,i) = pow(res1/9.0/sqrt(3.0), 2);


    // start values for the fit, e.g. from pattern recognition
    TVector3 pos(0.0,0.0,0.0);
    pos = 0.1*pos;
    TVector3 mom(first_hit_pos.X(), first_hit_pos.Y(), first_hit_pos.Z());
    mom.SetMag(2.0);

    stateSmeared.setPosMomCov(pos, mom, covM);

    // create track
    TVectorD    seedState(6);
    TMatrixDSym seedCov(6);
    stateSmeared.get6DStateCov(seedState, seedCov);
    genfit::Track fitTrack(rep, seedState, seedCov);


    for( auto thit : (*g4hits) ) {

      TVector3 point     = {thit->position.X()*0.1, thit->position.Y()*0.1, thit->position.Z()*0.1};
      TVector3 zdir      = {0,0,1};
      TVector3 planeNorm = {thit->position.X()*0.1, thit->position.Y()*0.1,0.0};//point.Unit();//clas12::geo::Convert(svt_geometry.GetChannelNorm(chan)/CLHEP::cm);
      planeNorm = planeNorm.Unit();
      genfit::SharedPlanePtr plane(new genfit::DetPlane(point, zdir.Cross(planeNorm) , zdir) );

      // add some planar hits to track with coordinates
      hitCoords[0] = 0.0;//thit->position.X()*0.1;
      hitCoords[1] = 0.0;//thit->position.Y()*0.1;

      //genfit::AbsMeasurement* scint_measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, 1, n_scint_hits, nullptr);
      //static_cast<genfit::PlanarMeasurement*>(scint_measurement)->setPlane(plane, n_scint_hits);
      //fitTrack.insertPoint(new genfit::TrackPoint(scint_measurement, &fitTrack));

      genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, nullptr);
      static_cast<genfit::PlanarMeasurement*>(measurement)->setPlane(plane, thit->cellID);

      fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));
    }
    for( auto thit : (*g4hits2) ) {

      TVector3 point     = {thit->position.X()*0.1, thit->position.Y()*0.1, thit->position.Z()*0.1};
      TVector3 zdir      = {0,0,1};
      TVector3 planeNorm = {thit->position.X()*0.1, thit->position.Y()*0.1,0.0};//point.Unit();//clas12::geo::Convert(svt_geometry.GetChannelNorm(chan)/CLHEP::cm);
      planeNorm = planeNorm.Unit();
      genfit::SharedPlanePtr plane(new genfit::DetPlane(point, zdir.Cross(planeNorm) , zdir) );
      //genfit::SharedPlanePtr plane(new genfit::DetPlane(point, planeNorm.Cross(zdir) ,  (planeNorm.Cross(zdir)).Cross(zdir) ) );

      // add some planar hits to track with coordinates
      hitCoords[0] = 0.0;//thit->position.X()*0.1;
      hitCoords[1] = 0.0;//thit->position.Y()*0.1;

      //genfit::AbsMeasurement* scint_measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, 1, n_scint_hits, nullptr);
      //static_cast<genfit::PlanarMeasurement*>(scint_measurement)->setPlane(plane, n_scint_hits);
      //fitTrack.insertPoint(new genfit::TrackPoint(scint_measurement, &fitTrack));

      //auto pv = vol_man.lookupPlacement( thit->cellID );

      genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, nullptr);
      measurement->setPlane(plane, thit->cellID);

      fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));
    }
   
    //check
    fitTrack.checkConsistency();

    // do the fit
    fitter->processTrack(&fitTrack);

    // print fit result
    fitTrack.getFittedState().Print();

    //check
    fitTrack.checkConsistency();

    display->addEvent(&fitTrack);

    //genfit::SharedPlanePtr test_plane(new genfit::DetPlane({0,0,0}, {1,0,0} , {0,1,0}) );
    //std::vector<genfit::MeasurementOnPlane*> constructMeasurementsOnPlane  ...
    //delete fitter;
  }

  //    hits.push_back( {thit->position.X(), thit->position.Y(), thit->position.Z(), 0.0} );
  //    std::cout << thit->position.X() << " , " <<  thit->position.Y() << std::endl;
  //    hxy->Fill(thit->position.X(), thit->position.Y());
  //    hxy2->Fill(thit->position.X(), thit->position.Y());
  //  }

  //  for( auto thit : (*g4hits2) ) {
  //    hits.push_back( {thit->position.X(), thit->position.Y(), thit->position.Z(), 0.0} );
  //    std::cout << thit->position.X() << " , " <<  thit->position.Y() << std::endl;
  //    hxy->Fill(thit->position.X(), thit->position.Y());
  //    hxy2->Fill(thit->position.X(), thit->position.Y());
  //  }
  //}


  //hitCoords[0] = -0.15;
  //hitCoords[1] = 0;
  //measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, nullptr);
  //measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,10), TVector3(1,0,0), TVector3(0,1,0))), ++planeId);
  //fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

  //hitCoords[0] = -0.4;
  //hitCoords[1] = 0;
  //measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, nullptr);
  //measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,20), TVector3(1,0,0), TVector3(0,1,0))), ++planeId);
  //fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

  std::cout << " |P| = " << (*mc_particles)[0]->psx << ", " << (*mc_particles)[0]->psy << ", " << (*mc_particles)[0]->psz << std::endl ;


  // open event display
  display->open();

  //// approximate covariance
  //TMatrixDSym covM(6);
  //double resolution = 0.01;
  // for (int i = 0; i < 3; ++i)
  //    covM(i,i) = resolution*resolution;
  // for (int i = 3; i < 6; ++i)
  //    covM(i,i) = pow(resolution/9.0/sqrt(3.0), 2);


  //std::vector<std::tuple<XYZTVector,XYZTVector,XYZTVector>> master_hits ;
  //std::vector<XYZTVector> all_chits ;

  //std::cout << "DERP\n";
  //for(const auto& fhit : hits){
  ////auto fhit = hits.at(0);
  //  auto chits = ht->GetConformalCoordinates(fhit, hits);
  //  all_chits.insert(all_chits.end(), chits.begin(), chits.end());
  //  //std::transform(hits.begin(), hits.end(), chits.begin(), master_hits.end(),
  //  //    [&](const auto& a, const auto& b){ return make_tuple(a,b,fhit);});
  //  for(int ihit = 0; ihit < hits.size() ; ihit++){
  //    master_hits.push_back(make_tuple(hits[ihit], chits[ihit], fhit));
  //  }
  //}
  //std::cout << " Master " << master_hits.size() << std::endl;


  //for(auto ahit : all_chits) {
  //  std::cout << ahit.X() << " , " <<  ahit.Y() << std::endl;
  //  huv->Fill( ahit.X(), ahit.Y() );
  //  hphi->Fill(TMath::ATan(ahit.Y()/ahit.X()));
  //  for(int i_theta = 0; i_theta<180; i_theta++) {
  //    double theta = 2.0*double(i_theta)*degree;
  //    hrphi->Fill(theta/degree, ahit.X()*TMath::Cos(theta) + ahit.Y()*TMath::Sin(theta) );
  //  }
  //}
  //auto mg = ht->FillHoughTransforms(all_chits);

  //std::cout << " n peaks : " << ht->FindPeaks(all_chits) << std::endl;

  //auto pretrack_hits = ht->GetPreTrackHits(master_hits,4.0);


  //std::vector<TH1F*> hists;
  //for( auto atrack : pretrack_hits ) {
  //  TH1F* fPhi2 = (TH1F*)ht->fPhi->Clone();
  //  fPhi2->Reset();
  //  fPhi2->SetLineColor(2+hists.size());
  //  for( auto thit : atrack ) {
  //    hpeaky->Fill( std::get<0>(thit).X(), std::get<0>(thit).Y());
  //    hpeaky2->Fill(std::get<0>(thit).X(), std::get<0>(thit).Y());
  //    fPhi2->Fill(TMath::ATan(std::get<1>(thit).Y()/std::get<1>(thit).X())/degree);
  //  }
  //  hists.push_back(fPhi2);
  //}

  //for(auto t: ht->fHTRoots) {
  //  htheta->Fill(std::get<2>(t));
  //}

  //TCanvas * c = new TCanvas();

  //c->Divide(2,2);
  //c->cd(1);
  //hxy->Draw("box");
  //hpeaky->SetLineColor(2);
  //hpeaky->Draw("box,same");
  //c->cd(2);
  //huv->Draw("box");
  //c->cd(3);
  //hxy2->Draw("box");
  //hpeaky2->SetLineColor(2);
  //hpeaky2->Draw("box,same");
  //c->cd(4);
  //hphi->Draw();
  ////hrphi->Draw("lego2");

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

  //c = new TCanvas();
  //c->Divide(2,2);
  //c->cd(1);
  //ht->fRhoTheta0->Draw("colz");
  //c->cd(2);
  //ht->fRhoTheta1->Draw("colz");
  //c->cd(3);
  //ht->fRhoTheta2->Draw("colz");

}
