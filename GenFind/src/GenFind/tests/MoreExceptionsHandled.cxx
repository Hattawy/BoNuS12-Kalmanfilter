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

void MoreExceptionsHandled(
    int i_event = 20,
    int Ntracks = 1)
{

  using namespace genfind;
  using namespace ROOT::Math;

  int increment_event = i_event;

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

  genfit::AbsKalmanFitter* fitter = new genfit::DAF();
  
  // trackrep
  //genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

  // smeared start state
  //genfit::MeasuredStateOnPlane stateSmeared(rep);

  // approximate covariance
  TMatrixDSym covM(6);
  double res1 = 0.01;
  for (int i = 0; i < 3; ++i)
    covM(i,i) = res1*res1;
  for (int i = 3; i < 6; ++i)
    covM(i,i) = pow(res1/9.0/sqrt(3.0), 2);


  const int detId(0); // detector ID
  int planeId(0); // detector plane ID
  int hitId(0); // hit ID

  double detectorResolution(0.0001); // resolution of planar detectors
  TMatrixDSym hitCov(2);
  hitCov.UnitMatrix();
  hitCov *= detectorResolution*detectorResolution;
  TVectorD hitCoords(2);
  
 //----------------------------------------------------------------------------------------------
  
  std::vector<XYZTVector> hits ;
  std::vector<int> cellID_hits ;
  std::vector<int>  tracktruth ;
  increment_event = i_event;

  for(int i_track = 0; i_track < Ntracks; i_track++){ // made tracks by combining events
 
    t->GetEntry(increment_event);

    for( auto thit : (*g4hits) ) {
      hits.push_back( {thit->position.X(), thit->position.Y(), thit->position.Z(), 0.0} );
      cellID_hits.push_back( thit->cellID );
      tracktruth.push_back(i_track);
      std::cout << thit->position.X() << " , " <<  thit->position.Y() << std::endl;
    }

    for( auto thit : (*g4hits2) ) {
      hits.push_back( {thit->position.X(), thit->position.Y(), thit->position.Z(), 0.0} );
      cellID_hits.push_back( thit->cellID );
      tracktruth.push_back(i_track);
      std::cout << thit->position.X() << " , " <<  thit->position.Y() << std::endl;
    }
    std::cout << "added one original event to combined 'event' " << std::endl;
    increment_event++;
  }
  std::cout << "Created 'event' by combining " << Ntracks << " real events from file\n";

  HoughTransform* ht = new HoughTransform();
  // run the standard GenFind method over composit event to disaggregate via algorithm
  std::vector<ROOT::Math::XYZTVector> conf_hits = ht->UseConformalMap(hits);

  auto mg = ht->FillHoughTransforms(conf_hits);
  std::vector<std::vector<std::tuple<int,double,double>>> FoundTracks = ht->FindTracks(conf_hits);
  std::cout << "FindTracks found " << FoundTracks.size() << " tracks " << endl;

 //-------------------------------------------------------------------------------------------
 
  bool exceptionbool=false;

  std::vector<genfit::Track*> all_tracks; 
  for( std::vector<std::tuple<int,double,double>> atrack : FoundTracks ) {

    t->GetEntry(i_event); // might have problems with multiple tracks this way?
 
    if( (*g4hits).size() == 0 ) {
      continue;
    }
    // get first hit position
    double firstx = hits[std::get<0>(atrack[0])].X();
    double firsty = hits[std::get<0>(atrack[0])].Y();
    double firstz = hits[std::get<0>(atrack[0])].Z();
    TVector3 first_hit_pos(firstx,firsty,firstz);
    //auto first_hit_pos  = (*((*g4hits)[0])).position;
    std::cout << "first hit pos x " << first_hit_pos.X() << std::endl;
    std::cout << "first hit pos y " << first_hit_pos.Y() << std::endl;
    std::cout << "first hit pos z " << first_hit_pos.Z() << std::endl;

    // trackrep
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

    // smeared start state
    genfit::MeasuredStateOnPlane stateSmeared(rep);

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
    //genfit::Track fitTrack(rep, seedState, seedCov);
    auto aTrack = new genfit::Track(rep, seedState, seedCov);
    auto& fitTrack = *aTrack;
    all_tracks.push_back(aTrack);

    for( auto thit : atrack ) {
      std::cout << " did a hit! " << std::endl;
      std::cout << "x position is " << hits[std::get<0>(thit)].X() << std::endl;
      std::cout << "y position is " << hits[std::get<0>(thit)].Y() << std::endl;
      std::cout << "z position is " << hits[std::get<0>(thit)].Z() << std::endl;
      TVector3 point     = {hits[std::get<0>(thit)].X()*0.1, hits[std::get<0>(thit)].Y()*0.1, hits[std::get<0>(thit)].Z()*0.1};
      TVector3 zdir      = {0,0,1};
      TVector3 planeNorm = {hits[std::get<0>(thit)].X()*0.1, hits[std::get<0>(thit)].Y()*0.1,0.0};
      planeNorm = planeNorm.Unit();
      genfit::SharedPlanePtr plane(new genfit::DetPlane(point, zdir.Cross(planeNorm) , zdir) );

      // add some planar hits to track with coordinates
      hitCoords[0] = 0.0;//thit->position.X()*0.1;
      hitCoords[1] = 0.0;//thit->position.Y()*0.1;


      genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, nullptr);
      measurement->setPlane(plane, cellID_hits[std::get<0>(thit)]);
      std::cout << "cellID is " << cellID_hits[std::get<0>(thit)] << std::endl;

      fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));
   


    }     
  
    //check
    try {
      fitTrack.checkConsistency();
    } catch (...) {
      std::cout << "These statements never run regardless, dunno why." << std::endl;
      exceptionbool=true;
    } 

    // do the fit
    try {
      //these must be together because the second one can avoid a segfault from just the first
      //but it does create an exception that must be handled or the script exits.
      std::cout << "doing the fit" << std::endl;
      fitter->processTrack(&fitTrack);
      std::cout << "checking fit results: " << std::endl;
      fitTrack.checkConsistency();
    } catch (...) {
      std::cout <<  "I guess I threw an exception" << std::endl;
      exceptionbool=true;
      break;
    }
    

    // print fit result
    std::cout << " printing fit results: " << std::endl;
    //fitTrack.getFittedState().Print();

    //check
    //std::cout << " checking fit results: " << std::endl;
    //fitTrack.checkConsistency();

    //display->addEvent(&fitTrack);
  
  
  }

    std::cout << " |P| = " << (*mc_particles)[0]->psx << ", " << (*mc_particles)[0]->psy << ", " << (*mc_particles)[0]->psz << std::endl ;
  
    display->addEvent(all_tracks);
    // open event display
    display->open();

}
