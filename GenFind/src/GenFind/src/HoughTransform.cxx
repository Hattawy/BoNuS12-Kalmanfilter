#include "HoughTransform.h"
#include "TMath.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include <algorithm>
#include "TRandom3.h"

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "TGraph.h"
#include "TList.h"
#include "TObject.h"

namespace genfind {

  HoughTransform::HoughTransform()
  {
    bool status = TH1::AddDirectoryStatus();
    TH1::AddDirectory(false);
    fPhi       = new TH1F{"hPhi","hPhi",90,-90, 90};
    fRhoTheta0 = new TH2F{"hRhoTheta0","hRhoTheta0",100,-1010, 1010, 100,-1010,1010};
    fRhoTheta1 = new TH2F{"hRhoTheta1","hRhoTheta1",90, 0, 180, 100,-0.1,0.1};
    fRhoTheta2 = new TH2F{"hRhoTheta2","hRhoTheta2",40, 0, 180, 10,-0.1,0.1};
    fIntersect = new TH2F{"hIntersect","hIntersect",90,0, 180, 75,-150,150 };
    fIntReduced= new TH2F{"hIntReduced","hIntReduced",90,0, 180, 75,-150,150 };
    fIntLeft   = new TH2F{"hIntLeft","hIntLeft",90,0, 180, 75,-150,150 };
    fSpec  = new TSpectrum2(2*fMaxTracks);
    fSpec1 = new TSpectrum(2*fMaxTracks);
    TH1::AddDirectory(status);
  }
  //__________________________________________________________________________

  HoughTransform::~HoughTransform()
  {
    delete fPhi;
    delete fRhoTheta0;
    delete fRhoTheta1;
    delete fRhoTheta2;
    delete fIntersect;
    if(fSpec) delete fSpec; fSpec = nullptr;
    if(fSpec1) delete fSpec1; fSpec1 = nullptr;
  }
  //__________________________________________________________________________
  //
     // the conformal map, 1/z, for complex value z, transforms a circle into a line
     // However, it only does this for a circle positioned on (radius,0).
     // The method finds the center and radius for a hit and a reference hit,
     // assuming the track segment comes from the vertex
     //
     // the conformal map is only saved if there exists a fourth hit, bhit,
     // which also lies on the defined circle.
     //
     // If the points are such that the radius is too large, the points are checked for linearity.

  std::vector<ROOT::Math::XYZTVector> HoughTransform::GetConformalCoordinates(
      const ROOT::Math::XYZTVector ref_point,
      const std::vector<ROOT::Math::XYZTVector>& hits) const
  {
    std::vector<ROOT::Math::XYZTVector> res;
    //std::for_each(hits.begin(), hits.end(), [&](const auto& ahit) {
    
    for(const auto& ahit : hits) {
      for(const auto & bhit : hits) {
        // we are assuming an x0,y0 that is the vertex.
        double x1 = ahit.X();
        double x2 = ref_point.X();
        double x3 = bhit.X();
        double y1 = ahit.Y(); 
        double y2 = ref_point.Y();
        double y3 = bhit.Y();

        double x = x1-x2;
        double y = y1-y2;
        double R = x*x+y*y;

        if( !(ahit == ref_point) ) {
          res.push_back(ROOT::Math::XYZTVector{ 2.0*x/R, 2.0*y/R, ahit.Z(), ahit.T() });
        }
      }
    };
    return res;
  }
  //_________________________________________________________________________
  //
  
  //std::vector<double> HoughTransform::UseConformalMap(

  std::vector<ROOT::Math::XYZTVector> HoughTransform::UseConformalMap(
     const std::vector<ROOT::Math::XYZTVector>& hits) const
   {
      std::vector<ROOT::Math::XYZTVector> res;
      std::for_each(hits.begin(), hits.end(), [&](const auto& ahit) {
         double x = ahit.X();
         double y = ahit.Y();
         double R = x*x + y*y;
         res.push_back(ROOT::Math::XYZTVector{ 2000.0*x/R, -2000.0*y/R, ahit.Z(), ahit.T() });
      });
     return res;    
  }
  //_________________________________________________________________________
  //
// merged in, maybe not used
//  void HoughTransform::FillHistograms(const std::vector<genfind::ConformalHit>& chits)
//  {
//    double degree = TMath::Pi()/180.0;
//    double phi_min = 0.0;
//    double phi_max = 180.0;
//
//    TH1F* temp0 = (TH1F*)fRhoTheta0->Clone();
//    TH1F* temp1 = (TH1F*)fRhoTheta1->Clone(); 
//    TH1F* temp2 = (TH1F*)fRhoTheta2->Clone(); 
//    fRhoTheta0->Reset();
//    fRhoTheta1->Reset();
//    fRhoTheta2->Reset();
//    for(const auto& ahit : chits) {
//      temp0->Reset();
//      temp1->Reset();
//      temp2->Reset();
//      for(int i_theta = 0; i_theta<1000; i_theta++) {
//
//        double theta = 180.0*degree*double(i_theta)/double(500);//rand.Uniform(0.0,180.0)*degree;
//
//        int bin = temp0->FindBin(theta/degree, ahit.X()*TMath::Cos(theta) + ahit.Y()*TMath::Sin(theta) );
//        temp0->SetBinContent(bin, 1);
//
//        bin = temp1->FindBin(theta/degree, ahit.X()*TMath::Cos(theta) + ahit.Y()*TMath::Sin(theta) );
//        temp1->SetBinContent(bin, 1);
//
//        bin = temp2->FindBin(theta/degree, ahit.X()*TMath::Cos(theta) + ahit.Y()*TMath::Sin(theta) );
//        temp2->SetBinContent(bin, 1);
//      }
//      fRhoTheta0->Add(temp0);
//      fRhoTheta1->Add(temp1);
//      fRhoTheta2->Add(temp2);
//    }
//
//    delete temp0;
//    delete temp1;
//    delete temp2;
//  }
//  //__________________________________________________________________________
//
//  std::vector<double> HoughTransform::FindPhiPeaks(const std::vector<genfind::ConformalHit>& chits)
//  {
//    double degree = TMath::Pi()/180.0;
//    // First fill the phi histogram
//    fPhi->Reset();
//    for(const auto& ahit : chits){
//      double phi  = TMath::ATan(ahit.Y()/ahit.X())/degree;
//      fPhi->Fill(phi);
//    }
//    double sigma  = 2;
//    double thresh = 0.1;
//    int    npeaks = fSpec1->Search(fPhi, sigma, "nodraw", thresh );
//    std::vector<double> peaks;
//    for(int i = 0; i<npeaks; i++) {
//      peaks.push_back( (fSpec1->GetPositionX())[i] );
//    }
//    return peaks;
//  }
//  //__________________________________________________________________________
//
  TMultiGraph* HoughTransform::FillHoughTransforms(const std::vector<ROOT::Math::XYZTVector>& hits)
  {
    double degree = TMath::Pi()/180.0;
    double phi_min = 0.0;
    double phi_max = 180.0;

    // This part is not really needed
    TMultiGraph* mg = new TMultiGraph();
    //std::vector<TF1*> funcs;
    double avg = 0;

    
    for(const auto& ahit : hits) {
      TF1 * f = new TF1("f",[&](double*x, double *p){ return p[0]*TMath::Cos(x[0]*degree) + p[1]*TMath::Sin(x[0]*degree); }, phi_min, phi_max, 2);
      f->SetLineColor(1);
      f->SetLineWidth(1);
      f->SetParameter(0, ahit.X());
      f->SetParameter(1, ahit.Y());
      fFuncs.push_back(f);
      mg->Add(new TGraph(f),"l");

    }
    return mg;
    //fFuncs as member changes in ways I don't understand.
    //Not constant, can't be trusted in other class methods.
  }
//_______________________________________________________________________________

  std::vector<TF1*> HoughTransform::FillFunctions(const std::vector<ROOT::Math::XYZTVector>& hits)
  {
    double degree = TMath::Pi()/180.0;
    double phi_min = 0.0;
    double phi_max = 180.0;
    std::vector<TF1*> funcs;
    double avg = 0;

    
    for(const auto& ahit : hits) {
      TF1 * f = new TF1("f",[&](double*x, double *p){ return p[0]*TMath::Cos(x[0]*degree) + p[1]*TMath::Sin(x[0]*degree); }, phi_min, phi_max, 2);
      f->SetParameter(0, ahit.X());
      f->SetParameter(1, ahit.Y());
      funcs.push_back(f);
    }
    return funcs;

    // for some reason, the returned vector from this is not constant. 
    // Changes in ways I don't understand.
  }
//______________________________________________________________________________
//
    // this function finds the crossings for every pair of hit-functions.
    // These are only kept in the list if a third hit-function passes very near.
 

  std::vector<std::tuple<int,int,double,double>> HoughTransform::FindIntersections( const std::vector<ROOT::Math::XYZTVector>& hits ) {

    //create a list of rho/theta functions from the hit parameters
    double degree = TMath::Pi()/180.0;
    double phi_min = 0.0;
    double phi_max = 180.0;
    std::vector<TF1*> funcs;
    std::vector<std::tuple<int,int,double,double>> Intersections;

    for(const auto& ahit : hits) {
      TF1 * f = new TF1("f",[&](double*x, double *p){ return p[0]*TMath::Cos(x[0]*degree) + p[1]*TMath::Sin(x[0]*degree); }, phi_min, phi_max, 2);
      f->SetParameter(0, ahit.X());
      f->SetParameter(1, ahit.Y());
      funcs.push_back(f);
    }
    int n_funcs = funcs.size();


    //Find the crossings for every pair, by solving roots of difference functions.
    int n_roots = 0;
    for(int i_func = 0; i_func < n_funcs; i_func++) {
      for(int j_func = i_func+1; j_func < n_funcs; j_func++) {
        auto f1 = funcs[i_func];
        auto f2 = funcs[j_func];
        TF1 * fdiff = new TF1("fdiff",[&](double*x, double *p){
            double val1 = f1->Eval(x[0]);
            double val2 = f2->Eval(x[0]);
            double dif = val1-val2;
            return dif;
            }, phi_min, phi_max, 0);
        ROOT::Math::WrappedTF1 wf1(*fdiff);
        ROOT::Math::BrentRootFinder brf;
        brf.SetFunction( wf1, phi_min, phi_max);
        brf.Solve();
        double root = brf.Root();
        fHTRoots.push_back( std::make_tuple(i_func, j_func, root) );
        n_roots++;

        double rhoroot = f1->Eval(root);

        //For better effect, make sure that there is a triple coincidence
        //fThresh = 3;        // this works well for regular, not conformal hits
        fThresh = 1;
        for(int k_func = 0; k_func < n_funcs; k_func++) {
          if (k_func==i_func || k_func==j_func){continue;}
          auto f3 = funcs[k_func];
          double checkval = f3->Eval(root);
          if (abs(checkval-rhoroot)<fThresh) {
            fIntersect->Fill(root,rhoroot);
            fHTIntersections.push_back( std::make_tuple(i_func, j_func, root, rhoroot) );
            Intersections.push_back( std::make_tuple(i_func, j_func, root, rhoroot) );
            break;
          }
        }
      }
    }
    return Intersections;
  }

//_______________________________________________________________________________
//
    // This function selects one reference hit, using a histogram
    // This hit is in the bin with the most intersections
    // It isn't actually a needed function if FindTrack just uses the intersection method itself...
 
  int HoughTransform::GetReferenceHit( const std::vector<ROOT::Math::XYZTVector>& hits, std::vector<std::tuple<int,int,double,double>> Intersections )  {


    //recreate the damn funcs vector since it's not stable otherwise...
    double degree = TMath::Pi()/180.0;
    double phi_min = 0.0;
    double phi_max = 180.0;
    std::vector<TF1*> funcs;
    for(const auto& ahit : hits) {
      TF1 * f = new TF1("f",[&](double*x, double *p){ return p[0]*TMath::Cos(x[0]*degree) + p[1]*TMath::Sin(x[0]*degree); }, phi_min, phi_max, 2);
      f->SetParameter(0, ahit.X());
      f->SetParameter(1, ahit.Y());
      funcs.push_back(f);
    }
    int n_funcs = funcs.size();

    //find the location of the most intersections, using hist bins.
    TH2F* hist_intersect = new TH2F{"hist_intersect","hist_intersect",90,0, 180, 75,-150,150 };
    for(auto t : Intersections ){
      hist_intersect->Fill(std::get<2>(t),std::get<3>(t));
    }
    double max = hist_intersect->GetMaximumBin();
    double value = hist_intersect->GetBinContent(max);
    Int_t binx,biny,binz;
    hist_intersect->GetBinXYZ(max, binx, biny, binz); 
    double xmaxbin = hist_intersect->GetXaxis()->GetBinCenter(binx);
    double ymaxbin = hist_intersect->GetYaxis()->GetBinCenter(biny);
    std::cout << "binx " << xmaxbin << "biny " << ymaxbin << std::endl;
    delete hist_intersect;


    //find one hit with function that goes closest to the maximum bin (
    int min = abs(funcs[0]->Eval(xmaxbin)-ymaxbin);
    int index = -1;
    for(int i = 0; i < n_funcs; i++) {
      auto f1 = funcs[i];
      double val = f1->Eval(xmaxbin);
      double dif = ymaxbin-val;
      std::cout << "dif is " << dif << std::endl;
  
      if (abs(dif) < min) {
        min = abs(dif);
        index = i;        
      }
    }
    
    std::cout << "reference " << index << std::endl;
    return index;
  }

//_____________________________________________________________________________
//
    //This finds a list of candidate hits constituting a single track.
    //It requires a list of hits and intersections, created by HoughTransform::FindIntersections
    //It also requires a reference hit. Suitable hit found by HoughTransform::GetReferenceHit
    //It doesn't actually need to require a reference hit, just does right now.
    //The track is the longest candidate track in the list of hit intersections given
    
  std::vector< int > HoughTransform::FindTrack( const std::vector<ROOT::Math::XYZTVector>& hits, std::vector<std::tuple<int,int,double,double>> Intersections, int ref ) {


    //recreate the damn funcs vector since it's not stable otherwise...
    double degree = TMath::Pi()/180.0;
    double phi_min = 0.0;
    double phi_max = 180.0;
    std::vector<TF1*> funcs;
    for(const auto& ahit : hits) {
      TF1 * f = new TF1("f",[&](double*x, double *p){ return p[0]*TMath::Cos(x[0]*degree) + p[1]*TMath::Sin(x[0]*degree); }, phi_min, phi_max, 2);
      f->SetParameter(0, ahit.X());
      f->SetParameter(1, ahit.Y());
      funcs.push_back(f);
    }
    int n_funcs = funcs.size();

    //check that the ref value isn't the -1 error value...
    if(ref<0){
      std::cout << "The reference is wrong " << std::endl;
      //also, this function doesn't actually use reference. rethink.
      //It uses the maximum histogram bin position instead
    }


    //find the location of the most intersections, using hist bins.
    TH2F* hist_intersect = new TH2F{"hist_intersect","hist_intersect",90,0, 180, 75,-150,150 };
    for(auto t : Intersections ){
      hist_intersect->Fill(std::get<2>(t),std::get<3>(t));
    }
    double max = hist_intersect->GetMaximumBin();
    double value = hist_intersect->GetBinContent(max);
    Int_t binx,biny,binz;
    hist_intersect->GetBinXYZ(max, binx, biny, binz); 
    double xmaxbin = hist_intersect->GetXaxis()->GetBinCenter(binx);
    double ymaxbin = hist_intersect->GetYaxis()->GetBinCenter(biny);
    std::cout << "binx " << xmaxbin << "biny " << ymaxbin << std::endl;
    delete hist_intersect;


    // find all points within a fairly relaxed threshold around max bin
    // Accumulate list of all the hits that make up this group
    //fThresh = 15;   // this works for regular, not conformal hits
    fThresh = 1;
    std::vector< int > trackhitlist; 
    for(int i_func = 0; i_func < n_funcs; i_func++) {
      auto f1 = funcs[i_func];
      double val = f1->Eval(xmaxbin);
      double dif = ymaxbin-val;
      std::cout << "index: " << i_func << " dif is " << dif << std::endl;
      if (abs(dif) < fThresh) {
        // include this hit in the pretrack for most likely track 
        trackhitlist.push_back(i_func);
        std::cout << "filling trackhitlist with " << i_func << std::endl;
      }
    }
    return trackhitlist;
  } 


  //__________________________________________________________________________________
  //
     //Taking an initial set of function intersections, and a candidate track list,
     //This function checks each intersection and culls it if either hit is on the track list.
     
  std::vector<std::tuple<int,int,double,double>> HoughTransform::ReduceIntersections( std::vector<std::tuple<int,int,double,double>> Intersections,  std::vector< int >  trackhitlist  ) { 

    std::vector<std::tuple<int,int,double,double>> ReducedInt; 
    for(auto t:Intersections) {
      int index1 = std::get<0>(t);
      int index2 = std::get<1>(t);
      //std::cout << "index " << index1 << " and " << index2 << std::endl;
      int addflag = 1;
      for( int i = 0; i<trackhitlist.size(); i++ ) {
         if ( (index1 == trackhitlist[i]) || (index2 == trackhitlist[i]) ){
           addflag = 0;
           //std::cout << trackhitlist[i] << std::endl;
         }
      }
      if ( addflag ==1 ) {
        fIntLeft->Fill(std::get<2>(t),std::get<3>(t));
        ReducedInt.push_back(t);
        std::cout << "filling reduced tuple " << std::get<0>(t) << " " << std::get<1>(t) << std::endl;  
      }
    } 
    return ReducedInt;
  }  
  
  
  //__________________________________________________________________________
  //
     //should return vector of vector of tuple of int and vector of xyztvectors
     //should continue calling until all tracks are found.
     
  std::vector<std::vector<std::tuple<int,double,double>>> HoughTransform::FindTracks( const std::vector<ROOT::Math::XYZTVector>& hits ){
     int num_tracks = 0;
     int max_tracks = 3;
     std::vector<std::vector<std::tuple<int,double,double>>> FoundTracks;

     std::vector<std::tuple<int,int,double,double>> fi = this->FindIntersections(hits);
 
     while (fi.size()>0&&num_tracks<max_tracks){
        int gr = this->GetReferenceHit(hits,fi);
        std::cout<< "intersections have tracks # " << fi.size() << std::endl;
        std::vector< int > trackhitlist = this->FindTrack( hits, fi, gr );

        std::vector<std::tuple<int,double,double>> TL;
        for(int i=0; i<trackhitlist.size(); i++){
           int tracknumber = trackhitlist[i];
           std::tuple<int,double,double> HTL = std::make_tuple(tracknumber,hits[tracknumber].X(),hits[tracknumber].Y());
               //push back std::tuple<int, HIT> ? or some other more clever solution here
           TL.push_back(HTL);
        }
        FoundTracks.push_back(TL); 
        num_tracks++;
        fi = this->ReduceIntersections( fi, trackhitlist );
        std::cout<< " intersections size " << fi.size() << std::endl;
     }
     //return num_tracks;
     return FoundTracks;
  }
  //__________________________________________________________________________
  //
  
  int HoughTransform::FindPeaks(const std::vector<ROOT::Math::XYZTVector>& hits)
  {
    double degree = TMath::Pi()/180.0;
    for(const auto& ahit : hits){
      double phi  = TMath::ATan(ahit.Y()/ahit.X())/degree;
      fPhi->Fill(phi);
    }
    int nfound = fSpec1->Search(fPhi,1,"col",0.2);
    std::cout << "number found peaks is " << nfound << std::endl;
    for(int n = 0; n < nfound; n++) {
      double xphi = (fSpec1->GetPositionX())[n];
      std::cout << "xphi is " << xphi << std::endl;
      fPhiPeaks.push_back(xphi);
    }

    return nfound;
  }
  //__________________________________________________________________________
  
  std::vector<std::vector<std::tuple<ROOT::Math::XYZTVector,ROOT::Math::XYZTVector,ROOT::Math::XYZTVector>>> HoughTransform::GetPreTrackHits(
      const std::vector<std::tuple<ROOT::Math::XYZTVector,ROOT::Math::XYZTVector,ROOT::Math::XYZTVector>>& hits,
      double dphi )
  {

    // note, this function is NOT apparently very consistent
    // even running the same input twice does not guarantee the same result.
    // Peak 1 has 16hits and Peak 1 has 1hits can come from same set of points.
    //
    // copy the hits 
    auto the_hits = hits;
    std::cout <<  " " << the_hits.size() << " hits start\n" ;
    std::vector<std::vector<std::tuple<ROOT::Math::XYZTVector,ROOT::Math::XYZTVector,ROOT::Math::XYZTVector>>> res;
    double degree = TMath::Pi()/180.0;
    for(auto& peak : fPhiPeaks){
      std::vector<std::tuple<ROOT::Math::XYZTVector,ROOT::Math::XYZTVector,ROOT::Math::XYZTVector>> peak_hits;
      the_hits.erase(
          std::remove_if(the_hits.begin(), the_hits.end(), [&](const auto& h){
            double phi  = TMath::ATan(std::get<1>(h).Y()/std::get<1>(h).X())/degree;
            if( TMath::Abs( peak-phi ) < dphi ) {
            peak_hits.push_back( h );
            //std::cout << " peak " << peak << ", " << "phi "  << phi << std::endl;
            //std::cout << " dphi " << dphi << ", " << "TMath::Abs( peak-phi ) "  << TMath::Abs( peak-phi ) << std::endl;
            return true;
            } else {
            return false;
            }
            }), the_hits.end());
      res.push_back(peak_hits);
      std::cout <<  " Peak " << res.size()  << " has " << peak_hits.size() << "hits\n";
    }
    //res.push_back(the_hits);
    //std::cout <<  " " << the_hits.size() << " hits remain\n" ;
    return  res;
  }
  //__________________________________________________________________________

}

