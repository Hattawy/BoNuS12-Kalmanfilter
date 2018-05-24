#ifndef genfind_HoughTransform_HH
#define genfind_HoughTransform_HH

#include "TObject.h"
#include "TMath.h"
#include <array>
#include <functional>
#include <vector>
#include <utility>
#include "TH2F.h"
#include "TH1F.h"
#include "Math/Vector4D.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TMultiGraph.h"
#include "GenFindHits.h"

namespace genfind {


  class HoughTransform : public TObject {
    public:

      TH1F*        fPhi       = nullptr;
      TH2F*        fRhoTheta0 = nullptr;
      TH2F*        fRhoTheta1 = nullptr;
      TH2F*        fRhoTheta2 = nullptr;
      TH2F*        fIntersect = nullptr;
      TH2F*        fIntReduced= nullptr;
      TH2F*        fIntLeft   = nullptr;


      TSpectrum*   fSpec1     = nullptr;
      TSpectrum2*  fSpec      = nullptr;

      int fThresh    = 3;
      int fMaxTracks = 20;

      std::vector<TF1*> fFuncs;
      std::vector<std::tuple<int,int,double>> fHTRoots;
      std::vector<std::tuple<int,int,double,double>> fHTIntersections;
      std::vector<double>    fPhiPeaks;

    public:
      HoughTransform();
      virtual ~HoughTransform();

      std::vector<ROOT::Math::XYZTVector> GetConformalCoordinates(
          const ROOT::Math::XYZTVector ref_point,
          const std::vector<ROOT::Math::XYZTVector>& hits) const;

      std::vector<ROOT::Math::XYZTVector> UseConformalMap( 
          const::std::vector<ROOT::Math::XYZTVector>& hits) const;

//====== delete the following from a merge if it conflicts with some of the other stuff.
//
//      /** @brief Method to get the result. 
//       *
//       * @param hits Vector of hit type T which.
//       *
//       * Returns vector of vector of hits group into possilbe track candidates. 
//       * Note: the reference hit for the conformal transform is at the origin.
//       *
//       */
//      template<class T>
//      std::vector<std::vector<T>> operator()(const std::vector<T>& hits){
//        std::vector<std::vector<T>> res; 
//        T ref = {0.0,0.0,0.0,0.0};
//        auto chits = genfind::compute_conformal_hits(hits, ref);
//        auto peaks = FindPhiPeaks(chits);
//        res = SelectPhiPeaks<T>(peaks, chits);
//        return res;
//      }
//
//      template<class T>
//      std::vector<std::vector<T>> SelectPhiPeaks(const std::vector<double>& peaks,
//                                                 const std::vector<genfind::ConformalHit>& chits)
//      {
//        double degree = TMath::Pi()/180.0;
//        std::vector<std::vector<T>> res;
//        for(auto& p : peaks){
//          std::vector<T> p_hits;
//          for(const auto& ahit : chits){
//            double phi  = TMath::ATan(ahit.Y()/ahit.X())/degree;
//            if( TMath::Abs(phi-p) < fPhiThresh ) {
//              p_hits.push_back(
//                {ahit.fImageHit->X(),
//                  ahit.fImageHit->Y(), 
//                  ahit.fImageHit->Z(),
//                  ahit.fImageHit->T()});
//            }
//          }
//          res.push_back(p_hits); 
//        }
//        return res;
//      }
//
//      /** @brief Fill the Hough transform histograms for searching.
//       *
//       * @param chits Vector of ConformalHits used to find lines.
//       *
//       */
//      void FillHistograms(const std::vector<genfind::ConformalHit>& chits);
//      
//      std::vector<double> FindPhiPeaks(const std::vector<genfind::ConformalHit>& chits);

// next line, end deletes if the merge caused a conflict.
//>>>>>>> d267ad98e7036905dfceb221e8b3f316819b0990

      TMultiGraph* FillHoughTransforms(const std::vector<ROOT::Math::XYZTVector>& hits);

      std::vector<TF1*> FillFunctions(const std::vector<ROOT::Math::XYZTVector>& hits);

      std::vector<std::tuple<int,int,double,double>> FindIntersections( const std::vector<ROOT::Math::XYZTVector>& hits);

      int GetReferenceHit( const std::vector<ROOT::Math::XYZTVector>& hits, std::vector<std::tuple<int,int,double,double>> Intersections );

      std::vector< int > FindTrack( const std::vector<ROOT::Math::XYZTVector>& hits, std::vector<std::tuple<int,int,double,double>> Intersections, int ref );

      std::vector<std::tuple<int,int,double,double>> ReduceIntersections(std::vector<std::tuple<int,int,double,double>> Intersections, std::vector< int > trackhitlist );

      int FindPeaks(const std::vector<ROOT::Math::XYZTVector>& hits);

      std::vector<std::vector<std::tuple<int,double,double>>> FindTracks( const std::vector<ROOT::Math::XYZTVector>& hits );

      //return type for find tracks
      //        innermost, hit element of a tracklist, the hit index, plus the hit: 
      //                std::tuple<int,std::vector<XYZTVector>> typename HTL
      //        Next, one full tracklist. Vector of HTLs
      //                std::vector<std::tuple<int,std::vector<XYZTVector>>> typename TL
      //        next, full findtracks output. Vector of tracklists.
      //                std::vector<std::vector<std::tuple<int,std::vector<XYZTVector>>>> typename FT                
      //
      // crude idea:
      //struct PreTrack {
      //  double phi_peak;
      //  ROOT::Math::XYZTVector ref_hit; // image or conformal space?
      //  std::vector<ROOT::Math::XYZTVector> hits; // image space
      //  std::vector<ROOT::Math::XYZTVector> chits; // conformal space
      //}

      //PreTrack GetPreTrack(
      //    const std::vector<std::tuple<ROOT::Math::XYZTVector,ROOT::Math::XYZTVector>>& hits,
      //    ROOT::Math::XYZTVector ref_hit,
      //    double dphi = 5.0 /*degrees*/);

      std::vector<std::vector<std::tuple<ROOT::Math::XYZTVector,ROOT::Math::XYZTVector,ROOT::Math::XYZTVector>>> GetPreTrackHits(
          const std::vector<std::tuple<ROOT::Math::XYZTVector,ROOT::Math::XYZTVector,ROOT::Math::XYZTVector>>& hits,
          double dphi = 5.0 /*degrees*/);

      ClassDef(HoughTransform,1);
  };

}

#endif


