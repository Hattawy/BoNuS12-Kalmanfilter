#ifndef genfind_GenFindHits_HH
#define genfind_GenFindHits_HH

#include <memory>
#include <functional>
#include <algorithm>
#include "TObject.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

namespace genfind {


  //__________________________________________________________________________

  /** Interface hit class.
   */
  class Hit {

    public:
      ROOT::Math::XYZTVector fPosition = {0,0,0,0};

      //DD4Hep::Simulation::Geant4Tracker::Hit* fG4TrackerHit;

    public:
      Hit();
      Hit(const ROOT::Math::XYZTVector& );
      Hit(double x, double y, double z, double t=0);
      Hit(const Hit&) = default;
      Hit(Hit&&) = default;
      Hit& operator=(const Hit&) = default;
      Hit& operator=(Hit&&) = default;
      virtual ~Hit();
      //virtual ~Hit() = default;
 
      double X() const { return fPosition.X(); }
      double Y() const { return fPosition.Y(); }
      double Z() const { return fPosition.Z(); }
      double T() const { return fPosition.T(); }

      ClassDef(Hit,1)
  };


  /** Conform space hits.
   */
  class ConformalHit {
    public:
      ROOT::Math::XYZTVector fPosition  = {0,0,0,0};
      std::shared_ptr<genfind::Hit>          fImageHit;
      std::shared_ptr<genfind::Hit>          fImageRef;

      double X() const { return fPosition.X(); }
      double Y() const { return fPosition.Y(); }
      double Z() const { return fPosition.Z(); }
      double T() const { return fPosition.T(); }

    public:
      ConformalHit();
      ConformalHit(const ConformalHit&) = default;
      ConformalHit(ConformalHit&&) = default;
      ConformalHit& operator=(const ConformalHit&) = default;
      ConformalHit& operator=(ConformalHit&&) = default;
      virtual ~ConformalHit();

      ClassDef(ConformalHit,1)
  };

  /** @brief Compute the conformal transform of a hit.
   *
   * @param hit The hit that provides the position space coordinates.
   * @param ref Hit providing the reference coordinates (x0,y0).
   *
   * \f$ u = \frac{x}{x^2+y^2} \f$
   * \f$ v = \frac{y}{x^2+y^2} \f$
   *
   */
  ConformalHit compute_conformal_hit(std::shared_ptr<Hit> hit, std::shared_ptr<Hit> ref); 

  ConformalHit compute_conformal_hit(std::shared_ptr<Hit> hit); 

  template<class T>
  ConformalHit compute_conformal_hit(const T& hit, const T& ref) {
    ConformalHit chit;
    chit.fImageHit = std::make_shared<genfind::Hit>(hit.X(), hit.Y(), 0.0) ;
    chit.fImageRef = std::make_shared<genfind::Hit>(ref.X(), ref.Y(), 0.0) ;
    double x = hit.X() - ref.X();
    double y = hit.Y() - ref.Y();
    double R = x*x + y*y;
    if( !(hit == ref) ) {
      chit.fPosition = { 2.0*x/R, 2.0*y/R, hit.Z(), hit.T() } ;
    }
    return chit;
  }

  template<class T>
  std::vector<ConformalHit> compute_conformal_hits(const std::vector<T>& hits, const T& ref)
  {
    // Local class  
    struct Sum {
      std::vector<ConformalHit> fChits;
      T                         fRef;
      Sum(const T& R) : fChits{}, fRef(R) { }
      void operator()(const T& h){
        fChits.push_back(compute_conformal_hit(h, fRef));
      }
      std::vector<ConformalHit> Get(){ return fChits; }
    };
    return std::for_each(hits.begin(), 
                         hits.end(), 
                         Sum(ref)
                        ).Get();
  } 


}


#endif
