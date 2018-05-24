#include "GenFindHits.h"

namespace genfind {

  Hit::Hit()
  { }
  //__________________________________________________________________________

  Hit::Hit(const ROOT::Math::XYZTVector& p) : fPosition(p)
  { }
  //__________________________________________________________________________
  
  Hit::Hit(double x, double y, double z, double t) : fPosition({x,y,z,t})
  { }
  //__________________________________________________________________________

  Hit::~Hit()
  { }
  //__________________________________________________________________________

  ConformalHit::ConformalHit()
  { }
  //__________________________________________________________________________

  ConformalHit::~ConformalHit()
  { }
  //__________________________________________________________________________

  ConformalHit compute_conformal_hit(std::shared_ptr<Hit> hit, std::shared_ptr<Hit> ref)
  {
    ConformalHit chit;
    chit.fImageHit = hit;
    chit.fImageRef = ref;
    double x = hit->X() - ref->X();
    double y = hit->Y() - ref->Y();
    double R = x*x + y*y;
    if( !(hit == ref) ) {
      chit.fPosition = { 2.0*x/R, 2.0*y/R, hit->Z(), hit->T() } ;
    }
    return chit;
  }
  //__________________________________________________________________________

  ConformalHit compute_conformal_hit(std::shared_ptr<Hit> hit)
  {
    // Use ref = {0,0};
    return compute_conformal_hit(hit, std::make_shared<Hit>());
  }
  //__________________________________________________________________________

}

