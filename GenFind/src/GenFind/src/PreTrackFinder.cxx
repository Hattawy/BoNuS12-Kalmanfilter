#include "PreTrackFinder.h"
#include <algorithm>

namespace genfind {

  PreTrackFinder::PreTrackFinder()
  { }
  //__________________________________________________________________________

  PreTrackFinder::~PreTrackFinder()
  { } 
  //__________________________________________________________________________


  Conformal2DFinder::Conformal2DFinder()
  {
    bool status = TH1::AddDirectoryStatus();
    TH1::AddDirectory(false);
    fPhi       = new TH1F{"hPhi","hPhi",90,-90, 90};
    fSpec1 = new TSpectrum(2*fMaxTracks);
    TH1::AddDirectory(status);
  } 
  //__________________________________________________________________________

  //Conformal2DFinder::Conformal2DFinder(const std::vector<genfind::ConformalHit*> h) :
  //   Conformal2DFinder()
  //{
  //  fInputHits = h;
  //}
  //__________________________________________________________________________

  Conformal2DFinder::~Conformal2DFinder()
  {
    delete fPhi;
    if(fSpec1) delete fSpec1; fSpec1 = nullptr;
  }
  //__________________________________________________________________________
  std::vector<std::shared_ptr<genfind::ConformalHit>> 
    Conformal2DFinder::GenerateConformalHits(
        const std::vector<genfind::Hit*>& hits
        )
    {
      fInputHits.clear();
      fConformalHits.clear();
      for(const auto& h : hits) {
        fInputHits.push_back( std::make_shared<genfind::Hit>(*h) );
      }

      auto ref = fInputHits[0];
      std::for_each(fInputHits.begin(), fInputHits.end(), [&,this](auto p) {
          this->fConformalHits.push_back( std::make_shared<genfind::ConformalHit>(
                compute_conformal_hit(p, ref)) );
          });
      return fConformalHits;
    }
  //__________________________________________________________________________

  std::vector<std::vector<std::shared_ptr<Hit>>> Conformal2DFinder::GroupHits() const
  {
    return {fInputHits};
  }
  //__________________________________________________________________________

  std::vector<std::vector<std::shared_ptr<Hit>>> Conformal2DFinder::GroupHits(const std::vector<genfind::Hit*>& hits)
  {
    fInputHits.clear();
    for(const auto& h : hits) {
      fInputHits.push_back( std::make_shared<genfind::Hit>(*h) );
    }
    return GroupHits();
  }

}

