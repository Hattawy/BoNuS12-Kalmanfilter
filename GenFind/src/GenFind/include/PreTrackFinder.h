#ifndef genfind_PreTrackFinder_HH
#define genfind_PreTrackFinder_HH

#include "TObject.h"
#include <array>
#include <vector>
#include <utility>
#include "GenFindHits.h"
#include "TH2F.h"
#include "TH1F.h"
#include "Math/Vector4D.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TMultiGraph.h"

namespace genfind {

  /** Pre-track Finder.
   *  Abstract base class.
   */
  class PreTrackFinder : public TObject {

    public:
      PreTrackFinder();
      PreTrackFinder(const PreTrackFinder&) = default;
      PreTrackFinder(PreTrackFinder&&) = default;
      PreTrackFinder& operator=(const PreTrackFinder&) = default;
      PreTrackFinder& operator=(PreTrackFinder&&) = default;
      virtual ~PreTrackFinder();

      virtual std::vector<std::vector<std::shared_ptr<Hit>>> GroupHits() const = 0;
      virtual std::vector<std::vector<std::shared_ptr<Hit>>> GroupHits(const std::vector<genfind::Hit*>& hits)  = 0;

      ClassDef(PreTrackFinder,1)
  };


  /**  Simple 2D conformal pre-track finder.
   */
  class Conformal2DFinder : public PreTrackFinder {
    protected:
      std::vector<std::shared_ptr<genfind::ConformalHit>> fConformalHits;
      std::vector<std::shared_ptr<genfind::Hit>>          fInputHits;

      double       fdphi = 4.0;
      TSpectrum*   fSpec1     = nullptr;
      TH1F*        fPhi       = nullptr;

      double fThresh = 0.05;
      int fMaxTracks = 20;

      std::vector<double>    fPhiPeaks;

    public:
      Conformal2DFinder();
      //Conformal2DFinder(const std::vector<genfind::Hit*> h);
      virtual ~Conformal2DFinder();

      //std::vector<std::shared_ptr<genfind::ConformalHit>> GenerateConformalHits();
      std::vector<std::shared_ptr<genfind::ConformalHit>> GenerateConformalHits(
            const std::vector<genfind::Hit*>& hits);

      virtual std::vector<std::vector<std::shared_ptr<Hit>>> GroupHits() const override;
      virtual std::vector<std::vector<std::shared_ptr<Hit>>> GroupHits(const std::vector<genfind::Hit*>& h) override;

      ClassDef(PreTrackFinder,1)
  };
}

#endif


