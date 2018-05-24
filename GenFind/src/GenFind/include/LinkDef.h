#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;

#pragma link C++ namespace genfind;

#pragma link C++ class genfind::Hit+;
#pragma link C++ class genfind::ConformalHit+;

#pragma link C++ class genfind::HoughTransform+;

#pragma link C++ class genfind::PreTrackFinder+;
#pragma link C++ class genfind::Conformal2DFinder+;

#pragma link C++ function genfind::compute_conformal_hit(const genfind::Hit&, const genfind::Hit&);


#pragma link C++ class Derp+;
//#pragma link C++ class Event+;
//#pragma link C++ class Track+;

#endif
