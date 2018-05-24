//*********************************************************************
// GenFind-- tool to select candidate tracks from a set of hit events. 
// Built to be generic with respect to detector geometry
// Meant to work together with GenFit, and Marlin, to find and reconstruct tracks
// 
// Version 0.0 Feb 21, 2017 by Sereres Johnston
// *****************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "GenFindConfig.h"
#include "utilities/resolution.h"

int main (){

  printf("I'm a stupid test right now-- cmake compiled\n");
  printf("I also have a version number %d.%d\n",GenFind_VERSION_MAJOR,GenFind_VERSION_MINOR);
  double res=1;
  auto res_obj = GenFindUtility::Resolution();
  res = res_obj.GetRes();
  printf("plus I can access a library? GetRes returns %f \n",res);

  return 0;
}



