/***************************************************************
 * Polygon Dilation / Buffering
 *
 * Based on clipper.cpp library
 *
 * C. Kirst, The Rockefeller University 2014
 *
 ***************************************************************/

// mex mexPolygonBuffer.cpp clipper.cpp 

#include <iostream>
using namespace std;

#include "mex.h"
#include "clipper.hpp"
#include "mexPolygonIO.h"

using namespace ClipperLib;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  
   if (nrhs < 2) mexWarnMsgTxt("mexPolygonBuffer: 2 input arguments expected!");
   if (nlhs < 1) mexWarnMsgTxt("mexPolygonBuffer: 1 output argument expected!");
   
   //Polygon 
   Paths poly;
   toPolygon(prhs[0], poly);
   
   double r = mxGetScalar(prhs[1]);

   //JoinType {jtSquare, jtRound, jtMiter}
   JoinType jT = jtMiter;
   if (nrhs > 2) jT = (JoinType) (int) mxGetScalar(prhs[2]);

   //EndType {etClosedPolygon, etClosedLine, etOpenButt, etOpenSquare, etOpenRound};
   EndType eT = etClosedPolygon;
   if (nrhs > 3) eT = (EndType) (int) mxGetScalar(prhs[3]); 


   ClipperOffset co;
   co.AddPaths(poly, jT, eT);
   
   if (nlhs < 2) {
      Paths solution;
      co.Execute(solution, r);
      plhs[0] = fromPolygon(solution);
      
   } else { //poly + tree
      PolyTree polytree;
      co.Execute(polytree, r);
      fromPolygonTree(&polytree, plhs[0], plhs[1]);
   }
   
   return;
};



