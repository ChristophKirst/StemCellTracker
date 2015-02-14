/***************************************************************
 * Polygon Dilation / Buffering
 *
 * Based on clipper.cpp library
 *
 * C. Kirst, The Rockefeller University 2014
 *
 ***************************************************************/


// mex mexPolygonBuffer.cpp clipper.cpp 

#include "mex.h"
#include "clipper.hpp"

#include <iostream>
using namespace std;

using namespace ClipperLib;

void toPolygon(const mxArray *mx, Paths &poly)
{
   int nc;
   int nx;
   long64 *xydat;
   const mxArray *xy;
   
   nc = mxGetNumberOfElements(mx);
   //cout << "nc: " << nc << endl; cout.flush();
   
   poly.resize(nc);
   for (unsigned i = 0; i < nc; i++){
      
      xy = mxGetCell(mx, i);
      unsigned nxy = mxGetNumberOfElements(xy);
      //cout << "nxy: " << nxy << " " << nxy / 2 <<  endl; cout.flush();
      
      poly[i].resize(nxy/2);
      
      xydat = (long64*)mxGetData(xy);
      for (unsigned j = 0; j < nxy/2; j++){
         poly[i][j].X = xydat[2*j];
         poly[i][j].Y = xydat[2*j+1];
      }
   }
}

mxArray* fromPolygon(Paths& poly)
{   
   int nc = poly.size();
   //cout << "nc: " << nc << endl; cout.flush();
   
   int ndim=2, dims[]={1, nc};
   mxArray* mx = mxCreateCellArray(ndim, dims);
 
   for (unsigned i = 0; i < nc; ++i)
   {
      int nxy = poly[i].size();
      //cout << "nxy: " << nxy << endl; cout.flush();

      mxArray* xy = mxCreateDoubleMatrix(2, nxy, mxREAL);
      double* xydat = mxGetPr(xy);

      unsigned k = 0;
      for (unsigned j = 0; j < nxy; j++)
      {
         xydat[k]   = poly[i][j].X;
         xydat[k+1] = poly[i][j].Y;
         k+=2;
      }
      mxSetCell(mx, i, xy);
   }
   
   return mx;
}
   
   
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  
   if (nrhs <  2) mexWarnMsgTxt("mexPolygonBuffer: 2 input arguments expected!");
   if (nlhs != 1) mexWarnMsgTxt("mexPolygonBuffer: 1 output argument expected!");
   
   //Polygon 
   Paths poly;
   toPolygon(prhs[0], poly);
   
   double r = mxGetScalar(prhs[1]);

   //JoinType {jtSquare, jtRound, jtMiter}
   JoinType jT = jtMiter;
   if (nrhs > 2) jT = (JoinType) mxGetScalar(prhs[2]);

   //EndType {etClosedPolygon, etClosedLine, etOpenButt, etOpenSquare, etOpenRound};
   EndType eT = etClosedPolygon;
   if (nrhs > 3) eT = (EndType) mxGetScalar(prhs[3]); 

   Paths solution;
   
   ClipperOffset co;
   co.AddPaths(poly, jT, eT);
   co.Execute(solution, r);

   plhs[0] = fromPolygon(solution);
   
   return;
};



