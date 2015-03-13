/***************************************************************
 * Polygon Union/Intersection/Difference/XOR
 *
 * Based on clipper.cpp library
 *
 * C. Kirst, The Rockefeller University 2014
 *
 ***************************************************************/

// clc; mex mexPolygonExecute.cpp clipper.cpp 

#include "mex.h"
#include "clipper.hpp"

#include <iostream>
using namespace std;

using namespace ClipperLib;

#include "mexPolygonIO.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   //Input: operation, fillType, subj, (subjclosed), clip1, (clip1closed), clip2, ...
   if (nrhs <  3) mexWarnMsgTxt("mexPolygonExecute: 3 input arguments expected: operation, fillType, subj, (subjclosed), (clip1), (clip1closed), (clip2), ...!");
   
   //Output: solultion (tree)
   if (nlhs != 1 && nlhs !=2) mexWarnMsgTxt("mexPolygonExecute: 1-2 output arguments expected!");
   
   //Clipper 
   Clipper c;
   
   //Operation
   ClipType clipT = (ClipType) (int) mxGetScalar(prhs[0]);

   //PolyFillType 
   //see http://glprogramming.com/red/chapter11.html
   //PolyFillType {pftEvenOdd, pftNonZero, pftPositive, pftNegative };
   
   PolyFillType subjFillType = pftEvenOdd; 
   PolyFillType clipFillType = pftEvenOdd;
   
   int nf = mxGetNumberOfElements(prhs[1]);
   double* type = mxGetPr(prhs[1]);

   if (nf > 0) {
      subjFillType = (PolyFillType) int(*type);
   }
   if (nf > 1) {
      type++;
      clipFillType = (PolyFillType) int(*type);
   }
   //cout << "fill types: " << (int) subjFillType << " " << (int) clipFillType << endl; cout.flush();

   //Polygon 
   Paths poly;
   int n = 2;
   bool closed = true;   

   //Subject Poly
   toPolygon(prhs[n], poly);
   n++;
   if (n < nrhs && !mxIsCell(prhs[n])) {
      closed = mxGetScalar(prhs[n]) > 0;
      n++;
   } 
   c.AddPaths(poly, ptSubject, closed);

   //Clipping Polys
   for(; n < nrhs; n++) {
      //cout << "clipping poly: " << n << endl; cout.flush();
      if (mxIsCell(prhs[n])) {
         toPolygon(prhs[n], poly);

         if (n + 1 < nrhs && !mxIsCell(prhs[n+1])) {
            n++;
            closed = mxGetScalar(prhs[n]) > 0;
         } else {
            closed = true;
         }
         
         c.AddPaths(poly, ptClip, closed);
         
      } else {
         break;
      }
   }
  
   if (n < nrhs) {
      mexWarnMsgTxt("mexPolygonExecute: inconsistent input arguments in polygons!");
   }

   
   //Execute
   //cout << "execute: " << n << endl; cout.flush();   
   if (nlhs < 2) {
      Paths solution;
      c.Execute(clipT, solution, subjFillType, clipFillType);
      plhs[0] = fromPolygon(solution);
      
   } else { //poly + tree
      
      PolyTree polytree;
      c.Execute(clipT, polytree, subjFillType, clipFillType);      
      fromPolygonTree(&polytree, plhs[0], plhs[1]);
   }
   
   return;
};







