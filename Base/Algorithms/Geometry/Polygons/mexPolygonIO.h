/***************************************************************
 * Polygon Clipper / Mex Interface
 *
 * C. Kirst, The Rockefeller University 2014
 *
 ***************************************************************/

#ifndef MEX_POLYGON_IO_H
#define MEX_POLYGON_IO_H

#include <iostream>
using namespace std;

#include "mex.h"
#include "clipper.hpp"

using namespace ClipperLib;

#define IJ(i,j,m) ((j)*m + (i))

void toPolygon(const mxArray *mx, Paths &poly)
{
   int nc;
   int nx;
   long64 *xydat;
   const mxArray *xy;
   
   nc = mxGetNumberOfElements(mx);
   //cout << "nc: " << nc << endl; cout.flush();
   
   poly.clear();
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

int fromPolygonTreeRecursive(PolyNode* pn, mxArray*& pol, double* treedat, int n, int id) 
{
   int parent = id;   
   //cout << "from id: " << id << " parent:" << parent << endl; cout.flush();

   for (unsigned i = 0; i < pn->ChildCount(); ++i)
   {
      id++;
      
      PolyNode* pc = pn->Childs[i];

      int nxy = pc->Contour.size();
      //cout << "from id:"<< id << " nxy: " << nxy << endl; cout.flush();

      mxArray* xy = mxCreateDoubleMatrix(2, nxy, mxREAL);
      double* xydat = mxGetPr(xy);

      unsigned k = 0;
      for (unsigned j = 0; j < nxy; j++) {
         xydat[k]   = pc->Contour[j].X;
         xydat[k+1] = pc->Contour[j].Y;
         k+=2;
      }

      mxSetCell(pol, id, xy);

      //construct the tree stucture
      if (parent >= 0) {
         treedat[IJ(id, parent, n)] = 1;
      }
            
      id = fromPolygonTreeRecursive(pc, pol, treedat, n, id);
   }
   
   return id;
}

void fromPolygonTree(PolyTree* polytree, mxArray*& pol, mxArray*& tree)
{   
   int nc = polytree->AllNodes.size();
   //cout << "from nc: " << nc << endl; cout.flush();

   int ndim=2, dims[]={1, nc};
   pol = mxCreateCellArray(ndim, dims);
   tree = mxCreateDoubleMatrix(nc,nc,mxREAL);
   double* treedat = mxGetPr(tree);

   fromPolygonTreeRecursive((PolyNode*) polytree, pol, treedat, nc, -1);
}

#endif
   