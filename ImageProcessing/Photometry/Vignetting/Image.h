/*  Image.h
 *
 *  representation of images
 *  interface with mex
 *
 *  C Kirst 2015 Rockefeller University
 *
 *  TODO: extend this to different types and color images
 */

#include "mex.h"

class image {

public:
   double * src;
   int n,m;
   
public:
   
   void image2D(double* s, int nn, int mm): src(s), n(nn), m(mm) {}:
      
   void ~image2D() {};

   TYPE operator() (const int i, const  int j) {
      return src[((j)*m + (i))];
   }

   void fromMxArray(mxArray* array) {
      const mwSize ndim = mxGetNumberOfDimensions(array);;
      if (2 != ndim){
         mexErrMsgTxt("input image must be 2d grayscale");
      }

      const mwSize* size = mxGetDimensions(array);

      const mxClassID  category = mxGetClassID(array);
      if (mxClassID != mxDOUBLE_CLASS){
         mexErrMsgTxt("input image class must be double");
      }
 
      m = size[0];
      n = size[1];
   }
   
   int width() {
      return m;
   }
   
   int height() {
      return n;
   } 
}