/*  Vignetting Correction Package
 *
 *  LUT - Lookup Table
 *
 *
 *  Based on Hugin photometric corrections
 *
 *  C Kirst 2015 Rockefeller University
 *
 */

#ifndef VIGNETTING_LUT_H
#define VIGNETTING_LUT_H

#include <assert.h>
#include <algorithm>

#include "EMoR.h"


namespace Photometry {
    
// Lookup table for the Camera Response
    
class LUT {
   
public:
   //look up table
   std::vector<double> lut;
   
public:
   LUT() {};
   ~LUT() {};
   
public:
      
   
   inline void fromEMoRLUT(const std::vector<double> & params)
   {
      const double s = (double) 1.0; //LUTTraits<VT>::max();
      
      // lookup tables
      size_t nDim = params.size();
      assert(nDim < 26);
      
      lut.resize(1024);
      for (int i=0; i<1024; ++i) {
         double t = Photometry::EMoR::f0[i];
         for (size_t j=0; j < nDim; j++) {
            t += params[j] * Photometry::EMoR::h[j][i];
         }
         lut[i] = t*s;
      }
   }
    
   void resize(int outsize) {
      
      assert(lut.size());
      assert(outsize);
      
            
      std::vector<double> lut_new;
      
      for(size_t oIdx = 0; oIdx < outsize; oIdx++) {
         double ix = oIdx/(outsize-1.0) * (lut.size()-1);
         unsigned iIdx = unsigned(ix);
         double deltaix = ix-iIdx;
         if (deltaix == 0.0) {
            // no interpolation
            lut_new[oIdx] = lut[iIdx];
         } else if (iIdx+1 <= lut.size()){
            // linear interpolation
            lut_new[oIdx] = (1-deltaix) * lut[iIdx] + deltaix * lut[iIdx+1];
         } else {
            lut_new[oIdx] = lut.back();
         }
      }
      lut = lut_new;
   }
   

   void enforceMonotonicity()
   {
      int lutsize = lut.size();

      if (lutsize) {
         double max = lut.back();
         for (int j=0; j < lutsize-1; j++)
         {
            if (lut[j+1] > max) {
               lut[j+1] = max;
            } else if (lut[j+1] < lut[j]) {
               lut[j+1] = lut[j];
            }
         }
		}
	}
   
   double apply(double v) const
   {
      assert(lut.size() > 0);
      
      if (v > 1) return lut.back();
      if (v < 0) return 0;
      
      double x = v * (lut.size()-1);
      unsigned i = unsigned(x);
      
      // interpolate
      x = x-i;
      if ( i+1 < lut.size()) {
         return (1-x) * lut[i] + x * lut[i+1];
      } else {
         return lut[i];
      }
   }
   
   double apply(std::vector<double> v) const {
      for (int i = 0; i < v.size(); i++) {
         v[i] = apply(v[i]);
      }
   }
   
   double apply(double* v, int n) const {
      for (int i = 0; i < n; i++) {
         v[i] = apply(v[i]);
      }
   }
   
   
   
   double applyInverse(double v) const {
      // assume float is scaled 0..1
      assert(lut.size() > 0);
      
      if (v >= lut.back()) return lut.back();
      if (v < lut[0])      return 0;
      
      // find the lower bound, p will point to the first *p >= v
      std::vector<double>::const_iterator p = std::lower_bound(lut.begin(), lut.end(), v);
      
      int x = p-lut.begin();
      
      if (v == 1) {
         return 1;
      } else if (x == 0) {
         return 0;
      } else if (v == *p) {
         return x/(lut.size()-1.0);
      } else {
         // interpolate position.
         // p points to the first element > v
         double lower = *(p-1);
         double upper = *(p);
         double delta =  (v - lower) / (upper - lower);
         return (x-1 + delta) / (lut.size()-1.0);
      }
   }
   
   double applyInverse(std::vector<double> v) const {
      for (int i = 0; i < v.size(); i++) {
         v[i] = applyInverse(v[i]);
      }
   }
   
   double applyInverse(double* v, int n) const {
      for (int i = 0; i < n; i++) {
         v[i] = applyInverse(v[i]);
      }
   }
   
}; // class LUT

  
} // namespace Photometry

#endif
