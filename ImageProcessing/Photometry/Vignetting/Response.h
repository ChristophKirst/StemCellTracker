/*  Vignetting Correction Package
 *
 *  Response that maps between corrected and 
 *
 * 
 *  Based on Hugin photometric corrections
 *
 *  C Kirst 2015 Rockefeller University
 *
 */

#ifndef VIGNETTING_RESPONSE_H
#define VIGNETTING_RESPONSE_H

#include <math.h>
#include <vector>
#include "LUT.h"

namespace Photometry {
    
/*  Transform from light irradinance input to camrea response output 
 *  radiometric transformation, includes exposure, vignetting and white balance.
 */
    
class Response
{
   public:
      
      //camera response
      enum ResponseType {LINEAR_RESPONSE = 0, EMOR_RESPONSE};
      ResponseType type;
      std::vector<double> responseParameter;
      LUT lut;
      
      //exposure
      double exposure;
      
      //vignetting
      double radiusScale;
      std::vector<double> vignettingParameter;
      double vignettingCenterX;
      double vignettingCenterY;
      
      //white balance
      double whiteBalanceRed;
      double whiteBalanceBlue;
      
      //flatfield
      //image* flatfield;

   public:

      Response(double width, double height, const ResponseType& tp = EMOR_RESPONSE) {
         init(width, height, tp);
      }
      
      virtual ~Response() {};

   public:
         
      void init(double width, double height, const ResponseType& tp = EMOR_RESPONSE) {
         // exposure
         exposure = 1.0;
         
         //init vignetting
         radiusScale = 1.0/sqrt(width/2.0*width/2.0 + height/2.0*height/2.0);
         vignettingCenterX = width / 2.0;
         vignettingCenterY = height / 2.0;
         
         vignettingParameter[0] = 1.0;
         for (int i = 1; i < 4; i++) {
            vignettingParameter[i] = 0;
         }
         
         //response
         type = tp;
         if (type == EMOR_RESPONSE) {
            initEMoRLUT();
         }  
      }
      
      void initEMoRLUT() {
         responseParameter.resize(6);
         for (int i = 0; i < 6; i++) {
            responseParameter[i] = 0;
         }
         type = EMOR_RESPONSE;
         lutFromResponseParameter();
      }
      
      void initEMoRLUT(const std::vector<double>& respParameter) {
         responseParameter = respParameter;
         type = EMOR_RESPONSE;
         lutFromResponseParameter();
      }

         
      void lutFromResponseParameter() {
         // build response function lookup table
         lut.fromEMoRLUT(responseParameter);
      }

      double getVignettingFactor(double x, double y) const {
         
         x = x - vignettingCenterX;
         y = y - vignettingCenterY;
         
         x *= radiusScale;
         y *= radiusScale;
         
         double vig = vignettingParameter[0];
         double r2 = x * x + y * y;
         double r = 1.0;
         for (unsigned int i = 1; i < 4; i++) {
            r *= r2;
            vig += vignettingParameter[i] * r;
         }
         return vig;
      }


      double apply(double v, double x, double y) const {
         //vignetting
         double vig = getVignettingFactor(x,y);
         
         //exposure
         double ret = v * vig * exposure;
         
         //response
         if (type) {
            return lut.apply(ret);
         } else {
            return ret;
         }
      }
      
      double applyVector(int n, double* v, double* x, double* y, double* r) const {
         for (int i = 0; i < n; i++) {
            r[i] = apply(v[i], x[i], y[i]);
         }
      }

      double applyInverse(double v, double x, double y) const {
         // inverse response
         double ret = v;
         if (type) {
            ret = lut.applyInverse(v);
         }
         
         // inverse vignetting and exposure
         double vig = getVignettingFactor(x,y);
         return ret / (vig * exposure);
      }
      
      double applyInverseVector(int n, double* v, double* x, double* y, double *r) const {
         for (int i = 0; i < n; i++) {
            r[i] = applyInverse(v[i], x[i], y[i]);
         }
      }
      
      
};


} // namespace Photometry

#endif // VIGNETTING_RESPONSE_H
