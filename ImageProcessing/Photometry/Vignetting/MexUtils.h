/* Mex Utils
 *
 * C Kirst 2015 Rockefeller University
 *
 */

//TODO: Traits Class conversions / Template versions

#ifndef MEX_UTILS_H
#define MEX_UTILS_H

#include "mex.h"

mxArray * mexUtilsVectorToArray(const std::vector<double>& v){
   mxArray * mx = mxCreateDoubleMatrix(1, v.size(), mxREAL);
   std::copy(v.begin(), v.end(), mxGetPr(mx));
   return mx;
}

std::vector<double> mexUtilsArrayToVector(const mxArray* mx) {
   mwSize size = mxGetNumberOfElements(mx);
   double * data = mxGetPr(mx);
   return std::vector<double>(data, data + (int) size);
}

double mexUtilsArrayElement(const mxArray* mx, int n) {
   double * data = mxGetPr(mx);
   return data[n];
}

double mexUtilsArrayElement(const mxArray* mx, int i, int j) {
   int m = mxGetM(mx);
   double * data = mxGetPr(mx);
   return data[j*m + i];
}

#endif