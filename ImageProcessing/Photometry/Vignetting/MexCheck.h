/* 
 * 
 * Input Argument Checkking for Mex
 *
 *  C Kirst 2015 Rockefeller University
 *
 */

#ifndef MEX_CHECK_H
#define MEX_CHECK_H

#include "mex.h"

#include <sstream>
#include <string>

class MexCheckError {
   public:
      std::string routine;
      std::string variable;
      
   public:
      MexCheckError() {};
      
      MexCheckError(const std::string& rt) : routine(rt) {};
      
      MexCheckError(const std::string& rt, const std::string& vn) : routine(rt), variable(vn) {};

      ~MexCheckError() {};
         
      std::string error() const {
         std::stringstream str;
         if (routine.empty()) {
            str << "Error"; 
         }
         else {
            str << routine << ": error";
         }
         if (!variable.empty()) str << ": inconsistent input for " << variable;
         else str << "!";
         return str.str();
      }
      
      std::string classError() const {
         std::stringstream str;
         str << error() << " wrong class!";
         return str.str();
      }
      
      std::string classError(const mxClassID& cid) const {
         std::stringstream str;
         str << error() << "! Not of class: " << int(cid) << "!";
         return str.str();
      }
      
      std::string classError(const mxClassID& cid, const mxClassID& cidfound) const {
         std::stringstream str;
         str << classError(cid) << " Found: " << int(cid) << "!";
         return str.str();
      }
      
      std::string dimError() const {
         std::stringstream str;
         str << error() << " wrong dimension!";
         return str.str();
      }

      std::string dimError(int dim) const {
         std::stringstream str;
         str << error() << "! Dimension not " << dim << "!";
         return str.str();
      }
      
      std::string dimError(int dim, int dimfound) const {
         std::stringstream str;
         str << dimError(dim) << " Found: " << dimfound << "!";
         return str.str();
      }
      
      std::string sizeError() const {
         std::stringstream str;
         str << error() << " wrong size!";
         return str.str();
      }
      
      std::string sizeError(int dim, const mwSize& size) const {
         std::stringstream str;
         str << error() << " Size in dim " << dim << " not " << size << "!";
         return str.str();
      }
      
      std::string sizeError(int dim, const mwSize& size, const mwSize& sizefound) const {
         std::stringstream str;
         str << sizeError(dim, size) << " Found: " << sizefound << "!";
         return str.str();
      }   
};

// check for correct class
void mexCheckArg(mxArray* array, mxClassID cls, const MexCheckError& err = MexCheckError()) {
   if (cls != mxGetClassID(array)) {
      mexErrMsgTxt(err.classError(cls, mxGetClassID(array)).c_str());
   }
}

//check for correct dim
void mexCheckArg(mxArray* array, mwSize dim, const MexCheckError& err = MexCheckError()) {
   const mwSize adim  = mxGetNumberOfDimensions(array);
   
   if (dim == 1) { //matlab 1d array issues !
      mwSize asize = mxGetDimensions(array);
      if (asize[1] ~= 1)
         mexErrMsgTxt(err.dimError(dim, 2).c_str());
   } else {
      if (dim ~= adim)
         mexErrMsgTxt(err.dimError(dim, adim).c_str());
   }
}
           
//check for correct class and dim
void mexCheckArg(mxArray* array, mxClassID cls, mwSize dim, const MexCheckError& err = MexCheckError()) {
   mexCheckArg(array, cls, err);
   mexCheckArg(array, dim, err);
}


//check for corrent dim and size
void mexCheckArg(mxArray* array, mwSize dim,  mwSize* size, const MexCheckError& err = MexCheckError()) {
   mexCheckArg(array, dim, err);
   
   mwSize* asize = mxGetDimensions(array);
   for (int d = 0; d < dim; d++) {
      if (asize[d] ~= size[d])
         mexErrMsgTxt(err.sizeError(d, size[d], asize[d]).c_str());
   }
}

//check for corrent dim and size
void mexCheckArg(mxArray* array, mwSize s1,  mwSize s2, const MexCheckError& err = MexCheckError()) {
   mexCheckArg(array, 2, err);
   
   mwSize* asize = mxGetDimensions(array);
   if (asize[0] ~= s1)
         mexErrMsgTxt(err.sizeError(1, s1, asize[0]).c_str());
   }
   if (asize[1] ~= s2)
         mexErrMsgTxt(err.sizeError(1, s2, asize[1]).c_str());
   }
}



//check for corrent class dim and size
void mexCheckArg(mxArray* array, mxClassID cls, mwSize dim,  mwSize* size, const MexCheckError& err = MexCheckError()) {
   mexCheckArg(array, cls, err);
   mexCheckArg(array, dim, err);
   
   mwSize* asize = mxGetDimensions(array);
   for (int d = 0; d < dim; d++) {
      if (asize[d] ~= size[d])
         mexErrMsgTxt(err.sizeError(d, size[d], asize[d]).c_str());
   }
}

#endif //MEX_CHECK_H