/***************************************************************
 * Propagation of Seeds
 * using Relative Intensity Changes
 *
 * C. Kirst, The Rockefeller University 2014
 *
 * propagate seeds until spatial distiance or intensity change
 * become to large
 *
 ***************************************************************/

#include <math.h>
#include "mex.h"
#include <queue>
#include <vector>
#include <map>
#include <iostream>
using namespace std;

/* Input Arguments */
#define IM_IN              prhs[0]
#define SEEDS_IN           prhs[1]
#define MASK_IN            prhs[2]
#define CUTOFF_INTE_IN     prhs[3]
#define CUTOFF_DIST_IN     prhs[4]
#define RADIUS_IN          prhs[5]
#define INTENSITY_REF_IN   prhs[6]


/* Output Arguments */
#define SEEDS_OUT         plhs[0]
#define DISTANCES_OUT     plhs[1]


// weight spatial vs intensity differences
double alpha;

// image sizes and strides
unsigned int m, n, l, mn;


inline unsigned int id(unsigned int i, unsigned int j, unsigned int k) {
   return (k*mn + j*m + i);
}

inline unsigned int id2i(unsigned int id) {
   return (id % m)
}

inline unsigned int id2j(unsigned int id) {
   return ((id - id/mn * mn) / m)
}

inline unsigned int id2k(unsigned int id) {
   return (id / mn)
}

inline void id2ijk(unsigned int id, unsigned int& i, unsigned int& j, unsigned int& k) {
   k = id /nm;
   j = (id -k * mn)/ m;
   i = id % m;
   return;
}


// average intensity around a point
static double average_intensity(double * image, unsigned int i, unsigned int j, unsigned int k, int radius) {
   
   double intensity = 0.0;
   int delta_i, delta_j, delta_k;
   int norm = 0;
   int idi, idj, idk;
   
   for (delta_i = -radius; delta_i <= radius; delta_i++) {
      idj = i + delta_i;
      if (idi >= 0 && idi < m) {
         for (delta_j = -radius; delta_j <= radius; delta_j++) {
            idj = j + delta_j;
            if (idj >=0 && idj <= n) {
               for (delta_k = -radius; delta_k <= radius; delta_k++) {
                  idk = k + delta_k;
                  if (idk >= 0 && idk <= l) {
                     intensity += image[id(idi,idj,idk)];
                     norm++;
                  }
               }
            }
         }
      }
   }
   return (intensity/norm);
};



//Pixel + info on distances
class Pixel { 
public:
   double spatial_distance;
   double intensity_distance;
   unsigned int ijk;
   double label;
   
public:
   Pixel (double id, double sd, unsigned int inijk, double l) :
      intensity_distance(d), spatial_distance(sd), ijk(inijk), label(l) {
      }
      
   unsigned int i {
      return id2i(ijk);
   }
   
   unsigned int j {
      return id2j(ijk);
   }
   
   unsigned int k {
      return id2k(ijk);
   }
};

struct PixelCompare { 
   bool operator() (const Pixel& a, const Pixel& b) const {
      return alpha * (a.spatial_distance - b.spatial_distance) + (1-alpha) * (a.intensity_difference - b.intensity_distance) > 0; 
   }
};

typedef priority_queue<Pixel, vector<Pixel>, PixelCompare> PixelQueue;



















/* difference between two pixels */
double update_distances(const Pixel& p1, Pixel& p2, double *image, int radius)
{
   
   // calculate intensity change
   int delta_i, delta_j, delta_k;
   int norm = 0;
   int idi, idj, idk;
   
   double intensity_diff = 0.0;
   
   for (delta_i = -radius; delta_i <= radius; delta_i++) {
      idj = i2 + delta_i;
      if (idi >= 0 && idi < m) {
         for (delta_j = -radius; delta_j <= radius; delta_j++) {
            idj = j2 + delta_j;
            if (idj >=0 && idj <= n) {
               for (delta_k    
  
   //here is space for taking into account gradient image / gradient crossings form the label center to the new pixel etc....
  
   return ((alpha-1)*pixel_diff + lambda * space_dist);= -radius; delta_k <= radius; delta_k++) {
                  idk = k2 + delta_k;
                  if (idk >= 0 && idk <= l) {
                     intensity_diff += fabs(image[id(idi, idj, idk)] - ref_intensity); 
                     norm++;
                  }
               }
            }
         }
      }
   }
   
   intensity_diff *= 1.0/norm;
   p2.intensity_distance = p1.intensity_distance + intensity_diff;
   
   // spatial change
   double d1 = (p1.i - p2.i); double d2 = (p1.j - p2.j); double d3 = (p1.k - p2.k);
   double space_diff = sqrt(d1*d1 + d2*d2 + d3*d3);  
  
   p2.space_distance = p1.space_distance + space_diff;
};






static void push_neighbors_on_queue(PixelQueue &pq, 
        
        
        
        double spatial_dist,
                        double *image,
                        unsigned int i, unsigned int j, unsigned int k,
                        int radius, double ref_intensity,  double lambda, double cutoff_dist,
                        double label, double *seeds_out, mxLogical* mask_in)
{
   
   double d, sd;
  
   /* 26-connected */
   unsigned int di, dj, dk, idi, jdj, kdk, ii;
   
   for(di = -1; di <= 1; di++) {
      idi = i +di; 
      if (idi>=0 && idi < m) {
         for (dj = -1; dj <= 1; dj++) {
            jdj = j + dj;
            if (jdj >= 0 && jdj < n) {
               for (dk = -1; dk <= 1; dk++) {
                  kdk = k + dk;
                  ii = id(idi,jdj, kdk);
                  
                  if (kdk >= 0 && kdk < l && mask_in[ii] &&  0 == seeds_out[ii]) {
                     //only push if neighbours are within resonable distance
                     d = difference(image, i, j, k, idi, jdj, kdk, m, n, l, radius, ref_intensity, lambda, spatial_dist, sd);
                     //cout << "i,j,k=" << i << "," << j << "," << k << endl;
                     //cout << "di,j,k=" << idi << "," << jdj << "," << kdk << endl;
                     //cout << "d=" << d << " sp_dist=" << spatial_dist << " sd=" << sd << endl;
                     //cout << "ref_intensity=" << ref_intensity << endl;
                     if (d < cutoff_dist) {
                        pq.push(Pixel(d, sd, idi, jdj, kdk, label));
                     }
                  }
               }
            }
         }
      }
   }          
};

static void propagate(double *seeds_in, double *im_in, mxLogical *mask_in, 
                      double *seeds_out, 
                      double *dists,
                      int radius, 
                      double *center_intensities, unsigned int nlabel,
                      double lambda, double cutoff_distance) {

   PixelQueue pixel_queue;
   //map<double, double> center_intensities;  / for auto center_intensities 
  
  
   /* initialize dist to Inf, read seeds_in and wrtite out to seeds_out */
   for (unsigned int i = 0; i < m; i++) {
      for (unsigned int j = 0; j < n; j++) {
         for (unsigned int k = 0; k < l; k++) {
            unsigned int ii = id(i,j,k);
            dists[ii] = mxGetInf();            
            seeds_out[ii] = seeds_in[ii];
         }
      }
   }
  
   /* if the pixel is already labeled (i.e., labeled in seeds_in) and within a mask, 
    * then set dist to 0 and push its neighbors for propagation */
   for (unsigned int i = 0; i < m; i++) { 
      for (unsigned int j = 0; j < n; j++) {
         for (unsigned int k = 0; k < l; k++) {
            unsigned int ii = id(i,j,k);
            double label = seeds_in[ii];
            if ((label > 0) && (mask_in[ii])) {
               
               dists[ii] = 0.0;
               push_neighbors_on_queue(pixel_queue, 0.0, im_in, i, j, k, radius, center_intensities[(int) label], lambda, cutoff_distance, label, seeds_out, mask_in);
        
            }
         }
      }
   }

   while (! pixel_queue.empty()) {
      Pixel p = pixel_queue.top();
      pixel_queue.pop();
      //cout << "popped " << p.i << " " << p.j << " " << p.k << endl;

      unsigned int ii = p.ijk;
      if ((dists[ii] > p.distance) /* && (mask_in[IJ(p.i,p.j, p.k)])*/) {
         dists[ii] = p.distance;
         seeds_out[ii] = p.label;      
         push_neighbors_on_queue(pixel_queue, /* p.distance,*/ p.spatial_distance, im_in, p.i, p.j, p.k, 
                                 radius, center_intensities[(int) p.label], lambda, cutoff_distance, p.label, seeds_out, mask_in);
      }
   }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ) { 
   
   double *seeds_in, *im_in; 
   mxLogical *mask_in;
   double *seeds_out, *dists;
   double *lambda;
   double *cutoff_distance;
   int radius /*, radius_center*/;    
   unsigned int m, n, l; 
   
   double *intensity_refs;
   unsigned int nlabel;
   mwSize  ndim;
    
    
   /* Check for proper number of arguments */  
   if (nrhs != 7) { 
      mexErrMsgTxt("7 input arguments required."); 
   } else if (nlhs !=1 && nlhs !=2) {
      mexErrMsgTxt("The number of output arguments should be 1-2."); 
   } 
    
   /* Size of Image */
   ndim = mxGetNumberOfDimensions(IM_IN);
   if (3 != ndim){
      mexErrMsgTxt("input image must be 3d grayscale image");
   }
   
   const mwSize* size = mxGetDimensions(IM_IN);
   m = size[0];
   n = size[1];
   l = size[2];
   
   //cout << m << " " << n << " " << l << endl;


   /* Size of Labeled Image */
   ndim = mxGetNumberOfDimensions(SEEDS_IN);
   size = mxGetDimensions(SEEDS_IN);
   
   /* Check for the correct number of indices  */
   if (3 != ndim){
      mexErrMsgTxt("input label image must be 3d grayscale image");
   }

   if ((m != size[0]) || (n != size[1]) || ( l != size[2])) {
      mexErrMsgTxt("sizes of image and labels do not agree");
   }

   if (! mxIsDouble(IM_IN)) {
      mexErrMsgTxt("First argument must be a double array.");
   }
   if (! mxIsDouble(SEEDS_IN)) {
      mexErrMsgTxt("Second argument must be a double array.");
   }
   if (! mxIsLogical(MASK_IN)) {
      mexErrMsgTxt("Third argument must be a logical array.");
   }

   /* Create matrices for the return arguments */    
    
   SEEDS_OUT = mxCreateNumericArray(ndim, size, mxDOUBLE_CLASS, mxREAL); 
   DISTANCES_OUT = mxCreateNumericArray(ndim, size, mxDOUBLE_CLASS, mxREAL);
    
   /* Assign pointers to the various parameters */ 
   seeds_in = mxGetPr(SEEDS_IN);
   im_in = mxGetPr(IM_IN);
   mask_in = mxGetLogicals(MASK_IN);
   lambda = mxGetPr(LAMBDA_IN);
   cutoff_distance = mxGetPr(CUTOFF_DIST_IN);
   intensity_refs = mxGetPr(INTENSITY_REF_IN);
   nlabel = mxGetM(INTENSITY_REF_IN);

   double * dptr = mxGetPr(RADIUS_IN);
   radius = (int) (*dptr);
   //dptr = mxGetPr(RADIUS_CENTER_IN);
   //radius_center = (int) (*dptr);
    
   seeds_out = mxGetPr(SEEDS_OUT);
   dists = mxGetPr(DISTANCES_OUT);

   //cout << "lambda = " << *lambda << " cutoff_difference = " << *cutoff_distance << endl;
   //cout << "radius = " << radius << " nlabel=" << nlabel << " (m,n,l)=(" << m << "," << n << "," << l << ")" << endl;
   
   /* Do the actual computations in a subroutine */
   propagate(seeds_in, im_in, mask_in, seeds_out, dists, m, n, l, radius, intensity_refs, nlabel, *lambda, *cutoff_distance); 

   if (nlhs == 1) {
      mxDestroyArray(DISTANCES_OUT);
   }      

   return;
};
