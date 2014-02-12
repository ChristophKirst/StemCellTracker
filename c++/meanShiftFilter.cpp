/*******************************************************************
 * meanShiftFilter.C
 *
 * Mean shift filtering of an 2d image
 *******************************************************************/

#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	IMAGE_IN     prhs[0]
#define	RADIUS_IN    prhs[1]
#define  DISTANCE_IN  prhs[2]


/* Output Arguments */

#define	IMAGE_OUT    plhs[0]


static void meanShiftFilter(double* img, int n, int m, double radius, double distance, double* out)
{
   for(int i = 0; i < n; i++) {
       for(int j = 0; i < n; i++) {
      
   
       } 
   } 
   return;
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    double *image_in;
    size_t  m,n;
    
    double *image_out;
    
    /* Check for proper number of arguments */
    
    if (nrhs != 3) { 
	    mexErrMsgIdAndTxt( "MATLAB:meanShiftilter:invalidNumInputs",
                "Three input arguments required."); 
    }
    
    m = mxGetM(IMAGE_IN); 
    n = mxGetN(IMAGE_IN);
    
    
    /* Create a matrix for the return argument */ 
    IMAGE_OUT = mxCreateDoubleMatrix( (mwSize)m, (mwSize)n, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    image_out = mxGetPr(IMAGE_OUT);
    
    image_in = mxGetPr(IMAGE_IN); 
        
    /* Do the actual computations in a subroutine */
    meanShiftFilter(image_in, RADIUS_IN, DISTANCE_IN, image_out); 
    return;
    
}





/****************************************************************************************\
*                                         Meanshift                                      *
\****************************************************************************************/

CV_IMPL void
cvPyrMeanShiftFiltering( const CvArr* srcarr, CvArr* dstarr,
                         double sp0, double sr, int max_level,
                         CvTermCriteria termcrit )
{
    const int cn = 3;
    const int MAX_LEVELS = 8;

    if( (unsigned)max_level > (unsigned)MAX_LEVELS )
        CV_Error( CV_StsOutOfRange, "The number of pyramid levels is too large or negative" );

    std::vector<cv::Mat> src_pyramid(max_level+1);
    std::vector<cv::Mat> dst_pyramid(max_level+1);
    cv::Mat mask0;
    int i, j, level;
    //uchar* submask = 0;

    #define cdiff(ofs0) (tab[c0-dptr[ofs0]+255] + \
        tab[c1-dptr[(ofs0)+1]+255] + tab[c2-dptr[(ofs0)+2]+255] >= isr22)

    double sr2 = sr * sr;
    int isr2 = cvRound(sr2), isr22 = MAX(isr2,16);
    int tab[768];
    cv::Mat src0 = cv::cvarrToMat(srcarr);
    cv::Mat dst0 = cv::cvarrToMat(dstarr);

    if( src0.type() != CV_8UC3 )
        CV_Error( CV_StsUnsupportedFormat, "Only 8-bit, 3-channel images are supported" );

    if( src0.type() != dst0.type() )
        CV_Error( CV_StsUnmatchedFormats, "The input and output images must have the same type" );

    if( src0.size() != dst0.size() )
        CV_Error( CV_StsUnmatchedSizes, "The input and output images must have the same size" );

    if( !(termcrit.type & CV_TERMCRIT_ITER) )
        termcrit.max_iter = 5;
    termcrit.max_iter = MAX(termcrit.max_iter,1);
    termcrit.max_iter = MIN(termcrit.max_iter,100);
    if( !(termcrit.type & CV_TERMCRIT_EPS) )
        termcrit.epsilon = 1.f;
    termcrit.epsilon = MAX(termcrit.epsilon, 0.f);

    for( i = 0; i < 768; i++ )
        tab[i] = (i - 255)*(i - 255);

    // 1. construct pyramid
    src_pyramid[0] = src0;
    dst_pyramid[0] = dst0;
    for( level = 1; level <= max_level; level++ )
    {
        src_pyramid[level].create( (src_pyramid[level-1].rows+1)/2,
                        (src_pyramid[level-1].cols+1)/2, src_pyramid[level-1].type() );
        dst_pyramid[level].create( src_pyramid[level].rows,
                        src_pyramid[level].cols, src_pyramid[level].type() );
        cv::pyrDown( src_pyramid[level-1], src_pyramid[level], src_pyramid[level].size() );
        //CV_CALL( cvResize( src_pyramid[level-1], src_pyramid[level], CV_INTER_AREA ));
    }

    mask0.create(src0.rows, src0.cols, CV_8UC1);
    //CV_CALL( submask = (uchar*)cvAlloc( (sp+2)*(sp+2) ));

    // 2. apply meanshift, starting from the pyramid top (i.e. the smallest layer)
    for( level = max_level; level >= 0; level-- )
    {
        cv::Mat src = src_pyramid[level];
        cv::Size size = src.size();
        uchar* sptr = src.data;
        int sstep = (int)src.step;
        uchar* mask = 0;
        int mstep = 0;
        uchar* dptr;
        int dstep;
        float sp = (float)(sp0 / (1 << level));
        sp = MAX( sp, 1 );

        if( level < max_level )
        {
            cv::Size size1 = dst_pyramid[level+1].size();
            cv::Mat m( size.height, size.width, CV_8UC1, mask0.data );
            dstep = (int)dst_pyramid[level+1].step;
            dptr = dst_pyramid[level+1].data + dstep + cn;
            mstep = (int)m.step;
            mask = m.data + mstep;
            //cvResize( dst_pyramid[level+1], dst_pyramid[level], CV_INTER_CUBIC );
            cv::pyrUp( dst_pyramid[level+1], dst_pyramid[level], dst_pyramid[level].size() );
            m.setTo(cv::Scalar::all(0));

            for( i = 1; i < size1.height-1; i++, dptr += dstep - (size1.width-2)*3, mask += mstep*2 )
            {
                for( j = 1; j < size1.width-1; j++, dptr += cn )
                {
                    int c0 = dptr[0], c1 = dptr[1], c2 = dptr[2];
                    mask[j*2 - 1] = cdiff(-3) || cdiff(3) || cdiff(-dstep-3) || cdiff(-dstep) ||
                        cdiff(-dstep+3) || cdiff(dstep-3) || cdiff(dstep) || cdiff(dstep+3);
                }
            }

            cv::dilate( m, m, cv::Mat() );
            mask = m.data;
        }

        dptr = dst_pyramid[level].data;
        dstep = (int)dst_pyramid[level].step;

        for( i = 0; i < size.height; i++, sptr += sstep - size.width*3,
                                          dptr += dstep - size.width*3,
                                          mask += mstep )
        {
            for( j = 0; j < size.width; j++, sptr += 3, dptr += 3 )
            {
                int x0 = j, y0 = i, x1, y1, iter;
                int c0, c1, c2;

                if( mask && !mask[j] )
                    continue;

                c0 = sptr[0], c1 = sptr[1], c2 = sptr[2];

                // iterate meanshift procedure
                for( iter = 0; iter < termcrit.max_iter; iter++ )
                {
                    uchar* ptr;
                    int x, y, count = 0;
                    int minx, miny, maxx, maxy;
                    int s0 = 0, s1 = 0, s2 = 0, sx = 0, sy = 0;
                    double icount;
                    int stop_flag;

                    //mean shift: process pixels in window (p-sigmaSp)x(p+sigmaSp)
                    minx = cvRound(x0 - sp); minx = MAX(minx, 0);
                    miny = cvRound(y0 - sp); miny = MAX(miny, 0);
                    maxx = cvRound(x0 + sp); maxx = MIN(maxx, size.width-1);
                    maxy = cvRound(y0 + sp); maxy = MIN(maxy, size.height-1);
                    ptr = sptr + (miny - i)*sstep + (minx - j)*3;

                    for( y = miny; y <= maxy; y++, ptr += sstep - (maxx-minx+1)*3 )
                    {
                        int row_count = 0;
                        x = minx;
                        #if CV_ENABLE_UNROLLED
                        for( ; x + 3 <= maxx; x += 4, ptr += 12 )
                        {
                            int t0 = ptr[0], t1 = ptr[1], t2 = ptr[2];
                            if( tab[t0-c0+255] + tab[t1-c1+255] + tab[t2-c2+255] <= isr2 )
                            {
                                s0 += t0; s1 += t1; s2 += t2;
                                sx += x; row_count++;
                            }
                            t0 = ptr[3], t1 = ptr[4], t2 = ptr[5];
                            if( tab[t0-c0+255] + tab[t1-c1+255] + tab[t2-c2+255] <= isr2 )
                            {
                                s0 += t0; s1 += t1; s2 += t2;
                                sx += x+1; row_count++;
                            }
                            t0 = ptr[6], t1 = ptr[7], t2 = ptr[8];
                            if( tab[t0-c0+255] + tab[t1-c1+255] + tab[t2-c2+255] <= isr2 )
                            {
                                s0 += t0; s1 += t1; s2 += t2;
                                sx += x+2; row_count++;
                            }
                            t0 = ptr[9], t1 = ptr[10], t2 = ptr[11];
                            if( tab[t0-c0+255] + tab[t1-c1+255] + tab[t2-c2+255] <= isr2 )
                            {
                                s0 += t0; s1 += t1; s2 += t2;
                                sx += x+3; row_count++;
                            }
                        }
                        #endif
                        for( ; x <= maxx; x++, ptr += 3 )
                        {
                            int t0 = ptr[0], t1 = ptr[1], t2 = ptr[2];
                            if( tab[t0-c0+255] + tab[t1-c1+255] + tab[t2-c2+255] <= isr2 )
                            {
                                s0 += t0; s1 += t1; s2 += t2;
                                sx += x; row_count++;
                            }
                        }
                        count += row_count;
                        sy += y*row_count;
                    }

                    if( count == 0 )
                        break;

                    icount = 1./count;
                    x1 = cvRound(sx*icount);
                    y1 = cvRound(sy*icount);
                    s0 = cvRound(s0*icount);
                    s1 = cvRound(s1*icount);
                    s2 = cvRound(s2*icount);

                    stop_flag = (x0 == x1 && y0 == y1) || abs(x1-x0) + abs(y1-y0) +
                        tab[s0 - c0 + 255] + tab[s1 - c1 + 255] +
                        tab[s2 - c2 + 255] <= termcrit.epsilon;

                    x0 = x1; y0 = y1;
                    c0 = s0; c1 = s1; c2 = s2;

                    if( stop_flag )
                        break;
                }

                dptr[0] = (uchar)c0;
                dptr[1] = (uchar)c1;
                dptr[2] = (uchar)c2;
            }
        }
    }
}

void cv::pyrMeanShiftFiltering( InputArray _src, OutputArray _dst,
                                double sp, double sr, int maxLevel,
                                TermCriteria termcrit )
{
    Mat src = _src.getMat();

    if( src.empty() )
        return;

    _dst.create( src.size(), src.type() );
    CvMat c_src = src, c_dst = _dst.getMat();
    cvPyrMeanShiftFiltering( &c_src, &c_dst, sp, sr, maxLevel, termcrit );
}

