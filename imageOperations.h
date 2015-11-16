#ifndef __imageOperations__
#define __imageOperations__
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cv.h>
#include <highgui.h>

struct POINTS2d
{
	int num;
	int* x;
	int* y;
};

void 		 getMatrix( IplImage* img, POINTS2d* points, CvMat* valMat );
void 		 cvShift( IplImage *src, IplImage*dst, int option );
void 		 cvmySqrt( IplImage *src );
void 		 cvAndDiff( IplImage *src, CvScalar val, IplImage* Index );
void 		 cvMyOverLapping( IplImage* src1, IplImage* src2, IplImage* src3, IplImage* src4, IplImage* Index, IplImage* dst );
void 		 cvAndDiffS( IplImage *src, IplImage *dst, IplImage* Index );
IplImage* 	 cvIm2GrayDouble( IplImage* I_image );
double 		 myMax( double data1, double data2 );
struct POINTS2d* cvImage2Sub( IplImage* img);
void 		 cvEquateIndex( IplImage *src, IplImage* equate, IplImage* index );
struct POINTS2d *cvCreatePoints(int* x, int*y, int num);
void 		 cvGetImage( CvMat* src, POINTS2d* idx, IplImage* dst );



#endif
