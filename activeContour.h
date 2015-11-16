#ifndef __activeContour__
#define __activeContour__


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cv.h>
#include <highgui.h>
#include <pthread.h>

#define eps 0.00000000000000000000002

void		cvGetCurvature( IplImage* phi, IplImage* idx, IplImage *dst );
void 		cvMask2phi( IplImage* init_a, IplImage* phi);
IplImage* 	region_seg( IplImage* I,IplImage* init_mask,int max_its,float alpha, int status );
void 		cvSussmanSign( IplImage *D, IplImage *S );
void 		cvSussman( IplImage *D, IplImage *dst, float dt );

#endif
