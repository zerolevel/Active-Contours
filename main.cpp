#include "activeContour.h"
#define __ERROR__ fprintf(stderr,"Main%d\n",__LINE__);
float alpha=0.6;
int iterations=200;
int main(int argc, char** argv)
{
	IplImage *img = 0, *mask=0,*mask2=0, *seg=0;
	
	if(argc < 3)
	{
		fprintf(stderr,"Usage: ma‭in <Image> <Mask>");
	}	
	
	
	img = cvLoadImage( argv[1] );
	
	cvNamedWindow( "segmented", CV_WINDOW_NORMAL );
	cvMoveWindow( "segmented", 100, 100);
	cvShowImage( "segmented", img);
	cvWaitKey(0); 
	
	mask = cvLoadImage( argv[2] );
	
	cvShowImage( "segmented", mask);
	cvWaitKey(0); 
	
	mask2 = cvCreateImage(cvGetSize(mask),img->depth,CV_8UC1);
	cvCvtColor(mask, mask2, CV_RGB2GRAY);
	cvReleaseImage( &mask);
	
	
	
	cvShowImage( "segmented", mask2);
	cvWaitKey(0); 
	
	fprintf(stderr, "hM:%d, wM:%d, dM:%d, hI:%d, wI:%d", mask2->height, mask2->width, mask2->depth, img->height,img->width );
	
	seg = region_seg( img, mask2, iterations, alpha, 1 );
	
	cvSaveImage( "Segmented.tiff", seg );	
	
	cvNamedWindow( "segmented", CV_WINDOW_NORMAL );
	cvMoveWindow( "segmented", 100, 100);
	cvShowImage( "segmented", seg);
	cvWaitKey(0); 


	return 0;
}

