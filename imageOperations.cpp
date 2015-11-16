#include "imageOperations.h"
#define __ERROR__ fprintf(stderr,"\t\tIO:%d\n",__LINE__);
void getMatrix( IplImage* img, POINTS2d *points, CvMat* valMat )
{
	
	int currentPoint = 0;
	int widthStep = (img->widthStep)/4;
	int baseIndex = -1*widthStep;
	int width = img->width;
	int height = img->height;
	float* imgdata = (float*) img->imageData;
	int j;
	int* x = points->x;
	int* y = points->y;
	float * vals = valMat->data.fl; 
	for( int i =0; i<height; i++ )
	{	
		baseIndex +=widthStep; 
		 
		while( (i+1)  == y[currentPoint] )
		{	
			j = x[currentPoint]-1;
			vals[currentPoint] = imgdata[ baseIndex+j ];
			currentPoint++;
			if(currentPoint == points->num)	return;
		}
	}		
}

void cvShift( IplImage *src, IplImage*dst, int option )
{

	int widthStep = (src->widthStep)/4;
	int width = src->width;
	int height = src->height;
	int nChannels = src->nChannels;
	int baseIndex = -1*widthStep;
	float *data = (float *)src->imageData;
	float *data_new = (float *)dst->imageData;
	
	int h,w,baseIndex1;
	
	baseIndex = -1*widthStep;

	switch(option)
	{
	case 1:
	//Shift_U
	for(  h = 0 ; h < height-1 ; h++ )
	{      
		baseIndex += widthStep;
		baseIndex1 = baseIndex + widthStep ;
		for(  w = 0 ; w < width ; w++ )
		{
			data_new[ baseIndex + w ] =   data[ baseIndex1 + w ];
		}

	}
	baseIndex += widthStep;
	for(  w = 0 ; w < width ; w++ )
	{
		data_new[ baseIndex + w ] =   data[ baseIndex + w];
	}
	break;
	
	case 2:
	//Shift_D	
	for(  h = 0 ; h < height ; h++ )
	{      
		baseIndex += widthStep;
		if( h != 0 )
		{
			 baseIndex1 = baseIndex - widthStep;
			for(  w = 0 ; w < width ; w++ )
			{
				data_new[ baseIndex + w ] =   data[ baseIndex1 + w];
			}
		}
		else
			for(  w = 0 ; w < width ; w++ )
			{
				data_new[ baseIndex + w ] =   data[ baseIndex + w ];
			}
	}
	break;
	
	case 3:
	//shift_R
	for(  h = 0 ; h < height ; h++ )
	{      
		baseIndex += widthStep;
		for(  w = 0 ; w < width ; w++ )
		{
			if ( w != 0  )	data_new[ baseIndex + w ] =   data[ baseIndex + w - 1];
			else data_new[ baseIndex + w ] =   data[ baseIndex + w ];
		}
	}
	break;
	
	case 4:
	//Shift_L

	for(  h = 0 ; h < height ; h++ )
	{      
		baseIndex += widthStep;
		for(  w = 0 ; w < width ; w++ )
		{
			if( w != width-1 ) 
				data_new[ baseIndex + w ] =   data[ baseIndex + w + 1];
			else data_new[ baseIndex + w ] =   data[ baseIndex + w ];
		}
	}
	break;
	}		
}

void cvmySqrt( IplImage *src )
{
	int widthStep = (src->widthStep)/4;
	int width = src->width;
	int height = src->height;
	int baseIndex = -1*widthStep;
	float *data = (float *)src->imageData;
	int h,w;
	for(  h = 0 ; h < height ; h++ )
	{      
		baseIndex += widthStep;
		for(  w = 0 ; w < width ; w++ )
		{	data[ baseIndex + w ] =   sqrt( data[ baseIndex + w ] );
		}
	}
} 

void cvAndDiff( IplImage *src, CvScalar val, IplImage* Index )
{
	cvSet( src, val, Index );
	return;
}

void cvMyOverLapping( IplImage* src1, IplImage* src2, IplImage* src3, IplImage* src4, IplImage* Index, IplImage* dst )
{
	uchar* idx = ( uchar* )Index->imageData;
	
	int widthStep1 = Index->widthStep;
	int widthStep = src1->widthStep/4;
	int baseIndex1 = -1*widthStep1;
	int baseIndex = -1*widthStep;
	int width = Index->width;
	int height = Index->height;
	
	float* src1Data = ( float* )src1->imageData;
	float* src2Data = ( float* )src2->imageData;
	float* src3Data = ( float* )src3->imageData;
	float* src4Data = ( float* )src4->imageData;
	float* dstData  =  ( float* )dst->imageData;
	double val1, val2, val3;
	int i, j;
	int k =0;
	uchar a = (uchar) 255;
	for( i =0; i<height; i++)
	{
		baseIndex1 +=widthStep1;
		baseIndex +=widthStep;
		for( j = 0 ; j<width; j++ )
		{	k++;
			if(idx[ baseIndex1+j ] == a )
			{
				val1 = src1Data[ baseIndex+j ] * src1Data[ baseIndex+j ];
				val2 = src2Data[ baseIndex+j ] * src2Data[ baseIndex+j ];
				val1 = myMax( val1, val2 );	
				val3 = src3Data[ baseIndex+j ] * src3Data[ baseIndex+j ];
				val2 = src4Data[ baseIndex+j ] * src4Data[ baseIndex+j ];
				val2 = myMax( val3, val2 );	
				val3 = sqrt(val1+val2) -1;
				//printf("%lf\n",val3);
				dstData[ baseIndex+j ] =val3;
			}
		}
		//printf("\n");	
	}
	//printf("%d", k);
}

void cvAndDiffS( IplImage *src, IplImage *dst, IplImage* Index )
{
	//BOTH Src, AND Dst ARE OF 32 BIT TYPE iMAGES.
	uchar* idx = ( uchar* )Index->imageData;
	int widthStep1 = Index->widthStep;
	int widthStep = (src->widthStep)/4;

	int baseIndex1 = -1*widthStep1;
	int baseIndex = -1*widthStep;
	int width = Index->width;
	int height = Index->height;
	
	float* srcData = (float*)src->imageData;
	float* dstData = (float*)dst->imageData;
	
	int i,j;
	uchar a = (uchar) 0xff;
	for( i =0; i<height; i++)
	{
		
		baseIndex +=widthStep;
		baseIndex1 +=widthStep1;
		for( j = 0 ; j<width; j++ )
		{	
			if( a == idx[ baseIndex1+j ] )
				dstData[ baseIndex+j ] = srcData[ baseIndex+j ];
		}	
	}

}

IplImage *cvIm2GrayDouble(IplImage* Iimage)
{	IplImage *IGray = NULL, *IGray32 = NULL;
	CvSize sizeImage = cvGetSize(Iimage);
	if(Iimage->nChannels!=1)
	{	IGray = cvCreateImage( sizeImage, Iimage->depth, 1 );	
		cvZero( IGray );
		cvCvtColor( Iimage, IGray, CV_RGB2GRAY );
		cvReleaseImage( &Iimage );
	}
	else 
	{ 
		IGray = Iimage; 
	}
				
	if(IGray->depth!=IPL_DEPTH_32F)
	{	
		IGray32  = cvCreateImage( sizeImage, IPL_DEPTH_32F, 1 );
		cvZero( IGray32 );
		cvConvertScale( IGray, IGray32 );
		cvReleaseImage( &IGray );
		
	}	
	else 
	{
		IGray32 = IGray;
	}
	return IGray32;	 
}

double myMax( double data1, double data2 )
{
	if(data1>data2)
		return data1;
	else return data2;
}


struct POINTS2d* cvImage2Sub( IplImage* img)  //works perfectly alright
{
	// img is a binary image with 8bit depth with 0xff for 1 and 0x00 f0r 0;
	CvScalar lengthCVS;
	lengthCVS = cvSum( img );
	int  length = lengthCVS.val[ 0 ]/255;
	struct POINTS2d* points=NULL;
	if(( points = ( POINTS2d*) malloc( sizeof(  POINTS2d ))) == NULL )
		exit( 1 );
	points->num = length;
	if(( points->x = ( int*) malloc( sizeof(  int ) * (length + 20) )) == NULL )
		exit( 1 );
		
	if(( points->y = ( int*) malloc( sizeof(  int ) * (length + 20) )) == NULL )
		exit( 1 );
	
	int baseIndex = -1*img->widthStep;
	int width = img->width;
	int height = img->height;
	int widthStep = img->widthStep;
	int i,j;
	uchar *imgData = (uchar *)img->imageData;  
	int currentPoint = 0;
	unsigned char a = 0xff;
	for( int i =0; i<height; i++ )
	{
		baseIndex +=widthStep; 
		for( int j = 0 ; j<width; j++ )
		{		
			if( imgData[ baseIndex+j ] == a )
			{
				points->y[ currentPoint ] = i+1;
				points->x[ currentPoint ] = j+1;
				currentPoint++;
			}
		}	
	} 
	return points;
}

void cvEquateIndex( IplImage *src, IplImage* equate, IplImage* index )
{
	/*
		when index(i,j) is 255 it puts value of equate into src
	*/	
	int widthStep = (src->widthStep)/4;
	int baseIndex  = -1*widthStep;
	int widthStep1 = index->widthStep;
	int baseIndex1 = -1*widthStep1;
	int width = src->width;
	int height = src->height;
	
	uchar* idx = (uchar*)index->imageData;
	float* srcData = (float *)src->imageData; 
	float* equateData = (float *) equate->imageData; 
	
	uchar a = 0xff;
	for(int i =0;i<height;i++)
	{	
		baseIndex  += widthStep;
		baseIndex1 += widthStep1;
		for( int j = 0 ; j<width; j++ )
		{		
			if( a == idx[baseIndex1+j] )
			{
				srcData[baseIndex+j] = equateData[baseIndex+j];
			}
			//else srcData[baseIndex+j] = srcData[baseIndex+j];
		}	
	}
}

struct POINTS2d* cvCreatePoints( int* y,  int*x, int num)
{
	struct POINTS2d *points = NULL;
	
	if((points = (struct POINTS2d* )malloc(sizeof( struct POINTS2d )) )== NULL) exit (1);
	
	points->num = num;
	points->x = x;
	points->y =y;
	
	return points;
}

void cvGetImage( CvMat* src, POINTS2d* idx, IplImage* dst )
{
	cvZero(dst);
	int currentPoint = 0;
	int widthStep = (dst->widthStep)/4;
	int baseIndex = -1*widthStep;
	int width = dst->width;
	int height = dst->height;
	cvZero( dst );
	float* dstdata = (float*) dst->imageData;
	int j;
	int* x = idx->x;
	int* y = idx->y;
	int totalPoints = idx->num;
	float * vals = src->data.fl; 
	//fprintf(stderr, "checking Called %d\n", points.num );	
	for( int i =0; i<height; i++ )
	{	
		baseIndex +=widthStep; 
		
		while( (i+1)==y[currentPoint] )
		{	
			j = x[currentPoint]-1;
			dstdata[ baseIndex + j ] = vals[currentPoint]; 
			currentPoint++;
			if(currentPoint == totalPoints)	return;
		}
	}		
}
