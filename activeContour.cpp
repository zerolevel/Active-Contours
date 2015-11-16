#include "activeContour.h"
#include "imageOperations.h"
#include "multiThread.h"
#define __ERROR__ fprintf(stderr,"\tAC:%d\n",__LINE__);

IplImage* region_seg( IplImage* I, IplImage* init_mask, int max_its, float alpha, int status )
{
	CvSize sizeImage = cvGetSize( I );
	
	IplImage* I_gray_double = cvIm2GrayDouble(I);  //releases image inside only cvReleaseImage( &I );
		
	
	I = I_gray_double;
	cvNamedWindow( "SEGMENTED IMAGE", CV_WINDOW_NORMAL );
	cvMoveWindow( "SEGMENTED IMAGE", 100, 100);


	IplImage* seg = cvCreateImage( sizeImage, IPL_DEPTH_8U, 1 );
	cvZero( seg );
	IplImage* phi =  cvCreateImage( sizeImage, IPL_DEPTH_32F, 1 );	
	cvZero( phi );
	cvMask2phi( init_mask, phi ); 
	IplImage *idx1 = cvCreateImage( sizeImage, IPL_DEPTH_8U, 1 );
	IplImage *idx2 = cvCreateImage( sizeImage, IPL_DEPTH_8U, 1 );
	IplImage *idx = cvCreateImage( sizeImage, IPL_DEPTH_8U, 1 );
	IplImage* U = cvCreateImage( sizeImage, phi->depth, 1);
	IplImage* f2 = cvCreateImage( sizeImage, IPL_DEPTH_32F, 1 );
	IplImage* f3 = cvCreateImage( sizeImage, IPL_DEPTH_32F, 1 );
	IplImage* f4 = cvCreateImage( sizeImage, IPL_DEPTH_32F, 1 );
	IplImage* V = cvCreateImage( sizeImage, phi->depth, 1);
	IplImage* upts = idx1;
	IplImage* vpts = idx2;
	IplImage* f = V; //using it so that neednot to be created again and again
	IplImage* f1 = U;
	IplImage* dphidt = f3;
	IplImage* curvature = f2;
	void** A1 = new void* [5];
	float uVal, vVal;
	A1[0] = (void*)f;
	A1[1] = (void*)f1;
	A1[2] = (void*)f2;
	A1[3] = (void*)f3;
	A1[4] = (void*)f4;
	A1[5] = (void*)idx;
	A1[6] = (void*)&uVal;
	A1[7] = (void*)&vVal;
	for( int i=0; i<max_its; i++)
	{	
		__ERROR__
		//idx = find(phi <= 1.2 & phi >= -1.2);
		//fprintf( stderr,"iterations = %d\n", i+1 );
		cvZero(idx1);
		__ERROR__
		cvZero(idx2);
		cvCmpS( phi, 1.2, idx1, CV_CMP_LE );
		cvCmpS( phi, -1.2, idx2, CV_CMP_GE );
		__ERROR__
		
		
		cvZero(idx);
		cvAnd( idx1, idx2, idx, NULL );
		__ERROR__
		
		//writeImageToFile(idx,"Idx.txt", 0);
		//cvReleaseImage( &idx1 );
		//cvReleaseImage( &idx2 );
		//writeImageToFile( phi,"phi.txt", 1);
		//writeImageToFile( idx,"Idx.txt", 0);
		// upts = find(phi<=0);     
		
		//puts( "phi" );
		//cvPrintParameters( phi );
		//puts( "idx" );
		//cvPrintParameters( idx );
		
		 //Instead of Releasing;
		cvZero( upts );
		cvCmpS( phi, 0, upts, CV_CMP_LE );
		
		//vpts = find(phi>0);   
		 //Instead of Releasing;
		cvZero( vpts );
		cvCmpS( phi, 0, vpts, CV_CMP_GT );
__ERROR__		
		// u = sum(I(upts))/(length(upts)+eps);
		
		cvZero( U );
		cvAndDiffS( I, U, upts );
__ERROR__		
		
		cvZero( V );
		cvAndDiffS( I, V, vpts );
__ERROR__		
		CvScalar USum = cvSum( U ); //U not needed after this
		CvScalar VSum = cvSum( V ); //V not needed after this
__ERROR__		
				
		CvScalar uOnes = cvSum( upts );
		   //upts not needed after this
		CvScalar vOnes = cvSum( vpts );
		  //vpts not needed after this
__ERROR__		
		float uSize = uOnes.val[0]/0xff;
		float vSize = vOnes.val[0]/0xff;
		
		uVal = USum.val[0]/(uSize + eps);
		vVal = VSum.val[0]/(vSize + eps);
__ERROR__	
		//fprintf( stderr ," %f, %f \n", uVal, vVal );
		
		//F = ( I(idx)- u ).^2-( I(idx)-v ).^2; 
		//F = ( f - u ).^2 - ( f - v ).^2;
		//F = 	f1.^2 - f2.^2
		//F =   f1 - f2
		//f3
		
		cvZero( f );
		cvAndDiffS( I, f, idx ); 
		//fprintf(stderr, "RG: %d\n", __LINE__);
__ERROR__		
		//void** A1 = (void**) malloc(5*sizeof(void*));
		
		
		//fprintf(stderr, "RG: %d\n", __LINE__);
		//fprintf(stderr, "RG: %d\n", __LINE__);
		
		multiThread1( A1 ); //Candidate to be run on parallel thread
__ERROR__		
		//fprintf(stderr, "RG: %d\n", __LINE__);
	//	cvZero( f1 );		
	//	cvZero( f4 );
	//	cvAndDiff( f4, cvScalar(uVal), idx ); 
	//	cvSub(f, f4, f1); //f1=f-f4
	//	cvPow( f1, f1, 2.0 );
	//
	//	cvZero( f2 );
	//	cvZero( f3 );
	//	cvAndDiff( f3, cvScalar(vVal), idx ); 
	//	cvSub(f, f3, f2);  //f2=f-f3
	//	cvPow( f2, f2, 2.0 );
		
		//f4 and f3 are not used after this
	//	cvSub( f1, f2, f );
			
		
		//curvature = get_curvature(phi,idx);
		cvZero( curvature );
__ERROR__		
		cvGetCurvature( phi, idx, curvature);
__ERROR__
		double minVal, maxVal; 
				
		// dphidt = F./max(abs(F)) + alpha*curvature;
		cvMinMaxLoc( f, &minVal, &maxVal, NULL, NULL, NULL );
		maxVal = myMax( abs(minVal), abs(maxVal) );
__ERROR__
		cvScale( f, f, 1/(maxVal) );
__ERROR__
		cvScale( curvature, curvature, alpha );
__ERROR__	
		cvZero( dphidt );
		
		cvAdd(f, curvature, dphidt);
__ERROR__
		cvMinMaxLoc( dphidt, &minVal, &maxVal, NULL, NULL, NULL );
__ERROR__
		float dt = 0.45/( maxVal + eps);
		
		//phi(idx) = phi(idx) + dt.*dphidt;
		cvZero( f1 );
		cvScale( dphidt, dphidt, dt );
__ERROR__		
		cvAndDiffS( phi, f1, idx );
__ERROR__
		cvAdd( f1, dphidt, f1 );
__ERROR__
		
		cvEquateIndex( phi, f1, idx );
__ERROR__			
		cvSussman( phi, phi, 0.5 );
__ERROR__	
		cvZero( seg ); 
__ERROR__		
		cvCmpS( phi, 0.00, seg, CV_CMP_LE );
__ERROR__		
		fprintf(stderr, "Iternation%d\n", i);
__ERROR__		
		if(status==1)
		{
			cvShowImage( "SEGMENTED IMAGE", seg );
			cvWaitKey(0);
		}	
		
	} 
		//cvShowImage( "SEGMENTED IMAGE", seg );
		//cvWaitKey(0);
		//cvShowImage( "SEGMENTED IMAGE", seg );
		//cvWaitKey(10);
	//cvCmpS( phi, 0, seg, CV_CMP_LE );
	//writeImageToFile(phi,"phi.txt", 1);
	//cvAdd( idx1, idx2, seg ); 
	/*
		cvReleaseImage( &vpts );
		cvReleaseImage( &upts );
		cvReleaseImage( &idx );
		cvReleaseImage( &f1  );
		cvReleaseImage( &f  );
		cvReleaseImage( &f4  );
		cvReleaseImage( &curvature  );
		cvReleaseImage( &dphidt );
	*/
	return seg;
}

void cvGetCurvature( IplImage* phi, IplImage* idx, IplImage *dst )
{	
	struct POINTS2d* points = cvImage2Sub( idx );
__ERROR__
	int *x = points->x;
	int *y = points->y;
	
__ERROR__	
	
	float length = points->num;
__ERROR__	
	CvMat* A[22];
__ERROR__	
	fprintf(stderr,"%d\n",length);
	for (int q = 0; q<22 ;q++)
	{
		A[q] = cvCreateMat( 1,length, CV_32FC1 );
		cvZero(A[q]);
	}
__ERROR__
	/*%-- get subscripts of neighbors
	ym1 = y-1; xm1 = x-1; yp1 = y+1; xp1 = x+1;
	
    	%-- bounds checking  
    	ym1(ym1<1) = 1; xm1(xm1<1) = 1;              
    	yp1(yp1>dimy)=dimy; xp1(xp1>dimx) = dimx;*/
    	
    	int height = idx->height;
    	int width = idx->width;
    	
    	int *ym1=NULL,*xm1=NULL,*yp1=NULL,*xp1=NULL;
	if(( ym1 = ( int*) malloc( sizeof( int )*(length + 20))) == NULL )	exit( 1 );
	if(( xm1 = ( int*) malloc( sizeof( int )*(length + 20))) == NULL )	exit( 1 );
	if(( yp1 = ( int*) malloc( sizeof( int )*(length + 20))) == NULL )	exit( 1 );
	if(( xp1 = ( int*) malloc( sizeof( int )*(length + 20 ))) == NULL )	exit( 1 );
__ERROR__		
	for( int i=0; i<points->num; i++ )
	{
		ym1[ i ] = y[ i ]-1; if( ym1[ i ] < 1)	ym1[ i ] = 1;
		xm1[ i ] = x[ i ]-1; if( xm1[ i ] < 1)	xm1[ i ] = 1;
		yp1[ i ] = y[ i ]+1; if( yp1[ i ] > width) yp1[ i ]  = width;
		xp1[ i ] = x[ i ]+1; if( xp1[ i ] > height) xp1[ i ] = height;
	}
	
	//CvSize sizeImage = cvGetSize(idx);	
	//struct POINTS2d idxPoints  = points;
	struct POINTS2d *idup = cvCreatePoints( yp1, x,   length);    
	struct POINTS2d *iddn = cvCreatePoints( ym1, x,   length);
	struct POINTS2d *idlt = cvCreatePoints( y,   xm1, length);
	struct POINTS2d *idrt = cvCreatePoints( y,   xp1, length);
	struct POINTS2d *idul = cvCreatePoints( yp1, xm1, length);
	struct POINTS2d *idur = cvCreatePoints( yp1, xp1, length);
	struct POINTS2d *iddl = cvCreatePoints( ym1, xm1, length);
	struct POINTS2d *iddr = cvCreatePoints( ym1, xp1, length);
__ERROR__	
	// phi_x  = -phi(idlt)+phi(idrt);
	
	CvMat* idxMat = A[0];
	getMatrix( phi, points, idxMat );
		
	CvMat* idupMat = A[1];
	getMatrix( phi, idup, idupMat );
__ERROR__	
	CvMat* iddnMat = A[2];
	getMatrix( phi, iddn, iddnMat );

	CvMat* idltMat = A[3];
	getMatrix( phi, idlt, idltMat );
		
	CvMat* idrtMat = A[4];
	getMatrix( phi, idrt, idrtMat );
__ERROR__		
	CvMat* idulMat = A[5];
	getMatrix( phi, idul, idulMat );
		
	CvMat* idurMat = A[6];
	getMatrix( phi, idur, idurMat );
	
	CvMat* iddlMat = A[7];
	getMatrix( phi, iddl, iddlMat );
	
	CvMat* iddrMat = A[8];
	getMatrix( phi, iddr, iddrMat );
__ERROR__	
	//phi_x  = -phi(idlt)+phi(idrt);
	CvMat* phi_x = A[9];
	cvSub( idrtMat, idltMat, phi_x );
	
	//phi_y  = -phi(iddn)+phi(idup);
	CvMat* phi_y = A[10];
	cvSub( idupMat, iddnMat, phi_y );	
	
	// phi_xx = phi(idlt)-2*phi(idx)+phi(idrt);
	CvMat* phi_xx = A[11];
	cvAdd( idltMat, idrtMat, phi_xx );	
	cvSub( phi_xx, idxMat, phi_xx );
	cvSub( phi_xx, idxMat, phi_xx );
	
	//phi_yy = phi(iddn)-2*phi(idx)+phi(idup);
__ERROR__	
	CvMat* phi_yy = idltMat;
	cvAdd( idupMat, iddnMat, phi_yy );	
	cvSub( phi_yy, idxMat, phi_yy );
	cvSub( phi_yy, idxMat, phi_yy );
__ERROR__	
	//phi_xy =  -0.25*phi(iddl)-0.25*phi(idur)+0.25*phi(iddr)+0.25*phi(idul);
	
	CvMat* phi_xy = A[12];
    	CvMat* f3 = A[13];	
__ERROR__    	
    	cvAddWeighted( iddlMat, 0.25, idurMat, 0.25, 0, f3);
    	
    	cvAddWeighted( iddrMat, 0.25, idulMat, 0.25, 0, phi_xy);
    	cvSub(phi_xy,f3,phi_xy);
    	
    	CvMat *phi_x2 = A[14];
    	CvMat *phi_y2 = A[15];
    	
    	//phi_x2 = phi_x.^2;
    	cvPow( phi_x, phi_x2, 2 );
	//phi_y2 = phi_y.^2;
	cvPow( phi_y, phi_y2, 2 );
__ERROR__	
	
	//curvature = ((phi_x2.*phi_yy + phi_y2.*phi_xx - 2*phi_x.*phi_y.*phi_xy)./...
             //(phi_x2 + phi_y2 +eps).^(3/2)).*(phi_x2 + phi_y2).^(1/2);        
  
	CvMat *arg1 = A[16];
	CvMat *arg2 = A[17];
	CvMat *arg3 = A[18];
	
   	
    	
	CvMat* f1 = A[20];
	cvMul( phi_x2, phi_yy, f1 ); //f1 = phi_x2.*phi_yy
	CvMat* f2 = A[19];
	
	cvMul( phi_y2, phi_xx, f2 ); //f2 = phi_y2.*phi_xx 
	
	
	
	cvMul( phi_x, phi_y, phi_x );//phi_x = phi_x*phi_y
	cvMul( phi_x, phi_xy, f3 ); // f3 = phi_x * phi_xy
	
	//CvMat* f3 = A[21];
	
	cvScale( f3, f3, 2 );	    // f3 = 2*f3;	
	cvAdd( f1, f2, arg1 );      // arg1 = f1 + f2;	
	cvSub( arg1, f3, arg1 );   // arg1 = arg1 - f3	

	cvAddWeighted( phi_x2, 1, phi_y2, 1, eps, arg2 ); // arg2 = phi_x2 + phi_y2 +eps

	cvPow( arg2, arg2, ( float ) 3/2 ); // arg2 = arg2.^(3/2)
	
	cvAdd( phi_x2, phi_y2, arg3 );      //arg3 = phi_x2 + phi_y2
	cvPow( arg3, arg3, ( float )1/2 );  // arg3 = arg3.^(1/2)	
	
	cvDiv( arg1, arg2, arg1 );
	cvMul( arg1, arg3, arg1 );
	float * answer = arg1->data.fl;
	

	cvGetImage( arg1, points, dst );
	
	for (int as = 0; as<22 ;as++)
	{
		cvReleaseMat( &A[as] );
	}
//	free(x);
//	free(y);
	free( points );
//	free( xp1 );	free( yp1 );	free( xm1 );	free( ym1 );
//	free( idup ); free( iddn  ); free( idlt ); free( idrt ); free( idul ); free( idur );
//	free( iddl ); free( iddr );
	
}

void cvMask2phi( IplImage* init_a, IplImage* phi)
{
	/*
	Matlab code and process 
	function phi = mask2phi(init_a)
	phi=bwdist(init_a)-bwdist(1-init_a)+im2double(init_a)-.5;	
	phi=dstB - dstA + dstC ;	
	*/
		
	/*bwdist == void cvDistTransform(const CvArr* src, CvArr* dst, int distance_type=CV_DIST_L2, int mask_size=3, const float* mask=NULL, CvArr* labels=NULL)*/
	
	CvSize sizeImage = cvGetSize( init_a );
	IplImage* dstA = cvCreateImage( sizeImage, IPL_DEPTH_32F, 1 );
	cvZero(dstA);
	cvDistTransform( init_a, dstA, CV_DIST_L1,CV_DIST_MASK_PRECISE);

	IplImage* init_aNOT = cvCreateImage( sizeImage, IPL_DEPTH_8U, 1 );
	cvZero( init_aNOT );
	IplImage* dstB = cvCreateImage( sizeImage, IPL_DEPTH_32F, 1 );	
	cvZero( dstB );
	cvSubRS( init_a, cvScalar(1), init_aNOT, 0 );
	
	//cvNot( init_a, imageNot );	
	cvDistTransform( init_aNOT, dstB, CV_DIST_L1, CV_DIST_MASK_PRECISE );
	cvReleaseImage( &init_aNOT );
	
	IplImage* dstC = cvCreateImage( sizeImage, IPL_DEPTH_32F, 1 );
	cvZero( dstC );
	cvConvertScale( init_a, dstC );
	CvScalar S = cvScalar( 0.5 );

	cvSubS( dstC, S, dstC );
	cvSub( dstB, dstA, phi );
	cvAdd( phi, dstC, phi );
	cvReleaseImage( &dstA );
	cvReleaseImage( &dstB );
	cvReleaseImage( &dstC );
}



/*
IplImage *cvIm2GrayDouble(IplImage* Iimage)
{

	IplImage *IGray32= cvCreateImage(cvGetSize(Iimage),IPL_DEPTH_32F,1);
	cvConvertScale(Iimage,IGray32,1.0,0);
	cvReleaseImage( &Iimage );
	return IGray32;	 
}*/

void cvSussman( IplImage *D, IplImage *dst, float dt )
{
	CvSize sizeImage = cvGetSize(D);
	int depth = D->depth;
	int nChannels = D->nChannels;
	
	IplImage *f2 = cvCreateImage(sizeImage, depth, 1 );	
	cvZero( f2 );
	IplImage *f1 = cvCreateImage(sizeImage, depth, 1 );	
	cvZero( f1 );	
		
	IplImage *a = cvCreateImage(sizeImage, depth, 1 );
	cvZero( a );
	cvShift( D, a, 3 );
	cvSub( D, a, a );
	
	IplImage *b = cvCreateImage(sizeImage, depth, 1 );
	cvZero( b );
	cvShift( D, b, 4 );
	cvSub( b, D, b );
	
	IplImage *c = cvCreateImage(sizeImage, depth, 1 );
	cvZero( c );
	cvShift( D, c, 2 );
	cvSub( D, c, c );
	
	IplImage *d = cvCreateImage(sizeImage, depth, 1 );
	cvZero( d );
	cvShift( D, d, 1 );
	cvSub( d, D, d );
	
	
	IplImage* Index = cvCreateImage(sizeImage, IPL_DEPTH_8U, 1 );
	cvZero( Index );
	cvZero( f1 );

	IplImage *a_n = cvCloneImage( a );
	cvCmpS( a, 0 , Index, CV_CMP_GT );
	//cvZero( f1 );
	cvAndDiff( a_n, cvScalar(0), Index );

		
	IplImage *a_p = cvCloneImage( a );
	cvCmpS( a, 0 , Index, CV_CMP_LT );
	cvAndDiff( a_p, cvScalar(0), Index );
	
	IplImage *b_n = cvCloneImage( b );
	cvCmpS( b, 0 , Index, CV_CMP_GT );
	cvAndDiff( b_n, cvScalar(0), Index );
	
	IplImage *b_p = cvCloneImage( b );
	cvCmpS( b, 0 , Index, CV_CMP_LT );
	cvAndDiff( b_p, cvScalar(0), Index );

	IplImage *c_n = cvCloneImage( c );
	cvCmpS( c, 0 , Index, CV_CMP_GT );
	cvAndDiff( c_n, cvScalar(0), Index );
	
	IplImage *c_p = cvCloneImage( c );
	cvCmpS( c, 0 , Index, CV_CMP_LT );
	cvAndDiff( c_p, cvScalar(0), Index );
	
	IplImage *d_n = cvCloneImage( d );
	cvCmpS( d, 0 , Index, CV_CMP_GT );
	cvAndDiff( d_n, cvScalar(0), Index );
	
	IplImage *d_p = cvCloneImage( d );
	cvCmpS( d, 0 , Index, CV_CMP_LT );
	cvAndDiff( d_p, cvScalar(0), Index );	
	
	IplImage *IndexPos = Index; //reusing with new name
	cvZero( IndexPos );
	IplImage *IndexNeg = cvCreateImage( sizeImage, IPL_DEPTH_8U, 1 );
	cvZero( IndexNeg );
	
	IplImage *Dd = cvCreateImage( sizeImage, depth, 1 );
	cvZero( Dd );

	cvCmpS( D, 0 , IndexPos, CV_CMP_GT );
	cvCmpS( D, 0 , IndexNeg, CV_CMP_LT );


	IplImage* img = IndexPos;
	
	cvMyOverLapping( a_p, b_n, c_p, d_n, IndexPos, Dd );
	cvMyOverLapping( a_n, b_p, c_n, d_p, IndexNeg, Dd );
	
	//  D = D - dt .* sussman_sign(D) .* dD;
	
	cvZero( f1 );
	
	cvSussmanSign( D, f1 );
	cvScale( f1, f1, dt );
	cvMul( f1, Dd, f1 );
	
	cvSub( D, f1, dst );
	
	cvReleaseImage( &Dd ); //error
	cvReleaseImage( &f1 );
	cvReleaseImage( &f2 );
	cvReleaseImage( &IndexPos );
	cvReleaseImage( &IndexNeg );
	cvReleaseImage( &a );
	cvReleaseImage( &a_n );
	cvReleaseImage( &a_p );
	cvReleaseImage( &b );
	cvReleaseImage( &b_n );
	cvReleaseImage( &b_p );
	cvReleaseImage( &c );
	cvReleaseImage( &c_n );
	cvReleaseImage( &c_p );	
	cvReleaseImage( &d );
	cvReleaseImage( &d_n );
	cvReleaseImage( &d_p );
}

//function S = sussman_sign(D)
 // S = D ./ sqrt(D.^2 + 1);
 
void cvSussmanSign( IplImage *D, IplImage *S ) 
{	
	cvMul( D, D, S );
	cvAddS( S, cvScalar(1), S );
	cvmySqrt( S ); 
	cvDiv( D, S, S, 1);		
}

void writeImageToFile(IplImage *src, char *fname, int option)
{
	FILE *fp =fopen(fname,"w+");
	
	if(fp == NULL)
	{	
		fprintf(stderr,"file error");
		return;
	}
	int widthStep, baseIndex, width, height;
	float * srcData1; 
	uchar* srcData;
	switch(option)
	{
	case 0:
	widthStep = src->widthStep;
	baseIndex  =-1*widthStep;
	width = src->width;
	height = src->height;
	srcData = (uchar *)src->imageData; 
		
	for(int i =0;i<height;i++)
	{	
		baseIndex  += widthStep;
		
		for( int j = 0 ; j<width; j++ )
		{		
			fprintf(fp,"%u ", srcData[baseIndex+j]);
			
		}	
		fprintf(fp,"\n");
		
	}
	break;
	case 1:
	widthStep = (src->widthStep)/4;
	baseIndex  =-1*widthStep;
	width = src->width;
	height = src->height;
	srcData1 = (float *)src->imageData; 
		
	for(int i =0;i<height;i++)
	{	
		baseIndex  += widthStep;
		
		for( int j = 0 ; j<width; j++ )
		{		
			fprintf(fp,"%f ", srcData1[baseIndex+j]);
			
		}	
		fprintf(fp,"\n");
		
	}
	break;
	}
	fclose(fp);
}

