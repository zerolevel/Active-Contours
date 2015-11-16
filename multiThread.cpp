#include "multiThread.h"
#include "imageOperations.h"

void* multiThread1_1(void* Arg1)
{
	void** Arg = (void**) Arg1;
	cvZero( (IplImage*)Arg[1] );		
	cvZero( (IplImage*)Arg[4] );
	cvAndDiff( (IplImage*)Arg[4], cvScalar(*(float*)Arg[6]), (IplImage*)Arg[5] ); 
	cvSub( (IplImage*)Arg[0], (IplImage*)Arg[4], (IplImage*)Arg[1] ); //f1=f-f4
	cvPow( (IplImage*)Arg[1], (IplImage*)Arg[1], 2.0 );
	return NULL;
}

void* multiThread1_2(void* Arg1)
{
	void** Arg = (void**) Arg1;
	//fprintf(stderr, "MT %d\n", __LINE__);
	cvZero( (IplImage*)Arg[2] );
	cvZero( (IplImage*)Arg[3] );
	//fprintf(stderr, "MT %d\n", __LINE__);
	cvAndDiff( (IplImage*)Arg[3], cvScalar(*(float*)Arg[7]), (IplImage*)Arg[5] ); 
	cvSub( (IplImage*)Arg[0], (IplImage*)Arg[3], Arg[2] );  //f2=f-f3
	cvPow( (IplImage*)Arg[2], (IplImage*)Arg[2], 2.0 );
	return NULL;
}
void multiThread1(void** Arg)
{
	//fprintf( stderr, "MT %d\n", __LINE__ );
	pthread_t th1, th2;
	pthread_create( &th1, NULL, &multiThread1_1, (void*)Arg );
	pthread_create( &th2, NULL, &multiThread1_2, (void*)Arg );
	pthread_join( th1, NULL );	
	pthread_join( th2, NULL );
	//fprintf(stderr, "MT %d\n", __LINE__);
	cvSub( (IplImage*)Arg[1], (IplImage*)Arg[2], (IplImage*)Arg[0] );
	return;
}

