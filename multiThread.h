#ifndef __multiThread__
#define __multiThread__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cv.h>
#include <highgui.h>
#include <pthread.h>

void* multiThread1_1(void* Arg1);
void* multiThread1_2(void* Arg1);
void multiThread1(void** Im);


#endif
