// rkf45.h

#ifndef _RKF45
#define _RKF45

void RKF45(void (*F)(float t1, float* y1, float* dy1),
	int NEQN, float* Y, float* T, float* TOUT,
	float* RELERR, float* ABSERR,
	int* IFLAG, float* WORK, int* IWORK);

#endif