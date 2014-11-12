
/*---------------------------------------------------------------
* Copyright NOV 2010 CSE@TAMU
* All rights reserved.
*
* other_functions.h
*
* This script implements the functions used in turbo coding and decoding.
*
* Version: 1.0
* Programmed By Qing(Tsing) Li
* Last updated date: Nov 2010
---------------------------------------------------------------*/#ifndef	OTHER_FUNCTIONS_H
#define	OTHER_FUNCTIONS_H

#include <stdlib.h>
#include "turbo_code_Log_MAP.h"

void gen_source(int *data, int length);

void AWGN(float *send, float *r, float sigma, int totallength);

void mgrns(double mean, double sigma, double seed, int n, double *a);

#endif