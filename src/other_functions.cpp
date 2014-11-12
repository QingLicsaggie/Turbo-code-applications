/*---------------------------------------------------------------
* Copyright NOV 2010 CSE@TAMU
* All rights reserved.
*
* other_functions.cpp
*
* This script implements the functions used in turbo coding and decoding.
*
* Version: 1.0
* Programmed By Qing(Tsing) Li
* Last updated date: Nov 2010
---------------------------------------------------------------*/
#include "turbo_code_Log_MAP.h"
#include "other_functions.h"
#include "huffmancode.h"

extern string text;                                 /*original file*/
extern string codes;                                /*encoded file*/
extern vector<node> code;	                                         
extern vector<node> result;  
/*---------------------------------------------------------------
FUNCTION: 
	void gen_source(int *data, int length)
	
DESCRIPTION:
	This function generate the source bits for simulation.

PARAMETERS:
	INPUT:
		length - Length of needed data.
	OUTPUT:
		data - Contains pointer to source data sequence.

RETURN VALUE:
	None
---------------------------------------------------------------*/
void gen_source(int *data, int length)
{
	int i;

	for (i=0; i<length; i++)
		*(data+i) = (int)codes[i] - 48;
}
/*---------------------------------------------------------------
FUNCTION: 
	void AWGN(int *send, float *r, float sigma, int totallength)
	
DESCRIPTION:
	This function simulate a AWGN channel.

PARAMETERS:
	INPUT:
		send - Input bit sequence need to add noise.
		sigma - Standard deviation of AWGN noise
		totallength - Length of "send".
	OUTPUT:
		r - Contains pointer to the data sequence added with gaussian white noise.

RETURN VALUE:
	None
---------------------------------------------------------------*/
void AWGN(float *send, float *r, float sigma, int totallength)
{
	int i;

	double *noise = (double *)malloc(sizeof(double)*totallength);
	 
	srand((int)time(0));
	
	double seed =  3.0 - (double)random(100)/100;
	
	mgrns(0,sigma,seed,totallength,noise);
	
	for(i=0; i<totallength; i++)
	{
		*(r+i) = (float)( *(send+i) + *(noise+i) );
	}
	ofstream output( "code.txt", ios::out );
	output<<r;
	output.close();
	free(noise);
}

/*-----------------------------------------
* FUNCTION: mgrns(double mean,double sigma,double seed,int n,double *a)
* DESCRIPTION:
*	INPUT   mean
*           sigma 
*			seed: a random seat
*  OUTPUT
*           a£º Gaussian sequence of length n
* 
*  RETURN VALUE
*           NONE
*------------------------------------------*/
void mgrns(double mean,double sigma,double seed,int n,double *a)
{  
	int i,k,m;
    double s,w,v,t;
    s=65536.0; w=2053.0; v=13849.0;
    for (k=0; k<=n-1; k++)
	{
		t=0.0;
		for (i=1; i<=12; i++)
        { 
			seed=seed*w+v; m=(int)(seed/s);
            seed=seed-m*s; t=t+(seed)/s;
        }/*According to the theory we learn in mathematical course*/
        *(a+k)=mean+sigma*(t-6.0);
    }
    return;
}
