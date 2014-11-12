/*---------------------------------------------------------------
* Copyright NOV 2010 CSE@TAMU
* All rights reserved.
*
* turbo_code_Log_MAP.h
*
* This script implements the functions used in turbo coding and decoding.
*
* Version: 1.0
* Programmed By Qing(Tsing) Li
* Last updated date: Nov 2010
---------------------------------------------------------------*/

#ifndef	TRUBO_CODE_LOG_MAP_H
#define	TRUBO_CODE_LOG_MAP_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#define random(x) (rand()%x)                /*we use this to generate a random number*/
/*==================================================*/
/*==================================================*/
/*	parameters for simulation	*/
#define N_ITERATION			18	/* number of iteration times */
#define TURBO_PUNCTURE		1				/* puncture or not
												1 for punctured
												0 for unpunctured */
/*==================================================*/
/* parameters for the code generator */
#define COLUMN_OF_G		4	
#define G_ROW_1			13
#define G_ROW_2			15

/*==================================================*/
/* structures */

/* code generator */
typedef struct
{
	int N_num_row;
	int K_num_col;
	int *g_matrix;
} TURBO_G;

/* trellis */
typedef struct
{
	int *mx_nextout;
	int *mx_nextstat;
	int *mx_lastout;
	int *mx_laststat;
	
} TURBO_TRELLIS;

/*==================================================*/
/* interfaces to outside */

void TurboCodingTraffic(int *trafficflow_source, float *coded_trafficflow_source,
						int *traffic_source_length);

void TurboDecodingTraffic(float *trafficflow_for_decode, int *trafficflow_decoded,
						  int *trafficflow_length, float EbN0dB);


void TurboCodingInit(int length_total);


void TurboCodingRelease();

/*==================================================*/
/*==================================================*/
/* functions used inside */
int gen_g_matrix(int k_column, int g_row1, int g_row2, int *mx_g_turbo);

void gen_trellis(int length_total);

void int2bin(int intstat, int *tempstat, int length);

int bin2int(int *binseq, int length);

float random_turbo();

void gen_rand_index(int length);

void randominterleaver_int(int *data_unintlvr, int *interleaverddata, int length, int *index);

void randominterleaver_float(float *data_unintlvr, float *interleaverddata, int length, int *index);

void random_deinterlvr_int(int *data_unintlvr, int *interleaverddata, int length, int *index);

void random_deinterlvr_float(float *data_unintlvr, float *interleaverddata, int length, int *index);

void interleave_int(int *data_unintlvr, int *interleaveddata, int typeofinterleave, int nturbo, int* index);

void interleave_float(float *data_unintlvr, float *interleaveddata, int typeofinterleave, int nturbo, int* index);

void de_interleave_int(int *data_unintlvr, int *interleaveddata, int typeofinterleave, int nturbo, int *index);

void de_interleave_float(float *data_unintlvr, float *interleaveddata, int typeofinterleave, int nturbo, int *index);

void encoderm_turbo(int *source, int *send_turbo, int len_info);

void rsc_encode(int *source, int *rsc, int len_info, int *state);

int encode_bit(int inbit, int *stat);

void gen_mx_punc();

void puncture(int *data_unpunc, int length_unpunc, int *data_punctured);

void depuncture(float *receive_punc, int length_punc, float *receive_depunced);

void demultiplex(float *rec_turbo, int len_total, float *yk_turbo);

void Log_MAP_decoder(float *recs_turbo, float *La_turbo, float *LLR_all_turbo, int len_total);

float E_algorithm(float x, float y);

float E_algorithm_seq(float *data_seq, int length);

void decision(float *LLR_seq, int length, int *output);

/*==================================================*/

#endif