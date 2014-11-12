
/*---------------------------------------------------------------
* Copyright NOV 2010 CSE@TAMU
* All rights reserved.
*
* main.cpp
*
* This script implements the functions used in turbo coding and decoding.
*
* Version: 1.0
* Programmed By Qing(Tsing) Li
* Last updated date: Nov 2010
---------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include<iostream>
#include <math.h>
#include <time.h>
#include<iomanip>
#include<cstdio>
#include<memory>
#include <fstream>
#include<cstring>
#include<queue>
#include<string> 
#include<vector>

#include "turbo_code_Log_MAP.h"
#include "other_functions.h"
#include "huffmancode.h"

int  FRAME_LENGTH;
extern float rate_coding;
extern  TURBO_G turbo_g;
extern int *mx_puncture_turbo;
extern string text;                                 /*original file*/
extern string codes;                                /*encoded file*/
//extern vector<node> code;	                                         
extern vector<node> result;  
extern int M_num_reg;
extern float *LLR_all_turbo;
extern int *reverse_index_randomintlvr;
void main()
{

	 FILE *stream;	/* point to "simu_report.txt" which is used to record the results */
	 if ((stream=fopen( "simu_report.txt", "w" )) == NULL)
    {
	  printf("\nError! Can not open file simu_report.txt\n");
	  exit(1);
    } 

	 FILE *fp;	/* point to "wrong_LLR.txt" which is used to record the LLR results */
	 if ((fp=fopen( "wrong_LLR.txt", "w++" )) == NULL)
    {
	  printf("\nError! Can not open file wrong_LLR.txt\n");
	  exit(1);
    } 

	 FILE *filepointer;	/* point to "right_LLR.txt" which is used to record the LLR results */
	 if ((filepointer=fopen( "right_LLR.txt", "w++" )) == NULL)
    {
	  printf("\nError! Can not open file right_LLR.txt\n");
	  exit(1);
    } 

	 fprintf(stream, "\n================================================================= \n");
	 fprintf(stream, "\n1. For a given text file, use Huffman coding algorithm to compress it \n");
	 fprintf(stream, "\n2. Adding Error Correcting Code, here we use Turbo code. \n");
	 fprintf(stream, "\n3. Simulate the real situtation, meaning adding AWGN noise \n");
	 fprintf(stream, "\n4. Using Turbo decoding skill to restore the impaired file \n");
	 fprintf(stream, "\n5. Using Huffman decoding skill to decode the file\n");
	 fprintf(stream, "\n=========================================================== \n ");


	 fprintf(stream,"\n                                                                \n");
	 fprintf(stream, "\n================================================================= \n");
	 fprintf(stream, "\n  HUFFMAN CODE PART \n");
	 fprintf(stream, "\n =================================================================\n");

	 
	 encode();

	 FRAME_LENGTH = codes.length();

	 fprintf(stream,"\n                                                                \n");
	 fprintf(stream, "\n================================================================= \n");
	 fprintf(stream, "\n  TURBO CODE PART \n");
	 fprintf(stream, "\n =================================================================\n");
	
/*====================================================================================================================*/	 

	int			*trafficflow_source = NULL,traffic_source_length, *trafficflow_decoded = NULL, err_bit_num_traffic[8],i, j, ien,total;
	float		*coded_trafficflow_source = NULL,*trafficflow_for_decode = NULL,EbN0dB = 0,en, sgma,err_bit_rate_traffic[8], max_LLR = 0.0;


	//double Eb_N0_dB[8] = {0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.05};
    // double Eb_N0_dB[8] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4};
   //double Eb_N0_dB[8] = {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75};
	//double Eb_N0_dB[8] = {0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1};
    //double Eb_N0_dB[8] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8};
	double Eb_N0_dB[8] = {2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8};
	
	clock_t start, end;


	trafficflow_source=(int *)malloc(FRAME_LENGTH*sizeof(int));
	
	coded_trafficflow_source=(float *)malloc(2*FRAME_LENGTH*sizeof(float));

	trafficflow_for_decode=(float *)malloc(2*FRAME_LENGTH*sizeof(float));
	
	trafficflow_decoded=(int *)malloc(FRAME_LENGTH*sizeof(int));

	traffic_source_length = FRAME_LENGTH;
	
	/*-----------------------------------------------------------------*/	
	/*-----------------------------------------------------------------*/
      TurboCodingInit(traffic_source_length);
	/*-----------------------------------------------------------------*/
	/*-----------------------------------------------------------------*/	

		/*====   output the simulation parameters to screen======*/
		printf("\n======================== Turbo code simulation :========================\n");
		printf("\n Some parameters are as follows: \n");
		printf("\nlog map decoder\n");
		printf("frame length : %d \n", traffic_source_length);
		printf("generators g = \n");

		for (i=0; i<turbo_g.N_num_row; i++)
		{
			for (j=0; j<turbo_g.K_num_col; j++)
			{
				printf(" %d  ", *(turbo_g.g_matrix+i*turbo_g.K_num_col+j));
			}
			printf("\n");
		}


		/*	to "simu_report.txt"*/
	    fprintf(stream,"\n======================== Turbo code simulation ========================\n");
		fprintf(stream, "\n Some parameters are as follows: \n");
		fprintf(stream, "/*------------------------------*/\n");

		fprintf(stream, "log map decoder\n");
		fprintf(stream, "frame length : %d \n", traffic_source_length);
		fprintf(stream, "generators g = \n");
		for (i=0; i<turbo_g.N_num_row; i++)
		{
			for (j=0; j<turbo_g.K_num_col; j++)
			{
				fprintf(stream, "%d ", *(turbo_g.g_matrix+i*turbo_g.K_num_col+j));
			}
			fprintf(stream, "\n");
		}

/*=============================================================================================================*/
	for (ien=0; ien<8; ien++)
	{
		total = 0;
        	                                                   /*printf out SNR to screen or text file*/
		printf("\n /*---------------------------------------------*/\n");
		printf("\n Simulation %d \n", ien + 1);
        
		EbN0dB = (float)Eb_N0_dB[ien];
		
		en = (float)pow(10,(EbN0dB)/10);

		sgma = (float)(1.0/sqrt(2*rate_coding*en));
		
		err_bit_num_traffic[ien] = 0;
		
		err_bit_rate_traffic[ien] = 0.0;
		
	
		printf("\n Eb/N0 = %f \n", EbN0dB);

		if (TURBO_PUNCTURE)	
     		printf(" punctured to rate 1/2\n");
		else
		    printf(" unpunctured\n");
	

		printf(" iteration number : %d \n", N_ITERATION);


		fprintf(stream, "\n Eb/N0 = %f \n", EbN0dB);

		if (TURBO_PUNCTURE)
			fprintf(stream, "\n punctured to rate 1/2\n");
		else
			fprintf(stream, "\n unpunctured\n");


		fprintf(stream, "iteration number : %d \n", N_ITERATION);
		/*=================================================================
                     Probess begins!
		===================================================================*/
			start=clock();/*staring time*/

			traffic_source_length = FRAME_LENGTH;

			gen_source(trafficflow_source, traffic_source_length);
		/*-----------------------------------------------------------------
			 Encoding process
		-----------------------------------------------------------------*/
			TurboCodingTraffic(trafficflow_source, coded_trafficflow_source, &traffic_source_length);
		/*-----------------------------------------------------------------
            With noise
		-----------------------------------------------------------------*/ 
			AWGN(coded_trafficflow_source, trafficflow_for_decode, sgma,traffic_source_length);
		/*-----------------------------------------------------------------
           Decoding process
		-----------------------------------------------------------------*/	
			TurboDecodingTraffic(trafficflow_for_decode, trafficflow_decoded, &traffic_source_length, EbN0dB);
		/*-----------------------------------------------------------------
                  Analysis process
		-----------------------------------------------------------------------*/
			printf("\n Analysis:\n");
			fprintf(stream,"\n /*---------------------------------------------*/\n");
			fprintf(stream,"\n Analysis:\n");

		    printf("\n The whole process takes ");
			fprintf(stream,"\n The whole process takes ");
			end=clock();
		    fprintf(stream,"%5f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
		    printf( "%5f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
			fprintf(fp,"The corresponding LLR of error bits for simulation %d :\n\n", ien + 1);
			fprintf(filepointer,"The corresponding LLR of error bits for simulation %d :\n\n", ien + 1);

			for (i=0; i< traffic_source_length ; i++)
			{
				if (*(trafficflow_source+i) != *(trafficflow_decoded+i))
				{
					err_bit_num_traffic[ien] = err_bit_num_traffic[ien]+1;
					fprintf(fp, "LLR[%d] =  %f\n\n",*(reverse_index_randomintlvr + i), *(LLR_all_turbo + *(reverse_index_randomintlvr + i) ));
					if(fabs(*(LLR_all_turbo + *(reverse_index_randomintlvr + i) )) > max_LLR)
						max_LLR = fabs(*(LLR_all_turbo + *(reverse_index_randomintlvr + i)));
				}
				else
				{
					fprintf(filepointer, "LLR[%d] =  %f\n\n",*(reverse_index_randomintlvr + i), *(LLR_all_turbo + *(reverse_index_randomintlvr + i) ));
              
				}
			}
			
		
			for (i=0; i< traffic_source_length ; i++)
			{


				if (*(trafficflow_source+i) == *(trafficflow_decoded+i)&&fabs(*(LLR_all_turbo +*(reverse_index_randomintlvr + i))) < max_LLR)
				{
					total ++;
				}
			}

			float wrong = (float) total/(float) traffic_source_length;
			err_bit_rate_traffic[ien] = (float)err_bit_num_traffic[ien]/(traffic_source_length) ;
			printf("\n The Error bit probability is %f  \n",err_bit_rate_traffic[ien]);
			fprintf(stream,"\n The Error bit probability is %f \n",err_bit_rate_traffic[ien]);
			fprintf(fp,"maxmum LLR of error bit is %f and SNR is %f \n\n", max_LLR,EbN0dB);
			fprintf(fp,"Among correct bits, there are %f % of them whose LLR less than max_LLR\n", wrong);
			fprintf(fp,"===========================================================================\n\n");
			fprintf(filepointer,"===========================================================================\n\n");
	/*=================================================================
                  HUFFMAN DECODE
	===================================================================*/
			decode( ); 
		}
			fprintf(fp,"maxmum LLR of error bit is %f\n\n", max_LLR);
	
	/*-----------------------------------------------------------------*/	
	/*-----------------------------------------------------------------*/	
		 //TurboCodingRelease();
	/*-----------------------------------------------------------------*/
	/*-----------------------------------------------------------------*/	

	fclose(stream);
	fclose(fp);
	fclose(filepointer);
	free(LLR_all_turbo);
	/*free(trafficflow_source);
	free(coded_trafficflow_source);
	free(trafficflow_for_decode);
	free(trafficflow_decoded);*/

}
