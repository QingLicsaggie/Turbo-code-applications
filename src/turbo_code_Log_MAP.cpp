/*---------------------------------------------------------------
* Copyright NOV 2010 CSE@TAMU
* All rights reserved.
*
* turbo_code_Log_MAP.cpp
*
* This script implements the functions used in turbo coding and decoding.
*
* Version: 1.0
* Programmed By Qing(Tsing) Li
* Last updated date: Nov 2010
---------------------------------------------------------------*/
#include<iostream>
#include<iomanip>
#include<cstdio>
#include<memory>
#include <fstream>
#include<cstring>
#include<queue>
#include<string> 
#include<vector>
#include "turbo_code_Log_MAP.h"

using namespace std;
using std::fixed;
using std::setprecision;

/*==================================================*/
/* lookup table used in Log-MAP decoder */
const double lookup_index_Log_MAP[16] = {0.0, 0.08824, 0.19587, 0.31026, 0.43275, 0.56508,
								0.70963, 0.86972, 1.0502, 1.2587, 1.5078, 1.8212,
								2.2522, 2.9706, 3.6764, 4.3758};
const double lookup_table_Log_MAP[16] = {0.69315, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35,
								0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.025, 0.0125};
/*==================================================*/
/* parameters on g matrix */
int M_num_reg = COLUMN_OF_G - 1;		/* number of registers */
const int n_states = 8;						/* number of states : pow(2, M) */

/*==================================================*/
/* A number used in Log-MAP decoder */
#define INFTY 1E20
#define SATELITE 5555555555.55555555555555555555
#define THRESHOLD 75.000
/*==================================================*/

/* global memory */

int *index_randomintlvr;		/* index of random intleaver */
int *reverse_index_randomintlvr;  /* reverse of random intleaver*/
/* puncture matrix */
int *mx_puncture_turbo;

TURBO_G turbo_g;				/* code generator struct */

TURBO_TRELLIS turbo_trellis;	/* turbo_trellis struct */	

float rate_coding;				/* the rate of coding */
float *LLR_all_turbo;
int *up_state;                  /*memory state for the up rsc encoder*/
int *down_state;                 /*memory state for the down rsc encoder*/
vector<int> locations_of_erasure;

extern string decoded_codes;    
/*=====================================================*/


/*---------------------------------------------------------------
FUNCTION: 
	TurboCodingInit(int *traffic_source_length)
	
DESCRIPTION:
	Initialize the parameters of coding.

PARAMETERS:
	INPUT:
		traffic_source_length - The pointer that give out the length of the source bits with terminating bits.

RETURN VALUE:
	None.
---------------------------------------------------------------*/
void TurboCodingInit(int length_info)
{
	turbo_g.N_num_row = 2;		/* initialize number of rows and columns of g */
	turbo_g.K_num_col = COLUMN_OF_G;

	int length_total = length_info + M_num_reg  + 1;

	/* caculate the coding rate */
	rate_coding = TURBO_PUNCTURE? (float)(0.5) : (float)(0.333333);

	turbo_g.g_matrix=(int *)malloc(turbo_g.N_num_row*turbo_g.K_num_col*sizeof(int));

	/* generate the up_state and down_state*/
	up_state = (int *)malloc(M_num_reg*sizeof(int));

	down_state=(int *)malloc(M_num_reg*sizeof(int));

	/*initialize the up_state and down_state*/
	for( int i = 0 ; i < M_num_reg ; i++)
		*(up_state + i) = *(down_state + i) = 0;

	/* generate the code generator matrix */
	if (!gen_g_matrix(turbo_g.K_num_col, G_ROW_1, G_ROW_2, turbo_g.g_matrix))
	{
		printf("error number of G\n");
		exit(1);
	}

	/* generate puncture matrix */
	if (TURBO_PUNCTURE)
	{
		mx_puncture_turbo=(int *)malloc(sizeof(int)*3*2);

		gen_mx_punc();	/*	generate the matrix for puncture */

	}

	/* initialize the trellis struct */
	turbo_trellis.mx_lastout=(int *)malloc(sizeof(int)*(length_total + 1)*n_states*4);

	turbo_trellis.mx_laststat=(int *)malloc(sizeof(int)*(length_total+1)*n_states*2);
	 
	turbo_trellis.mx_nextout=(int *)malloc(sizeof(int)*(length_total+1)*n_states*4);
	
	turbo_trellis.mx_nextstat=(int *)malloc(sizeof(int)*(length_total+1)*n_states*2);

	gen_trellis(length_total);

	/* malloc memory for the index of random interleaver */
	index_randomintlvr=(int *)malloc(length_total*sizeof(int));
	reverse_index_randomintlvr = (int *)malloc(length_total*sizeof(int));
}
/*---------------------------------------------------------------
FUNCTION: 
	TurboCodingTraffic(int *trafficflow_source, int *coded_trafficflow_source,
						int *traffic_source_length)
	
DESCRIPTION:
	This function encodes the traffic flow bits.

PARAMETERS:
	INPUT:
		trafficflow_source - Contains pointer to the source bits sequence.
		traffic_source_length - The pointer that give out the length of the source bits.
	OUTPUT:
		coded_trafficflow_source - Contains pointer to the coded bits sequence.

RETURN VALUE:
	None.
---------------------------------------------------------------*/
void TurboCodingTraffic(int *trafficflow_source, float *coded_trafficflow_source,
						int *traffic_source_length)
{
	int i;
	int *temp_send = NULL;
	int *send = NULL;
	int *send_punc = NULL;

	int length_info = *traffic_source_length;

	send = (int *)malloc((3*length_info+3*(M_num_reg +1))*sizeof(int));

	send_punc = (int *)malloc(2*(length_info+M_num_reg +1)*sizeof(int));
	
	gen_rand_index(length_info+M_num_reg +1);

	encoderm_turbo(trafficflow_source, send, length_info);/* encode and BPSK modulate */

	temp_send = send;

	if (TURBO_PUNCTURE)			/* punture the coded bits */
	{  
		//printf("\n After puncturing, the bits are as follows:\n");
		puncture(send, 3*length_info+3*(M_num_reg+1), send_punc);
		temp_send = send_punc;
	}

	for (i=0; i<(3-TURBO_PUNCTURE)*(length_info+M_num_reg +1); i++)
	{
		*(coded_trafficflow_source+i) = (float) *(temp_send+i);
	}

	*traffic_source_length = (3-TURBO_PUNCTURE)*(length_info+M_num_reg+1);

	free(send);
	free(send_punc);
}
/*---------------------------------------------------------------
FUNCTION: 
	TurboDecodingTraffic(float *trafficflow_for_decode, int *trafficflow_decoded,
							int *trafficflow_length, float EbN0dB)
DESCRIPTION:
	This function decodes the received traffic flow bits.

PARAMETERS:
	INPUT:
		trafficflow_for_decode - Contains pointer to the received bits sequence.
		trafficflow_length - The pointer that give out the length of the received bits.
		EbN0dB - Eb/N0 in dB.
	OUTPUT:
		trafficflow_decoded - Contains pointer to the decoded bits sequence.

RETURN VALUE:
	None.
---------------------------------------------------------------*/
void TurboDecodingTraffic(float *trafficflow_for_decode, int *trafficflow_decoded,
						  int *trafficflow_length, float EbN0dB)
{
	int i;
	int length_total, length_info;
	int iteration;

	float *receive_punc = NULL;						/*	receiving data	*/
	float *yk_turbo = NULL;							/* include ys & yp */
	float *La_turbo, *Le_turbo;		               /*extrinsic information & LLR*/

	int *tempout;
                         
	float en_rate = (float)pow(10, EbN0dB*0.1);
	float Lc_turbo = 4*en_rate;			/* reliability value of the channel */


	if (TURBO_PUNCTURE)
	{
		length_info = (*trafficflow_length)/2 - (M_num_reg+1) ;
	}
	else
	{
		length_info = (*trafficflow_length)/3 - (M_num_reg+1);
	}
	
	length_total = length_info +  M_num_reg +1 ;

	receive_punc=(float *)malloc(3*length_total*sizeof(float));

	yk_turbo=(float *)malloc(4*length_total*sizeof(float));

	La_turbo=(float *)malloc(length_total*sizeof(float));
	
	Le_turbo=(float *)malloc(length_total*sizeof(float));
	
	LLR_all_turbo=(float *)malloc(length_total*sizeof(float));

	tempout=(int *)malloc(length_total*sizeof(int));

	if (TURBO_PUNCTURE)		/* fill in the punctured bits and demultiplex */
	{   
		depuncture(trafficflow_for_decode,*trafficflow_length, receive_punc);
		demultiplex(receive_punc, length_total, yk_turbo);
	}
	else
	{
		demultiplex(trafficflow_for_decode, length_total , yk_turbo);
	}
	
	/*	scale the data	*/
	for (i=0; i<4*length_total; i++)
	{
		*(yk_turbo+i) = (float)( *(yk_turbo+i) * Lc_turbo );
		//printf("%f \n", *(yk_turbo+i) );
	}
	/*initialize the La, Le and LLR*/
	for (i=0; i<length_total; i++)
	{
		*(La_turbo+i) = *(Le_turbo+i) = *(LLR_all_turbo+i) = 0;
	}
    

	for (iteration=0; iteration<N_ITERATION; iteration++ )		/* start iteration */
	{   

		//printf("\n iteration %d \n ", iteration +1 );
	
		/*============================== decoder one:===================================== */

		/* get extrinsic information from decoder two */
		//printf("\nAfer the deinterleaving\n");

		random_deinterlvr_float(La_turbo,  Le_turbo, length_total,index_randomintlvr);
		

		//printf("\nLOG_MAP_DECODER 1\n");
		Log_MAP_decoder(yk_turbo, La_turbo, LLR_all_turbo, length_total);
		
		
		/* caculate the extrinsic information */
		for (i=0; i<length_total ; i++)
		{
			*(Le_turbo+i) = *(LLR_all_turbo+i) - *(La_turbo+i) - *(yk_turbo + 2*i);
			// printf( "Le_turbo[%d] = %f \n",i,  *(Le_turbo+i));
		}

		/* =====================================decoder two:============================== */
		/* get extrinsic information from decoder one */

		 randominterleaver_float(Le_turbo, La_turbo, length_total,index_randomintlvr);
	
		
		 //printf("\nLOG_MAP_DECODER 2\n");
		 Log_MAP_decoder(yk_turbo+2*(length_total), La_turbo, LLR_all_turbo,length_total);
	

		/* caculate the extrinsic information */
		for(i=0; i<length_total; i++)
		{
			*(Le_turbo+i) = *(LLR_all_turbo+i) - *(La_turbo+i) - *(yk_turbo+2*length_total+2*i);

		}
		/* end of decoder two */

	}
	
	/* get the decision bits from LLR gained from these decoder */
	//printf("\nMaking decisions.\n");

	decision(LLR_all_turbo, length_total, tempout);

	random_deinterlvr_int(trafficflow_decoded,tempout,length_total,index_randomintlvr);

	//printf("\n Finish de_interleaving! \n");
	
	decoded_codes = "";
	for( i = 0 ; i < length_info ; i++ )
		decoded_codes += *(trafficflow_decoded + i) + 48;

	*trafficflow_length = length_info;
	
	free(receive_punc);
	free(yk_turbo);
	free(La_turbo);
	free(Le_turbo);
	//free(LLR_all_turbo);
	free(tempout);
}

/*---------------------------------------------------------------
FUNCTION: 
	TurboCodingFinish()
	
DESCRIPTION:
	Free the memory of turbo coding.

PARAMETERS:
	None.

RETURN VALUE:
	None.
---------------------------------------------------------------*/
void TurboCodingRelease()
{
	free(turbo_g.g_matrix);

	if (TURBO_PUNCTURE)
	{		
		free(mx_puncture_turbo);
	}

	free(turbo_trellis.mx_lastout);
	free(turbo_trellis.mx_laststat);
	free(turbo_trellis.mx_nextout);
	free(turbo_trellis.mx_nextstat);
	free(index_randomintlvr);
	
}

/*---------------------------------------------------------------
FUNCTION: 
	int gen_g_matrix(int k_column, int g_row1, int g_row2, int *mx_g_turbo)
	
DESCRIPTION:
	This function generates the code generator.

PARAMETERS:
	INPUT:
		k - Number of columns of g matrix.
		g_row1 - An octal number indicate the first row of turbo_g.
		g_row2 - An octal number indicate the second row of turbo_g.
	OUTPUT:
		mx_g - Contains pointer to g matrix.

RETURN VALUE:
	1 - If the matrix was successfully generated.
	0 - Error occurred generating the g matrix.
---------------------------------------------------------------*/
int gen_g_matrix(int k_column, int g_row1, int g_row2, int *mx_g_turbo)
{
	int i, position;
	int high_num, low_num;
	
	/* for the first row */
	high_num = g_row1;
	position = 1;
	while (high_num>0)
	{
		low_num = high_num%10;
		if (low_num>7)
		{
			return 0;
		}
		high_num = high_num/10;

		for (i=k_column-(position-1)*3-1; i>=0 && i>=k_column-position*3; i--)
		{
			*(mx_g_turbo+i) = low_num%2;
			low_num = low_num/2;
		}
		position++;
		if (i<0)
		{
			break;
		}
	}

	/*for the second row */
	high_num = g_row2;
	position = 1;
	while (high_num>0)
	{
		low_num = high_num%10;
		if (low_num>7)
		{
			return 0;
		}
		high_num = high_num/10;

		for (i=k_column-(position-1)*3-1; i>=0 && i>=k_column-position*3; i--)
		{
			*(mx_g_turbo+k_column+i) = low_num%2;
			low_num = low_num/2;
		}
		position++;
		if (i<0)
		{
			break;
		}
	}
	return 1;
}
/*---------------------------------------------------------------
FUNCTION: 
	void gen_trellis(int length)
	
DESCRIPTION:
	This function generates the turbo_trellis structure.

PARAMETERS:
	INPUT 
	   length - the length of souce with terminating bits

RETURN VALUE:
	None
---------------------------------------------------------------*/
void gen_trellis(int length)
{
	int i, j, k, l;
	int dk_turbo, ak_turbo, outbit;

	int *tempstat;

	tempstat=(int *)malloc(sizeof(int)*M_num_reg);

	/*Initialize*/
	for( l = 0; l <= length ; l++)

		for( i = 0 ; i < n_states ; i ++)

			for(j = 0; j < 2; j++)

			{	

					*(turbo_trellis.mx_nextout + l*n_states*4 + i*4 + 2*j ) = 0;

						
					*(turbo_trellis.mx_nextout + l*n_states*4 + i*4 + 2*j + 1 ) = 0;


					*(turbo_trellis.mx_nextstat + l*n_states*2 + i*2 + j) = -1;


				    *(turbo_trellis.mx_lastout + l*n_states*4 + i*4 + 2*j) = 0;					
				

				   *(turbo_trellis.mx_lastout + l*n_states*4 + i*4 + 2*j + 1) = 0 ; 


        			*(turbo_trellis.mx_laststat + l*n_states*2 +  i*2 + j ) = -1;

			}
	/* ======================================generate the first 0=====================================*/
	for( l = 0; l < 1; l ++)
		for( i = 0; i < 1; i ++)
			for (j=0; j<2; j++)
			{
				int2bin(i, tempstat, M_num_reg);
				dk_turbo = j;

				ak_turbo = (*(turbo_g.g_matrix+0)) * dk_turbo;

				for (k=1; k<turbo_g.K_num_col; k++)
				{
					ak_turbo = ak_turbo + (*(turbo_g.g_matrix+k)) * (*(tempstat+k-1));
				}

				ak_turbo = ak_turbo % 2;

				outbit = encode_bit(ak_turbo, tempstat);

				*(turbo_trellis.mx_nextstat + l*n_states*2 + i*2 + j) = bin2int(tempstat,M_num_reg);

				*(turbo_trellis.mx_laststat + ( l +1)*n_states*2 + bin2int(tempstat, M_num_reg)*2 + j ) = i;

				*(turbo_trellis.mx_lastout + (l + 1)*n_states*4 + bin2int(tempstat, M_num_reg)*4 + 2*j) =
					
					*(turbo_trellis.mx_nextout + l*n_states*4 + i*4 + 2*j ) =
					
						2*dk_turbo - 1;

				*(turbo_trellis.mx_lastout + (l+1)*n_states*4 + bin2int(tempstat, M_num_reg)*4 + 2*j + 1 ) = 
					
					*(turbo_trellis.mx_nextout + l*n_states*4 + i*4 + 2*j + 1 ) = 
					
						2*outbit - 1;
			}
	/*====================================The first 1 to M_num_reg - 1=================================*/
	for( l = 1 ; l < length -3 ; l ++)
		for( i = 0 ; i < n_states ; i ++)
			for( j = 0 ; j <= 1; j ++)
			{	
					if(	*(turbo_trellis.mx_laststat + l*n_states*2 +  i*2 + 0 ) != -1 ||	*(turbo_trellis.mx_laststat + l*n_states*2 +  i*2 + 1 ) != -1)
					{	
						int2bin(i, tempstat, M_num_reg);
							
						dk_turbo = j;

						ak_turbo = (*(turbo_g.g_matrix+0)) * dk_turbo;

						for (k=1; k<turbo_g.K_num_col; k++)
						{
							ak_turbo = ak_turbo + (*(turbo_g.g_matrix+k)) * (*(tempstat+k-1));
						}

						ak_turbo = ak_turbo % 2;

						outbit = encode_bit(ak_turbo, tempstat);

						*(turbo_trellis.mx_nextstat + l* n_states *2 + i*2 + j) = bin2int(tempstat,M_num_reg);

						*(turbo_trellis.mx_laststat + (l+1)*n_states*2 + bin2int(tempstat, M_num_reg)*2 + j ) = i;
			
						*(turbo_trellis.mx_lastout + (l + 1)*n_states*4 + bin2int(tempstat, M_num_reg)*4 + 2*j)  = *(turbo_trellis.mx_nextout + l* n_states *4 + i*4 + 2*j ) =
					
							2*dk_turbo - 1;

						*(turbo_trellis.mx_lastout + ( l + 1)*n_states*4 + bin2int(tempstat, M_num_reg)*4 + 2*j + 1 ) = 
					
						*(turbo_trellis.mx_nextout + l* n_states *4 + i*4 + 2*j + 1 ) = 
					
						2*outbit - 1;
					}
			
			}
			for( l =  length - 3 ; l <=length; l ++)
				
				for( i = 0 ; i < n_states ; i ++)
					
					for( j = 0 ; j <= 1; j ++)
				     {	
						*(turbo_trellis.mx_nextstat + l* n_states *2 + i*2 + j) = *(turbo_trellis.mx_laststat +  (length-l) *n_states*2 + i*2 + j);

						*(turbo_trellis.mx_laststat + l*n_states*2 + i*2 + j ) =  *(turbo_trellis.mx_nextstat + ( length-l)* n_states *2 + i*2 + j) ;

						*(turbo_trellis.mx_lastout + l*n_states*4 + i*4 + 2*j) =  *(turbo_trellis.mx_nextout + ( length - l)* n_states *4 + i*4 + 2*j );
					
					    *(turbo_trellis.mx_lastout + l*n_states*4 + i*4 + 2*j + 1) =  *(turbo_trellis.mx_nextout + (length - l)* n_states *4 + i*4 + 2*j + 1);

						*(turbo_trellis.mx_nextout + l*n_states*4 + i*4 + 2*j ) =   *(turbo_trellis.mx_lastout + (length -l)* n_states *4 + i*4 + 2*j );
					
						*(turbo_trellis.mx_nextout + l*n_states*4 + i*4 + 2*j + 1 ) =   *(turbo_trellis.mx_lastout + (length -l)* n_states *4 + i*4 + 2*j  + 1);
				  }

					
			 FILE *stream;	
	      if ((stream=fopen( "turbo_trellis.txt", "w" )) == NULL)
			  {
				printf("\nError! Can not open file\n");
				exit(1);
				} 

			for( l = 0; l <= length; l ++)
				for( i  = 0; i < n_states; i++ )
					for( j = 0; j <= 1; j ++)
					{	
						fprintf(stream," turbo_trellis.mx_nextstat[%d][%d][%d] = %d    " ,l,i,j,*(turbo_trellis.mx_nextstat + l* n_states *2 + i*2 + j) );

						fprintf(stream," turbo_trellis.mx_nextout[%d][%d][%d]= %d%d   " ,l,i,j,*(turbo_trellis.mx_nextout + l* n_states *4 + i*4 + 2*j + 0 ),*(turbo_trellis.mx_nextout + l* n_states *4 + i*4 + 2*j + 1 ) );

						fprintf(stream," turbo_trellis.mx_laststat[%d][%d][%d] = %d   " ,l,i,j,*(turbo_trellis.mx_laststat + l* n_states *2 + i*2 + j) );

						fprintf(stream," turbo_trellis.mx_lastout[%d][%d][%d] = %d%d  \n \n \n",l,i,j,*(turbo_trellis.mx_lastout + l* n_states *4 + i*4 + 2*j + 0 ),*(turbo_trellis.mx_lastout + l* n_states *4 + i*4 + 2*j + 1 )  );


					}
	 //printf("Trellis succeeds!\n");
     fclose(stream);
	free(tempstat);
}
/*---------------------------------------------------------------
FUNCTION: 
	void int2bin(int intstat, int *tempstat, int length)
	
DESCRIPTION:
	This function transfers a decimal integer to a binary sequence.

PARAMETERS:
	INPUT:
		intstat - Decimal integer needed to be changed.
		length - Length of binary sequence.
	OUTPUT:
		tempstat - Contains pointer to the binary sequence.

RETURN VALUE:
	None
---------------------------------------------------------------*/
void int2bin(int intstat, int *tempstat, int length)
{
	int i, temp;

	temp = intstat;

	for (i=length-1; i>=0; i--)
	{
		*(tempstat+i) = temp%2;
		temp = temp/2;
	}
}
/*---------------------------------------------------------------
FUNCTION: 
	int bin2int(int *binseq, int length)
	
DESCRIPTION:
	This function transfers a binary sequence to a decimal integer.

PARAMETERS:
	INPUT:
		binseq - Contains pointer to the binary sequence.
		length - Length of binary sequence.

RETURN VALUE:
	The decimal integer corresponding to the binary sequence.
---------------------------------------------------------------*/
int bin2int(int *binseq, int length)
{
	int i, j, temp;
	int sum = 0;

	for (i=0; i<length; i++)
	{
		temp = 1;

		for (j=1; j<=i; j++)
		{
			temp = temp * 2;
		}
		sum = sum + temp * (*(binseq+length-1-i));
	}

	return sum;
}
/*---------------------------------------------------------------
FUNCTION: 
	random_turbo()
	
DESCRIPTION:
	Generate a random number between 0 and 1.

RETURN VALUE:
	The random number between 0 and 1.
---------------------------------------------------------------*/
float random_turbo()
{
   float random_number;
   srand((int)time(0));
   random_number = (float)random(100)/100;
   return random_number;
}
/*---------------------------------------------------------------
FUNCTION: 
	void gen_rand_index(int length)
	
DESCRIPTION:
	This function generate the index for the random interleaver of two kinds; one is with terminated bits and the other without .

PARAMETERS:
	INPUT:
		length - Length of interleaver.

RETURN VALUE:
	None
---------------------------------------------------------------*/
void gen_rand_index(int length)
{
	int *index_random;
	float *tempindex;
	float tempmax;
	int selectedscr;
	int i, j;

	tempindex=(float *)malloc((length)*sizeof(float));

	if(tempindex == NULL)
		printf("Wrong!\n");

	index_random = index_randomintlvr;

	//printf("\nThe interleaver matrix with terminating bits is as follows:\n ");
	
	srand((int)time(0));
	for (i=0; i<length; i++)
	{
		*(tempindex+i)  = (float)random(100)/100;
	}

	for (i=0; i<length; i++)
	{
		tempmax = 0.0;

		for (j=0; j<length; j++)
		{
			if (*(tempindex+j) >= tempmax)
			{
				tempmax = *(tempindex+j);
				selectedscr = j;
			}
		}

		*(index_random+i) = selectedscr;
		*(reverse_index_randomintlvr + selectedscr) = i;
	   // printf(" %d ---> %d ", i,*(index_random+i));
		*(tempindex+selectedscr) = -1.0;
	}
	free(tempindex);
}
/*---------------------------------------------------------------
FUNCTION: 
	void randominterleaver_int(int *data_unintlvr, int *interleaverddata,
								 int length, int *index)
	
DESCRIPTION:
	Random interleaver of int.

PARAMETERS:
	INPUT:
		data_unintlvr - Contains pointer to data_unintlvr before interleavered.
		length - Length of interleaver.
		index - the pointer to index
	OUTPUT:
		interleaverddata - Contains pointer to interleavered data_unintlvr.

RETURN VALUE:
	None
---------------------------------------------------------------*/
void randominterleaver_int(int *data_unintlvr, int *interleaverddata,
						   int length, int * index)
{
	int i;
	int *index_random = index;

	for (i=0; i<length; i++)
	{
		*(interleaverddata+i) = *(data_unintlvr+ (*(index_random+i)));
		//printf("%d",	*(interleaverddata+i) );
	}


}
/*---------------------------------------------------------------
FUNCTION: 
	void randominterleaver_float(float *data_unintlvr, float *interleaverddata, int length, int *index)
	
DESCRIPTION:
	Random interleaver of float.

PARAMETERS:
	INPUT:
		data_unintlvr - Contains pointer to data_unintlvr before interleavered.
		length - Length of interleaver.
		index - the pointer to index
	OUTPUT:
		interleaverddata - Contains pointer to interleavered data_unintlvr.

RETURN VALUE:
	None
---------------------------------------------------------------*/
void randominterleaver_float(float *data_unintlvr, float *interleaverddata, int length, int *index)
{
	int i;
	int *index_random =  index ;

	for (i=0; i<length; i++)
	{
		*(interleaverddata+i) = *(data_unintlvr+ (*(index_random+i)));
	}

	/*for( i = 0 ; i < length ; i++)
		printf("\n %f \n",*(interleaverddata+i));*/
}
/*---------------------------------------------------------------
FUNCTION: 
	void random_deinterlvr_int(int *data_unintlvr, int *interleaverddata, int length, int* index)
	
DESCRIPTION:
	Random deinterleaver of int.

PARAMETERS:
	INPUT:
		data_unintlvr - Contains pointer to data_unintlvr before interleavered.
		length - Length of interleaver.
		index - the pointer to index

	OUTPUT:
		interleaverddata - Contains pointer to interleavered data_unintlvr.

RETURN VALUE:
	None
---------------------------------------------------------------*/
void random_deinterlvr_int(int *data_unintlvr, int *interleaverddata, int length, int * index)
{
	int i;
	int *index_random =  index;

   //printf("\n Process of deinterleaving! \n ");
	for (i=0; i<length; i++)
	{
		*(data_unintlvr+(*(index_random+i))) = *(interleaverddata+i);
		//printf("%d ",	*(data_unintlvr+(*(index_random+i)))  );
	}
   
	/*printf("\n After deinterleaving: \n");
	for (i=0; i<length; i++)
		printf(" %d ",	*(data_unintlvr+ i));*/

}
/*---------------------------------------------------------------
FUNCTION: 
	void random_deinterlvr_float(float *data_unintlvr, float *interleaverddata, int length, int *index)
	
DESCRIPTION:
	Random interleaver of float.

PARAMETERS:
	INPUT:
		data_unintlvr - Contains pointer to data_unintlvr before interleavered.
		length - Length of interleaver.
		index - the pointer to index
	OUTPUT:
		interleaverddata - Contains pointer to interleavered data_unintlvr.

RETURN VALUE:
	None
---------------------------------------------------------------*/
void random_deinterlvr_float(float *data_unintlvr, float *interleaverddata, int length, int *index)
{
	int i;
	int *index_random =  index;

	for (i=0; i<length; i++)
	{
		*(data_unintlvr+(*(index_random+i))) = *(interleaverddata+i);
		
	}

	/*for (i=0; i<length; i++)
	{
		printf("\n%f\n", *(data_unintlvr + i));
		
	}*/

}

/*---------------------------------------------------------------
FUNCTION: 
	void encoderm_turbo(int *source, int *send_turbo, int len_info)
	
DESCRIPTION:
	This function encode the source data and modulate them using BPSK.

PARAMETERS:
	INPUT:
		source - Contains pointer to source bits.
		len_info - Length of the input information bits.

	OUTPUT:
		send_turbo - Contains pointer to the sequence after encoding and modulation.

RETURN VALUE:
	None
---------------------------------------------------------------*/

void encoderm_turbo(int *source, int *send_turbo, int len_info)
{
	int i;
	int len_total = len_info + M_num_reg + 1;

	int *rsc1, *rsc2;
	int *with_terminate_bits;
	
	int *input2;

    with_terminate_bits = (int *)malloc(len_total*sizeof(int));

    for( i = 0; i < len_total; i++) /*adding terminating bits*/
	{	
		if(i < len_info)
			*(with_terminate_bits + i) = *(source + i);
		else 
			*(with_terminate_bits + i) = 0;
		//printf("%d  ",*(with_terminate_bits + i));
	}
	 
	rsc1 = (int *)malloc(2*len_total*sizeof(int));

	rsc2 = (int *)malloc(2*len_total*sizeof(int));

	input2 = (int *)malloc(len_total*sizeof(int));  /*input2 is used to load the elements after interleaving*

	/* the first RSC encoder */
	//printf("\n After the first encoding: the bits are as follows.\n");

	rsc_encode(with_terminate_bits, rsc1, len_total, up_state);

	//printf("\n The end of first rsc encoder. \n");

	//printf("\n After interleaving: the bits are as follows.\n");
	
	randominterleaver_int(with_terminate_bits, input2,len_total, index_randomintlvr);

	/* the second RSC encoder */
	//printf("\n After the second encoding: the bits are as follows.\n");
	rsc_encode(input2, rsc2, len_total, down_state);

	//printf("\n After turbo encoding: the bits are as follows.\n");
	for (i=0; i<len_total; i++)
	{
		*(send_turbo+3*i) = *(rsc1+2*i) *2 - 1;
		//printf("%d",*(send_turbo+3*i));
		*(send_turbo+3*i+1) = *(rsc1+2*i+1) *2 - 1;
		//printf("%d",*(send_turbo+3*i+1));
		*(send_turbo+3*i+2) = *(rsc2+2*i+1) *2 - 1;
		//printf("%d",*(send_turbo+3*i+2));
	}

	free(rsc1);
	free(rsc2);
	free(input2);
	free(with_terminate_bits);
}
/*---------------------------------------------------------------
FUNCTION: 
	void rsc_encode(int *source, int *rsc, int len_info, int *state)

DESCRIPTION:
	Encodes a block of data to a recursive systematic convolutional code.

PARAMETERS:
	INPUT:
		source - Contains pointer to the input data sequence together with terminating bits.
		state - The memory state of registers.
		len_total - Length of the input data sequencd plus terminating bits.
	OUTPUT:
		rsc - Contains pointer to the output data sequence.

RETURN VALUE:
	None
---------------------------------------------------------------*/
void rsc_encode(int *source, int *rsc, int len_total, int * state )  
{
	int i, j;

	int dk_turbo, ak_turbo, outbit;

	for (i=0; i<len_total; i++)			/* encoding bit by bit */
	{
		dk_turbo = *(source+i);

		ak_turbo = *(turbo_g.g_matrix+0) * dk_turbo;
		
		for (j=1; j<turbo_g.K_num_col; j++)
		{
			ak_turbo = ak_turbo + (*(turbo_g.g_matrix+j))*(*(state+j-1));
		}                              /*compute the feeding back bit*/

		ak_turbo = ak_turbo%2;

		outbit = encode_bit(ak_turbo, state); /*using the feeding back bit to compute the out put bit*/

		*(rsc+2*i) = dk_turbo;      /*systematic code*/
		//printf("%d",*(rsc+2*i));
		*(rsc+2*i+1) = outbit;      /* the other code*/
		//printf("%d",*(rsc+2*i+1));
	}				
}


/*---------------------------------------------------------------
FUNCTION: 
	int encode_bit(int inbit, int *stat)
	
DESCRIPTION:
	This function can get a coded bit from the memory state and the input bit.
	The memory state is accordingly modified.

PARAMETERS:
	INPUT:
		inbit - The input bit.
	OUTPUT & INPUT:
		stat - Contains pointer to memory state.

RETURN VALUE:
	The coded bit.
---------------------------------------------------------------*/
int encode_bit(int inbit, int *stat)
{
	int j;
	int output;

	output = (*(turbo_g.g_matrix+turbo_g.K_num_col+0)) * inbit;

	for (j=1; j<turbo_g.K_num_col; j++)
	{
		output = (output + (*(turbo_g.g_matrix+turbo_g.K_num_col+j)) * (*(stat+j-1)))%2;
	}

	for (j=turbo_g.K_num_col-2; j>0; j--)
	{
		*(stat+j)=*(stat+j-1);
	}

	*(stat+0) = inbit;
 
//printf("\n Now the states of registers are :\n");
//for( j = 0; j <= turbo_g.K_num_col - 2 ; j ++)
//	printf(" %d ", *(stat+j) );

	return output;
}

/*---------------------------------------------------------------
FUNCTION: 
	void gen_mx_punc()

DESCRIPTION:                                         1 1
	This function generate the puncture matrix.      0 1
	                                                 1 0
	                                            
PARAMETERS:

RETURN VALUE:
	None
---------------------------------------------------------------*/
void gen_mx_punc()
{
	*(mx_puncture_turbo + 0*2 + 0) = *(mx_puncture_turbo + 0*2 + 1) = *(mx_puncture_turbo +1*2 +1) = *(mx_puncture_turbo + 2*2 + 0) = 1;
	*(mx_puncture_turbo + 1*2 + 0) = *(mx_puncture_turbo + 2*2 + 1) = 0;
	return ;
}

/*---------------------------------------------------------------
FUNCTION: 
	puncture(int *data_unpunc, int length_unpunc, int *data_punctured)

DESCRIPTION:
	This function puncture the sending bits.

PARAMETERS:
	INPUT:
		data_unpunc - Contains pointer to the input sequence.
		length_unpunc - Length of "data_unpunc".
		mx_puncture_turbo - Puncture matrix.

	OUTPUT:
		data_punctured - Contains pointer to the punctured sequence.

RETURN VALUE:
	None
---------------------------------------------------------------*/
void puncture(int *data_unpunc, int length_unpunc, int *data_punctured)
{
	int i, j, k=0, temp_time=0, time;
		
	time = length_unpunc/6;
	
	for (temp_time=0; temp_time< time; temp_time++)
	
		for (j=0; j<2; j++)

			for (i=0; i<3;i++)
			
				if (*(mx_puncture_turbo+i*2+j))
				{
					*(data_punctured+k) = *(data_unpunc+temp_time*6+j*3+i);
				    // printf("%d",*(data_punctured+k));
					k++;
				}
	if(length_unpunc %6 !=0)
	{
		for (j=0; j<1; j++)

			for (i=0; i<3;i++)
			
				if (*(mx_puncture_turbo+i*2+j))
				{
					*(data_punctured+k) = *(data_unpunc+temp_time*6+j*3+i);
				    //printf("%d",*(data_punctured+k));
					k++;
				}
	}

}
/*---------------------------------------------------------------
FUNCTION: 
	depuncture(float *receive_punc, int length_punc, float *receive_depunced)
	
DESCRIPTION:
	This function fills in the punctured sequence.

PARAMETERS:
	INPUT:
		receive_punc - Contains pointer to the punctured input sequence.
		length_punc - Length of "receive_punc".

	OUTPUT:
		receive_depunced - Contains pointer to the filled sequence.

RETURN VALUE:
	None
---------------------------------------------------------------*/
void depuncture(float *receive_punc, int length_punc, float *receive_depunced)
{
	int i, j, k=0, temp_time = 0; 

	while( temp_time < length_punc/4)
	{
		for (j=0; j<2; j++)
		{
			for (i=0; i<3; i++)
			{
				if (*(mx_puncture_turbo+i*2+j))
				{
					*(receive_depunced+temp_time*6+j*3+i) = *(receive_punc+k);
					//printf("%f ",*(receive_depunced+temp_time*6+j*3+i));
					k++;

				}
				else
				{
					*(receive_depunced+temp_time*6+j*3+i) = 0;
					//printf("%f ",*(receive_depunced+temp_time*6+j*3+i));

				}
			}
		}
		
		temp_time ++;
	}

	if(length_punc%4 != 0)
		
		for (j=0; j<1; j++)
		{
			for (i=0; i<3; i++)
			{
				if (*(mx_puncture_turbo+i*2+j))
				{
					*(receive_depunced+temp_time*6+j*3+i) = *(receive_punc+k);
					//printf("%f ",*(receive_depunced+temp_time*6+j*3+i));
					k++;

				}
				else
				{
					*(receive_depunced+temp_time*6+j*3+i) = 0;
					//printf("%f ",*(receive_depunced+temp_time*6+j*3+i));

				}
			}
		}

		
}
/*---------------------------------------------------------------
FUNCTION: 
	void demultiplex(double *rec_turbo, int len_info, double *yk_turbo)
	
DESCRIPTION:
	Demultiplex the receiving data.

PARAMETERS:
	INPUT:
		rec_turbo - Contains pointer to receiving data.
		len_total - Length of the total input bits.

	OUTPUT:
		yk_turbo - Contains pointer to the sequence after demultiplexing.

RETURN VALUE:
	None
---------------------------------------------------------------*/
void demultiplex(float *rec_turbo, int len_total, float *yk_turbo)
{
	int i;

    float *info2, *inted_info2;

	info2=(float *)malloc( len_total*sizeof(float));
	inted_info2=(float *)malloc(len_total*sizeof(float));
	
	/* for information bits */
	for(i=0; i<len_total; i++)                          /* The first 2*len_total is the received information bits(origianl infomation bits addeing AWGN) and the ECC from the first encoder*/
	{                                                   /*The rest one should be the sequence after the first decoding plus the ECC from the second encoder, however, here we only have ECC*/ 
		*(info2+i) = *(yk_turbo+2*i) = *(rec_turbo+3*i);
		*(yk_turbo+2*i+1) = *(rec_turbo+3*i+1);
		*(yk_turbo+2*len_total+2*i+1) = *(rec_turbo+3*i+2);
		//printf("\n yk_turbo[%d] =  %f,  yk_turbo[%d] = %f, yk_turbo[%d] = %f, yk_turbo[%d] = %f    \n ", 2*i,*(yk_turbo+2*i),2*i+1,*(yk_turbo+2*i+1),2*(len_total)+2*i,*(yk_turbo+2*(len_total)+2*i),2*(len_total)+2*i + 1,*(yk_turbo+2*(len_total)+2*i + 1));
	}
    
	randominterleaver_float(info2, inted_info2, len_total, index_randomintlvr);
	for (i=0; i<len_total; i++)
	{
		*(yk_turbo+2*(len_total)+2*i) = *(inted_info2+i);
		//printf("\n yk_turbo[%d] =  %f,  yk_turbo[%d] = %f, yk_turbo[%d] = %f, yk_turbo[%d] = %f    \n ", 2*i,*(yk_turbo+2*i),2*i+1,*(yk_turbo+2*i+1),2*(len_total)+2*i,*(yk_turbo+2*(len_total)+2*i),2*(len_total)+2*i + 1,*(yk_turbo+2*(len_total)+2*i + 1));
	}
    

	free(info2);
}
/*---------------------------------------------------------------
FUNCTION: 
	Log_MAP_decoder(float *recs_turbo, float *La_turbo, float *LLR_all_turbo, int len_total)

DESCRIPTION:
	Log-MAP decoder which caculate the LLR of input sequence.

PARAMETERS:
	INPUT:
		recs_turbo - Scaled received bits.
		La_turbo - A priori information for the current decoder, 
			scrambled version of extrinsic information of the previous decoder.
		len_total - Length of the input data sequence.
	OUTPUT:
		LLR_all_turbo - The caculated LLR information sequence.

RETURN VALUE:
	None
---------------------------------------------------------------*/
void Log_MAP_decoder(float *recs_turbo, float *La_turbo, float *LLR_all_turbo, int len_total)
{
	int i, j;

	float *alpha_Log, *beta_Log, *gama_Log;

	float *tempmax;
	float *temp0, *temp1;
	float tempx, tempy;

	alpha_Log=(float *)malloc(n_states*(len_total +1)*sizeof(float));
	
	beta_Log=(float *)malloc(n_states*(len_total + 1)*sizeof(float));
	
	gama_Log=(float *)malloc(n_states*(len_total)*2*sizeof(float));
	temp0=(float *)malloc(n_states*sizeof(float));
	
	temp1=(float *)malloc(n_states*sizeof(float));
	
	/*===== Initialization of alpha_Log and beta_Log =====*/
	*(alpha_Log + 0) = 0;
	*(beta_Log + len_total*n_states) = 0;

	for (i=1; i<n_states; i++)
	{
		*(alpha_Log + 0*n_states + i) = (float)-INFTY;

		*(beta_Log + len_total*n_states + i) = (float)-INFTY;
	
	}

	/*========compute Gama_Log========*/
	for (i=0; i<len_total; i++)			
	{
		if( i < len_total - (M_num_reg + 1))
		{
			for (j=0; j<n_states; j++)		/* j->k */
			{
				if(*(turbo_trellis.mx_nextstat + i*2*n_states + j*2 + 0 ) != -1)
					*(gama_Log + i*n_states*2 + j*2 +0) 
						= (*(recs_turbo+2*i)*(*(turbo_trellis.mx_nextout + i*4*n_states + j*4 + 0 )) + *(recs_turbo+2*i+1)*(*(turbo_trellis.mx_nextout + i*4*n_states + j*4 + 1)))/2 - *(La_turbo+i)/2;
				else
					*(gama_Log + i*n_states*2 + j*2 +0) = (float)SATELITE;

				if(*(turbo_trellis.mx_nextstat + i*2*n_states + j*2 + 1 ) != -1)
					*(gama_Log+ i*n_states*2 + j*2 + 1)
						= (*(recs_turbo+2*i)*(*(turbo_trellis.mx_nextout + i*4*n_states + j*4 + 2))  + *(recs_turbo+2*i+1)*(*(turbo_trellis.mx_nextout + i*4*n_states + j*4 + 3)))/2  +  *(La_turbo+i)/2;
				else
					
					 *(gama_Log+ i*n_states*2 + j*2 + 1)  = (float)SATELITE;
				//printf("\n gama_Log[%d][%d][%d] = %f,gama_Log[%d][%d][%d] = %f  \n", i,j,*(turbo_trellis.mx_nextstat + i*2*n_states + j*2 + 0 ),*(gama_Log+ i*n_states*2 + j*2 ),i,j,*(turbo_trellis.mx_nextstat + i*2*n_states + j*2 + 1),*(gama_Log+ i*n_states*2 + j*2 +1));
			 }
		}


		else 
		{	
			for (j=0; j<n_states; j++)		/* j->k */
			{
					if(*(turbo_trellis.mx_nextstat + i*2*n_states + j*2 + 0 ) != -1)
						*(gama_Log + i*n_states*2 + j*2 + 0) = ( *(recs_turbo + 2*i) * ( *(turbo_trellis.mx_nextout + i*4*n_states + j*4 + 0) )  +  *( recs_turbo + 2*i + 1)*( *(turbo_trellis.mx_nextout + i*4*n_states +j*4 + 1 ) ) )/2;
					else
						*(gama_Log + i*n_states*2 + j*2 +0) = (float)SATELITE;

					if(*(turbo_trellis.mx_nextstat + i*2*n_states + j*2 + 1 ) != -1)
						*(gama_Log+i*n_states*2 +j*2+1 )= (*(recs_turbo+2*i)*(*(turbo_trellis.mx_nextout + i*4*n_states + j*4 + 2)) + *(recs_turbo+2*i+1)*(*(turbo_trellis.mx_nextout + i*4*n_states + j*4 + 3)))/2;

					else	
						  *(gama_Log+ i*n_states*2 + j*2 + 1)  = (float)SATELITE;
					//printf("\n gama_Log[%d][%d][%d] = %f,gama_Log[%d][%d][%d] = %f  \n", i,j,*(turbo_trellis.mx_nextstat + i*2*n_states + j*2 + 0 ),*(gama_Log+ i*n_states*2 + j*2 ),i,j,*(turbo_trellis.mx_nextstat + i*2*n_states + j*2 + 1),*(gama_Log+ i*n_states*2 + j*2 +1));
					
			}
		}
	}
	/*========Trace forward, compute Alpha_Log========*/
	for (i=1; i<len_total + 1; i++)
	{
		for (j=0; j<n_states; j++)		/* centers at j*/
		{

			if(*(turbo_trellis.mx_laststat + i*2*n_states + 2*j + 0 ) != -1)
				tempx = *(alpha_Log+ (i-1)*n_states +  *(turbo_trellis.mx_laststat + i*2*n_states + j*2 + 0 ))  +  *(gama_Log + (i-1)*n_states*2 + *(turbo_trellis.mx_laststat + i*2*n_states + j*2 + 0 )*2  + 0 ) ;
			else 
				tempx  = (float) -INFTY;

			if( *(turbo_trellis.mx_laststat + i*2*n_states + 2*j + 1 ) != -1)
				tempy = *(alpha_Log+ (i-1)*n_states +  *(turbo_trellis.mx_laststat + i*2*n_states + j*2 + 1 )) +  *(gama_Log + (i-1)*n_states*2 + *(turbo_trellis.mx_laststat + i*2*n_states + j*2 + 1 )*2  + 1 ) ;
			else	
				tempy  = (float) -INFTY;

			if( *(turbo_trellis.mx_laststat + i*2*n_states + 2*j + 0 ) == -1 && *(turbo_trellis.mx_laststat + i*2*n_states + 2*j + 1 ) == -1)
				*(alpha_Log + i*n_states + j) = (float)SATELITE ;

			else 
				*(alpha_Log + i*n_states + j)  =   E_algorithm(tempx, tempy);
			//printf("\nalpha_Log[%d][%d] = %f \n", i,j,*(alpha_Log + i*n_states + j));
		}

	}

	/*========Trace backward, compute Beta_Log========*/
	for (i=len_total - 1; i>=0; i--)
	{
		for (j=0; j<n_states; j++)		/*centers at j*/
		{   
			if( *(turbo_trellis.mx_nextstat + i*2*n_states + 2*j + 0 ) != -1 )
					tempx =  *(beta_Log + (i + 1)*n_states + *(turbo_trellis.mx_nextstat + (i*n_states*2) + j*2 + 0) ) +  *(gama_Log + i*n_states*2 + j*2 + 0 );
			else  
				    tempx =  (float)-INFTY;
			
			if( *(turbo_trellis.mx_nextstat + i*2*n_states + 2*j + 1 ) != -1 )
		
              	   tempy = *(beta_Log + (i + 1)*n_states + *(turbo_trellis.mx_nextstat + (i*n_states*2) + j*2 + 1) ) + *(gama_Log + i*n_states*2 + j*2 + 1);
			else 
				   tempy = (float)-INFTY;

			if( *(turbo_trellis.mx_nextstat + i*2*n_states + 2*j + 0 ) == -1 && *(turbo_trellis.mx_nextstat + i*2*n_states + 2*j + 1 ) == -1)
					*(beta_Log + i*n_states + j) = (float)SATELITE ;

			else 
				*(beta_Log + i*n_states + j) = E_algorithm(tempx, tempy);

				//printf("\n beta_Log[%d][%d] = %f \n", i,j,*(beta_Log + i*n_states + j));
		
		}
	}

	/*===Compute the soft output,log-likelihood ratio of symbols in the frame===*/
	for (i=0; i<len_total; i++)	
	{
		for (j=0; j<n_states; j++)
		{
				
			if(*(turbo_trellis.mx_nextstat + i*n_states*2 + 2*j + 0) != -1 )
			
				 *(temp0+j) =  *(alpha_Log + i*n_states + j) + *(gama_Log +i*n_states*2 + j*2 + 0 ) + *(beta_Log + ( i + 1)*n_states + *(turbo_trellis.mx_nextstat + i*n_states*2 + 2*j + 0));
			else
				*(temp0+j) = (float)SATELITE;

        
			if(*(turbo_trellis.mx_nextstat + i*n_states*2 + 2*j + 1)!= -1 )
			
				*(temp1+j) =  *(alpha_Log + i*n_states + j) + *(gama_Log +i*n_states*2 + j*2 + 1 ) + *(beta_Log + ( i + 1)*n_states + *(turbo_trellis.mx_nextstat + i*n_states*2 + 2*j + 1));
			else	
				*(temp1+j) = (float)SATELITE;
		}
		*(LLR_all_turbo+i) = E_algorithm_seq(temp1, n_states) - E_algorithm_seq(temp0, n_states);
		//printf("\nLLR_all_turbo[%d] = %f\n",i,*(LLR_all_turbo+i) );
	}

	free(alpha_Log);
	free(beta_Log);
	free(gama_Log);
	free(temp0);
	free(temp1);

}
/*---------------------------------------------------------------
FUNCTION: 
	E_algorithm(float x, float y)

DESCRIPTION:
	Compute: log(exp(x) + exp(y)) = max(x,y)+log(1+exp(-|y-x|)) 
	where log(1+exp(-|y-x|)) can be implemented in a lookup table.
					
PARAMETERS:
	INPUT:
		x - One number.
		y - The other number.

RETURN VALUE:
	log(exp(x) + exp(y)).
---------------------------------------------------------------*/
float E_algorithm(float x, float y) // log(exp(x) + exp(y)) = max(x,y)+f(-|y-x|) 
{
	float temp = (y-x)>0? (y-x):(x-y);
	int i;

	if (temp>=4.3758)
	{
		temp = 0;
	}
	else
	{
		for (i=0; i<16 && temp>=lookup_index_Log_MAP[i]; i++)
		{
			;
		}
		temp = (float)lookup_table_Log_MAP[i-1];
	}
	
	return ( (x>y?x:y) + temp );
}	

/*---------------------------------------------------------------
FUNCTION: 
	E_algorithm_seq(float *data_seq, int length)

DESCRIPTION:
	E_algorithm for a data_seq sequence.

PARAMETERS:
	INPUT:
		data_seq - Contains pointer to the input sequence.
		length - Length of "data_seq".

RETURN VALUE:
	The result of E_algorithm for this data_seq sequence.
---------------------------------------------------------------*/
float E_algorithm_seq(float *data_seq, int length)
{
	int j = 0,key1 = 0, key2 = 0,find1 = 0, find2 = 0;
	
	float temp;

    while(j < length)
	{	
		if(*(data_seq + j ) != (float)SATELITE)
		{	
			if(find1 == 1)
			{
				find2 = 1;
				key2 = j;
			}
			else
			{ 
				find1 = 1;
				key1 = j;
			}
		}

		if(find1*find2 == 1)
		{	
			break;
			j++;
		}
		else j++;
	}

	if( find1*find2 == 0) 
		return *(data_seq + key1);
	
	else
	{

		temp = E_algorithm(*(data_seq+key1), *(data_seq+key2));

		while(j < length)
		{   

			find2 = 0;

			while(j < length)
			{	
				if(*(data_seq + j ) != (float)SATELITE)
				{	
					find2 = 1;
					key2 = j;
					j++;
					break;
				}

				j++;
			}
			if(find2 == 1)

				temp = E_algorithm(temp, *(data_seq+key2));
	
			else  return temp;
		}

	  return temp;
	}
   
}

/*---------------------------------------------------------------
FUNCTION: 
	void decision(double *LLR_seq, int length, int *output)
	
DESCRIPTION:
	Make the final decision.

PARAMETERS:
	INPUT:
		LLR_seq - The LLR_seq information sequence.
		length - Length of "LLR_seq".
	OUTPUT:
		output - Contains pointer to the output data sequence.

RETURN VALUE:
	None
---------------------------------------------------------------*/
void decision(float *LLR_seq, int length, int *output)
{
	int i;
   
	for (i=0; i<length; i++)
	{
		if (*(LLR_seq+i) < 0)
			*(output+i) = 0;
		else
			*(output+i) = 1;

		if(fabs(*(LLR_seq+i)) < THRESHOLD )//We need to store the locations of might-be erasures;
			locations_of_erasure.push_back(*(index_randomintlvr+i));
	}		
}