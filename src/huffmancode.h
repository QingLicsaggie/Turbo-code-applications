/*---------------------------------------------------------------
* Copyright NOV 2010 CSE@TAMU
* All rights reserved.
*
* huffmancode.h
*
* This script implements the functions used in turbo coding and decoding.
*
* Version: 1.0
* Programmed By Qing(Tsing) Li
* Last updated date: Nov 2010
---------------------------------------------------------------*/
#include<iostream>
using namespace std;
using std::fixed;
#include<iomanip>
using std::setprecision;


#include<cstdio>
#include<memory>
#include <fstream>
#include<cstring>
#include<queue>
#include<string> 
#include<vector>
/*-----------------------------------------------------------------*/
/*---------------Global value--------------------------------------*/
struct node
{
    int times;                                  /*frequency*/
    struct node *child1;
    struct node *child2;
	string code;                                /*huffmancode*/
	char c;                                     /*character*/
    
    friend bool operator<( const node &a,const node &b )       /*overload*/
    {
         return ( a.times < b.times );
    }
    friend bool operator>( const node &a1,const node &b1 )
    {
         return ( a1.times > b1.times );
    }
};
                                          
int findChar( char );    
void generate_code();
void encode();                                        
void decode();                            
