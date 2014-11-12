#include"huffmancode.h"

using namespace std;
using std::fixed;
using std::setprecision;


/*--------------------Global variables----------------------------------------*/
string text;                                 /*original file*/
string codes;
string  decoded_codes;                 
vector<node> code;	                          /*encoded file*/                     
vector<node> result;            
extern vector<int> locations_of_erasure;
/*---------------------------------------------------------------
* Copyright NOV 2010 CSE@TAMU
* All rights reserved.
*
* huffmancode.cpp
*
* This script implements the functions used in turbo coding and decoding.
*
* Version: 1.0
* Programmed By Qing(Tsing) Li
* Last updated date: Nov 2010
---------------------------------------------------------------*/

/*======================================================
FUNCTION: void generate_code()
DESCRIPTION: generating huffman codes for the input text file.
	INPUT:  none
	OUTPUT: The characters and corresponding codes are stored at the result. 
RETURN VALUE
    None
======================================================*/

void generate_code()
{
    priority_queue< node,vector<node>,greater<node> > que;     //Use priority_queue to construct huffman code
	int i;
	string pathName ;


	cout<<"Please input the text address"<<endl;
	cin>>pathName;

	ifstream in( pathName.c_str(), ios::in );                  //load the text.

	while( !in.eof() )
	{
		text +=in.get();
	}

	for( i=0; i<text.length(); i++ )                           //store all appearing character and accounting the frequency.
	{
		int sub=findChar(text[i]);                            //look up this code appearing or not
		if( sub == -1 )
		{
			node p;                                           //if not, add a new node
            p.times = 1; 
			p.c = text[i];
			p.code = "";
            p.child1 = NULL;
            p.child2 = NULL;
            code.push_back(p); 
		}
		else
		{
			code[sub].times++;                                //if so, times ++
		}
	}

	while( que.size() > 0 )
            que.pop();
        
        for( i = 0;i < code.size();i++ )                   //Initiate the priory queue. 
			que.push( code[i] );

	
        node *t1,*t2,*t3;
        while( que.size() >= 2 )                           //generating the tree
        {
            t1 = new node;
            t2 = new node;

            *t1 = que.top();                               //pop out two least frequent nodes.
            que.pop();

            *t2 = que.top();
            que.pop();   

            t3 = new node;
			t3->c = '~';
            t3->times = t1->times + t2->times;
            t3->child1 = t1;
            t3->child2 = t2;
            que.push( *t3 );                                //regenerate a new node.
            
        }

		node root = que.top();                              //the root.
		que.pop();

		queue<node*> queue;
		queue.push(&root);

		while( !queue.empty() )                              //begin encode, BFS
		{
			node *temp=queue.front();
			queue.pop();
			if( temp->child1!=NULL )
			{
				temp->child1->code=temp->code+"0";           //left one,code connects'0'
				queue.push(temp->child1);
			}
			if( temp->child2!=NULL )
			{
				temp->child2->code=temp->code+"1";           //right one£¬code connects'1'
				queue.push(temp->child2);
			}
			if( temp->child1==NULL && temp->child2==NULL )   
			{
				result.push_back(*temp);
			}
		}

    return ;
                
}
/*-----------------------------------------------------
FUNCTION: int findChar( char c )       

DESCRIPTION: find character c appearing or not

	INPUT:  character c

	OUTPUT:  if appearing   return the corresponding location.
	         otherwise, return the invalid location -1
RETURN VALUE
		
		NONE
------------------------------------------------------*/
int findChar( char c )                                      
{
	for( int i=0; i < code.size(); i++ )
	{
		if( code[i].c == c )
			return i;
	}
	return -1;
}


/*-----------------------------------------------------
FUNCTION: void encode()

DESCRIPTION  Huffman code function

	INTPUT:  Global parameter text
	
	OUTPUT:   Reuslts are stored at code.txt

RETURN VALUE
			NONE.
-----------------------------------------------------*/
void encode()                                   
{
	int i, j;

	generate_code();
	
	ofstream output( "code.txt", ios::out );

	for( i=0; i<text.length(); i++ )
	{
		for( j=0; j<result.size(); j++ )
		{
			if( result[j].c == text[i] )                      /*look up the corresponding code*/
			{
				codes+=result[j].code;
				break;
			}
		}

	output<<codes;
	}
}

/*---------------------------------------------------------------
FUNCTION£º void decode( string textName )   

DESCRIPTION£º huffman decoding function.

	INPUT: None.
	
	OUTPUT: the decoded text file

RETURN VALUE	
----------------------------------------------------------------*/
void decode( )                               
{
	int i,j=0, length = 0;
	string ch;
	bool flag = false;
	string temp;  

	ofstream output( "decoded.txt", ios::app);
	
	temp = decoded_codes[j];


	while( j < decoded_codes.length() )
	{

		for( i=0; i<result.size(); i++ )
		{                                                   /*looking for corresponding code*/
			if( result[i].code==temp )
			{
				ch +=result[i].c;                               
				flag=true;
				length ++;
				break;
			}
		}

		if( flag )
		{
			temp="";
             
			j++; 
			
			temp = decoded_codes[j];

	    	flag=false;
		}				/*successfully findly and reset the temp*/
		else 
		{  
			j++;
            
			temp += decoded_codes[j];
		}
	}
	if(!flag)
	{	
		cout<<endl<<" Sorry! We can not restore it completely."<<endl;

		cout<<endl<<" Part of the file is as follows are stored in decoded.txt:"<<endl;

		output<<ch;

	    output<<endl<<"========================================================="<<endl;
		return ;
	}
	
	else
	{
		output<<ch;
		cout<<ch;
		ch.~string();
		output<<endl<<"========================================================="<<endl;
	}
	output.close();

}