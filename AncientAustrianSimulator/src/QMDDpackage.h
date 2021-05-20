/**************************************************************************************

    package definitions and types
    
    Michael Miller
    University of Victoria
    mmiller@cs.uvic.ca
    
    Date: July 2008
    
**************************************************************************************/
#ifndef QMDDpackage_H
#define QMDDpackage_H

#include "external.h"

#define QMDDversion "QMDD Package V.R1 September 2015\n"

// problem parameter limits 

#define MAXSTRLEN 11
#define MAXN 300			// max no. of inputs
#define MAXRADIX 2    		// max logic radix                                		   
#define MAXNEDGE 4 			// max no. of edges = MAXRADIX^2
#define MAXNODECOUNT 2000000 	// max number of nodes in a QMDD for counting 			   
#define GCLIMIT1 250000   	// first garbage collection limit                          
#define GCLIMIT_INC 0 		// garbage collection limit increment                      
							// added to garbage collection limit after each collection
#define MAXND 8    			// max n for display purposes
#define MAXDIM 256           	// max dimension of matrix for printing, (should be 2^MAXND)
#define NBUCKET 8192     	// no. of hash table buckets; must be a power of 2 		   
#define HASHMASK 8191   	// must be nbuckets-1                 					   
#define CTSLOTS 16384  		// no. of computed table slots 							   
#define CTMASK  16383  		// must be CTSLOTS-1   
#define COMPLEXTSIZE 100000  	// complex table size
#define COMPLEXTMASK 127   	// complex table index mask   (not used anywhere?!)
#define TTSLOTS 2048		// Toffoli table slots
#define TTMASK 2047			// must be TTSLOTS-1
#define MAXREFCNT 4000000		// max reference count (saturates at this value)
#define MAXPL 65536			// max size for a permutation recording

#define DYNREORDERLIMIT 500	// minimum value for dynamic reordering limit
#define VERBOSE 0

#include <stdint.h>
#include <iostream>
#include "QMDDcomplex.h"
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// edge and node definitions 

typedef struct QMDDnode *QMDDnodeptr; 

typedef struct QMDDedge
{
   QMDDnodeptr p;  		// edge pointer 											 
   unsigned int w;          	// index of weight edge in complex value table 		
   int sentinel;		// in 64-bit architecture sizeof(QMDDedge) is a multiple of 8 (8 Byte = 64 bit).
				// So sizeof() is 16 with AND without this sentinel. thus when using memcmp identical nodes may not be identified. 
}  QMDDedge;

typedef struct QMDDnode
{
   QMDDnodeptr next;  // link for unique table and available space chain 
   unsigned int ref;  // reference count 												 
   unsigned char v;   // variable index (nonterminal) value (-1 for terminal)
   unsigned int renormFactor; // factor that records renormalization factor
   char ident,diag,block,symm,c01;        // flag to mark if vertex heads a QMDD for a special matrix
   char computeSpecialMatricesFlag;	  // flag to mark whether SpecialMatrices are to be computed
   //QMDDedge e[0]; 	  	// edges out of this node - variable so must be last in structure 
   QMDDedge e[MAXNEDGE];	// when calling malloc in QMDDgetnode
}  QMDDnode;


// list definitions for breadth first traversals (e.g. printing)  
typedef struct ListElement *ListElementPtr;

typedef struct ListElement
{
   int w,cnt;
   int line[MAXN];
   QMDDnodeptr p;
   ListElementPtr next;
}  ListElement;

// computed table definitions 

typedef enum{add,mult,kronecker,reduce,transpose,conjugateTranspose,transform,c0,c1,c2,none,norm,createHdmSign,findCmnSign,findBin,reduceHdm, renormalize} CTkind; // compute table entry kinds 

typedef struct CTentry// computed table entry defn 										 
{			
  QMDDedge a,b,r;     // a and b are arguments, r is the result 						 
  CTkind which;       // type of operation 												 
} CTentry;

typedef struct TTentry // Toffoli table entry defn
{
  int n,m,t,line[MAXN];
  QMDDedge e;
} TTentry;

typedef struct CircuitLine
{
  char input[MAXSTRLEN];
  char output[MAXSTRLEN];
  char variable[MAXSTRLEN];
  char ancillary;
  char garbage;

} CircuitLine;

typedef struct QMDDrevlibDescription // circuit description structure
{
  int n,ngates,qcost,nancillary,ngarbage;
  QMDDedge e,totalDC;
  CircuitLine line[MAXN];
  char version[MAXSTRLEN];
  char inperm[MAXN],outperm[MAXN];
  char ngate,cgate,tgate,fgate,pgate,vgate,kind[7],dc[5],name[32],no[8],modified;
} QMDDrevlibDescription;


typedef complex QMDD_matrix[MAXRADIX][MAXRADIX];

/* GLOBALS */

EXTERN QMDD_matrix Nm,Vm,VPm,Sm,Rm,Hm,Zm,ZEROm,Qm;



/***************************************

    Global variables

***************************************/
EXTERN void cleanCtable(std::vector<QMDDedge> save_edges);


#ifndef DEFINE_VARIABLES
EXTERN int Radix;  // radix (default is 2) 

EXTERN int Nedge;				// no. of edges (default is 4) 
#endif 

EXTERN QMDDnodeptr Avail;			// pointer to available space chain 

EXTERN ListElementPtr Lavail;		// pointer to available list elements for breadth first searchess

EXTERN QMDDnodeptr QMDDtnode;		// pointer to terminal node 

EXTERN QMDDedge QMDDone,QMDDzero; 	// edges pointing to zero and one QMDD constants 


EXTERN long QMDDorder[MAXN];		// variable order initially 0,1,... from bottom up | Usage: QMDDorder[level] := varible at a certain level
EXTERN long QMDDinvorder[MAXN];	// inverse of variable order (inverse permutation) | Usage: QMDDinvorder[variable] := level of a certain variable

EXTERN long QMDDnodecount;			// counts active nodes 
EXTERN long QMDDpeaknodecount;                 // records peak node count in unique table

EXTERN long Ncount;				// used in QMDD node count - very naive approach 
EXTERN QMDDnodeptr Nlist[MAXNODECOUNT];

EXTERN long Nop[6];				// operation counters

EXTERN long CTlook[20],CThit[20];	// counters for gathering compute table hit stats 

EXTERN long UTcol, UTmatch;			// counter for collisions / matches in hash tables   

EXTERN int GCcurrentLimit;			// current garbage collection limit 

EXTERN int ActiveNodeCount;		// number of active nodes 

EXTERN int Active[MAXN];			// number of active nodes for each variable 

#ifndef DEFINE_VARIABLES
EXTERN int GCswitch;           // set switch to 1 to enable garbage collection 

EXTERN int Smode;				// S mode switch for spectral transformation
							// Smode==1 0->+1 1->-1; Smode==0 0->0 1->1

EXTERN int RMmode;				// Select RM transformation mode
							// forces mod Radix arithmetic
							
						
EXTERN int MultMode;			// set to 1 for matrix - vector multiplication
#endif

EXTERN QMDDedge QMDDnullEdge;		    // set in QMDDinit routine

EXTERN int PermList[MAXPL];		// array for recording a permutation

// for sifting
#ifndef DEFINE_VARIABLES
EXTERN int RenormalizationNodeCount;	// number of active nodes that need renormalization (used in QMDDdecref) 
EXTERN int blockMatrixCounter;	        // number of active nodes that represent block matrices (used in QMDDincref, QMDDdecref)
EXTERN char globalComputeSpecialMatricesFlag; // default value for computeSpecialMatricesFlag of newly created nodes (used in QMDDmakeNonterminal)
EXTERN int dynamicReorderingTreshold;

EXTERN int largestRefCount;
#endif
/*******************************************

	Unique Tables (one per input variable)
	
*******************************************/

EXTERN QMDDnodeptr Unique[MAXN][NBUCKET];

/****************************************************

    Compute Table (only one for all operation types)
    
****************************************************/

EXTERN CTentry CTable[CTSLOTS];

/****************************************************

    Toffoli gate table

*****************************************************/

EXTERN TTentry TTable[TTSLOTS];

/****************************************************

    Identity matrix table
    
****************************************************/

EXTERN QMDDedge QMDDid[MAXN];

/****************************************************

Variable labels

****************************************************/

EXTERN int  Nlabel;		// number of labels
EXTERN char Label[MAXN][MAXSTRLEN];  // label table

/****************************************************

Output file for writing permutation 

****************************************************/

EXTERN FILE *outfile;



#ifdef DEFINE_VARIABLES

// initalize global variables
int Radix = 2;				// radix (default is 2) 
int Nedge=4;				// no. of edges (default is 4) 
int GCswitch = 1;           // set switch to 1 to enable garbage collection
int Smode = 1;				// S mode switch for spectral transformation
					// Smode==1 0->+1 1->-1; Smode==0 0->0 1->1
int RMmode =0;				// Select RM transformation mode
int MultMode = 0;			// set to 1 for matrix - vector multiplication

// for sifting

int RenormalizationNodeCount = 0;	// number of active nodes that need renormalization (used in QMDDdecref) 
int blockMatrixCounter = 0;	        // number of active nodes that represent block matrices (used in QMDDincref, QMDDdecref)
int globalComputeSpecialMatricesFlag = 1; // default value for computeSpecialMatricesFlag of newly created nodes (used in QMDDmakeNonterminal)
int dynamicReorderingTreshold = DYNREORDERLIMIT;

int largestRefCount = 0;

#else


#endif


/***************************************

	Public Macros

***************************************/


#define __NormA__ // must be __NormA__ for leftmost-nonzero normalization and __NormB__ for leftmost absolut-value normalization

#define QMDDterminal(e) (e.p==QMDDtnode) // checks if an edge points to the terminal node

#define QMDDedgeEqual(a,b) ((a.p==b.p)&&(a.w==b.w)) // checks if two edges are equal





/** Routines 
 * 
 */

void throwException(const char*, int);

void QMDDpause(void);
void QMDDdebugnode(QMDDnodeptr); // for debugging purposes - not normally used
ListElementPtr QMDDnewListElement(void);
void QMDDprint(QMDDedge,int);
void QMDD2dot(QMDDedge,int, std::ostream&, QMDDrevlibDescription);
QMDDedge QMDDmultiply(QMDDedge,QMDDedge);
QMDDedge QMDDadd(QMDDedge,QMDDedge);
QMDDedge QMDDkron(QMDDedge,QMDDedge);
void QMDDdecref(QMDDedge);
void QMDDincref(QMDDedge);
QMDDedge QMDDident(int,int);
QMDDedge QMDDmvlgate(complex[MAXRADIX][MAXRADIX],int ,int[]);
void TTinsert(int,int,int,int[],QMDDedge);
QMDDedge TTlookup(int,int,int,int[]);
void QMDDgarbageCollect(void);
QMDDedge QMDDtranspose(QMDDedge); //prototype
void QMDDmatrixPrint2(QMDDedge); // prototype
QMDDedge QMDDnormalize(QMDDedge);
void QMDDinitGateMatrices(void);
void QMDDcheckSpecialMatrices(QMDDedge);
QMDDedge CTlookup(QMDDedge,QMDDedge,CTkind);
void CTinsert(QMDDedge,QMDDedge,QMDDedge,CTkind);
void QMDDinitComputeTable(void);
QMDDedge QMDDutLookup(QMDDedge);
QMDDedge QMDDmakeNonterminal(short,QMDDedge[]);
void QMDDinit(int verbose);
void QMDDdotExport(QMDDedge basic, int n, char outputFilename[], QMDDrevlibDescription circ, int show);
void QMDDstatistics(void);
QMDDedge QMDDconjugateTranspose(QMDDedge a);
QMDDedge QMDDtrace(QMDDedge a, unsigned char var, char remove[], char all);
void QMDDprintActive(int n);
#endif
