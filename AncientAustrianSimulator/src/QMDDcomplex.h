#ifndef QMDDcomplex_H
#define QMDDcomplex_H
#include <gmp.h>
#include <mpfr.h>
#include <map>
#include <set>
#include <vector>
#include "mpreal.h"
/***************************************

    Basic type definition

***************************************/

using mpfr::mpreal;

typedef struct
{
//	long double r,i;
   mpreal r;
   mpreal i;
} complex;


#define PREC 500

#include "QMDDpackage.h"

#include <ostream>
#include <math.h>

#ifndef DEFINE_COMPLEX_H_VARIABLES
#define EXTERN_C	extern 
#else
#define EXTERN_C	
#endif

static mpreal Ctol; // = 1.0e-12;

EXTERN mpreal Pi;					// Pi is defined using asin function in QMDDinit routine

/***************************************

	Complex Value Lookup Table
	
	The table holds the unique values needed
	in the QQDD.  It stores one value for each
	conjugate pair (the positive iinary).

***************************************/
void Cprint(complex, std::ostream&);
void Cprint(complex); // print a complex value to STD_OUT


bool my_compare(complex& x, complex& y);

struct complex_cmp
{
    bool operator() ( complex x, complex y ) const {
    	return my_compare(x,y);
    }
};


EXTERN_C int Ctentries;					 // number of complex table entries
EXTERN_C std::map<unsigned int, complex> Ctable; //[COMPLEXTSIZE];    // value
EXTERN_C std::map<complex, unsigned int, complex_cmp> Ctable2;
EXTERN_C std::map<unsigned int, mpreal> Cmag; //mpfr_t /*long double*/ Cmag[COMPLEXTSIZE];  // magnitude to avoid repeated computation
EXTERN_C std::map<unsigned int, mpreal> Cangle; //mpfr_t /*long double*/ Cangle[COMPLEXTSIZE];// angle to avoid repeated computation
EXTERN_C int CTa[MAXRADIX];				 // complex table positions for roots of unity

/*********************************************

Complex computation tables.  These tables save
result of computations over complex values to
avoid recomputation later.

  Cta - addition
  Cts - subtraction
  Ctm - multiplication
  Ctd - division

*********************************************/

EXTERN_C std::map<unsigned long, unsigned int> cta,cts,ctm,ctd;


/**********************************************

Compute trig functions for angle factor*Pi/div
Note use of cosl and sinl for long double computation

**********************************************/


mpreal QMDDcos(int fac, long div);
mpreal QMDDsin(int fac, long div);

/***********************************************

Tolerance for testing equality of complex values

***********************************************/



#define Cvalue(x) Ctable[x]
#define Dzero 0.0

int Ceq(complex x, complex y);
/**************************************

    Routines
    
**************************************/


mpreal angle(int); // computes angle for polar coordinate representation

int Cgt(int,int); // greater than
int Clt(int, int); // analogous to Cgt

complex Cmake(mpreal, mpreal); // make a complex value
complex CmakeOne(void);  // make +1
complex CmakeZero(void); // make  0
complex CmakeMOne(void); // make -1

mpreal Qmake(int,int,int);
// returns the complex number equal to (a+b*sqrt(2))/c
// required to be compatible with quadratic irrational-based 
// complex number package

void QMDDinitCtable(void); // initialize the complex value table and complex operation tables to empty
void QMDDcomplexInit(void); // initialization


void QMDDcvalue_table_list(void); // print the complex value table entries
int Clookup(complex); // lookup a complex value in the complex value table; if not found add it

complex Conj(complex); /// return complex conjugate


// basic operations on complex values
// meanings are self-evident from the names
// NOTE arguments are the indices to the values 
// in the complex value table not the values themselves

int Cnegative(int);
int Cadd(int,int);
int Csub(int,int);
int Cmul(int,int);
int CintMul(int,int); // multiply by an integer
int Cdiv(int,int);
void QMDDmakeRootsOfUnity(void);
int CAbs(int); /// by PN: returns the absolut value of a complex number
int CUnit(int a); ///by PN: returns whether a complex number has norm 1

#endif
