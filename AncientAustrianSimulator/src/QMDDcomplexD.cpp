/***************************************************************

Complex number defnitions and routines for QMDD using
doubles for the real and imaginary part of a complex number.

January 28, 2008
Michael Miller
University of Victoria
Victoria, BC 
CANADA V8W 3P6
mmiller@cs.uvic.ca

****************************************************************/

/****************************************************************

The basic idea is that the required complex values are 
stored in a lookup table.	

The value 0 is always in slot 0 and the value 1 is always in slot
1 so that for those two values the index corresponds to the value.

Current implementation uses simple linear searching to find a value.

QMDDinit (in QMDDpackage.c which is the initialization routine that
must be called before the other package routines are used) invokes 
QMDDinitCtable()

*****************************************************************/

#define DEFINE_COMPLEX_H_VARIABLES
#include "QMDDcomplex.h"

using mpfr::mpreal;


void QMDDpause(void);

/**************************************

    Routines
    
**************************************/

mpreal QMDDcos(int fac, long div) {
	return cos(Pi * fac/div);
}

mpreal QMDDsin(int fac, long div) {
	return sin(Pi*fac/div);
}


mpreal tmp,tmp2;
complex tmp_c;

bool my_compare(complex& x, complex& y) {
	mpfr_sub(tmp.mpfr_ptr(), x.r.mpfr_ptr(), y.r.mpfr_ptr(), MPFR_RNDN);

	bool sign = mpfr_signbit(tmp.mpfr_ptr());
	mpfr_abs(tmp.mpfr_ptr(), tmp.mpfr_ptr(), MPFR_RNDN);
	if(mpfr_cmp(tmp.mpfr_ptr(), Ctol.mpfr_ptr()) < 0) {
		mpfr_sub(tmp.mpfr_ptr(), x.i.mpfr_ptr(), y.i.mpfr_ptr(), MPFR_RNDN);
		sign = mpfr_signbit(tmp.mpfr_ptr());
		mpfr_abs(tmp.mpfr_ptr(), tmp.mpfr_ptr(), MPFR_RNDN);
		if(mpfr_cmp(tmp.mpfr_ptr(), Ctol.mpfr_ptr()) < 0) {
			return 0;
		} else {
			return sign;
		}
	} else {
		return sign;
	}
}



int Ceq(complex& x, complex& y) {
	mpfr_sub(tmp.mpfr_ptr(), x.r.mpfr_ptr(), y.r.mpfr_ptr(), MPFR_RNDN);
	mpfr_abs(tmp.mpfr_ptr(), tmp.mpfr_ptr(), MPFR_RNDN);
	if(mpfr_cmp(tmp.mpfr_ptr(), Ctol.mpfr_ptr()) < 0) {
		mpfr_sub(tmp.mpfr_ptr(), x.i.mpfr_ptr(), y.i.mpfr_ptr(), MPFR_RNDN);
		mpfr_abs(tmp.mpfr_ptr(), tmp.mpfr_ptr(), MPFR_RNDN);
		if(mpfr_cmp(tmp.mpfr_ptr(), Ctol.mpfr_ptr()) < 0) {
			return 1;
		}
	}
	return 0;
}


void Cprint(complex c, std::ostream &os)
{
	if(c.r >= 0)
	  os << " ";
	os << c.r;
	if (c.i > 0)
	   os << "+" << c.i << "i";
	if (c.i < 0)
	  os << c.i << "i";
}

void Cprint(complex c)
// print a complex value
{
	std::ostringstream oss;
	Cprint(c, oss);
	std::cout << oss.str();  
}

mpreal angle(int a)
// computes angle for polar coordinate representation of Cvalue(a)
{
	complex ca;
  ca=Cvalue(a);

  tmp = - Ctol;

  std::map<unsigned int, mpreal>::iterator it = Cmag.find(a);

  mpreal res =ca.r / it->second;
  res = acos(res);
  if(ca.i < tmp) {
	  tmp = Pi * 2;
	  res = tmp - res;
  }   else {
  }
  return res;
}

int Cdiv(int ai,int bi); /* prototype */

int Cgt(int a, int b)
{  
  complex ca,cb;
  if(a==b) return(0);
  ca=Cvalue(a);
  cb=Cvalue(b);
  
  if (a == 0)
    return(1);
  if (b == 0)
    return(0);
  
  std::map<unsigned int, mpreal>::iterator it = Cmag.find(b);
  tmp = it->second + Ctol;

  it = Cmag.find(a);
  if(it->second > tmp) {
	  return 1;
  }
  tmp = it->second + Ctol;
  it = Cmag.find(b);
  if(it->second > tmp) {
	  return(0);
  }
  //CHANGED by pN 120831
  it = Cangle.find(a);

  tmp = it->second + Ctol;
  it = Cangle.find(b);
  int ret_val = (tmp < it->second);
  return ret_val;
}

int Cgt_new(int a, int b)
{  
  complex ca,cb;
  if(a==b) return(0);
  ca=Cvalue(a);
  cb=Cvalue(b);

  std::map<unsigned int, mpreal>::iterator it = Cangle.find(a);

  tmp = it->second + Ctol;
  it = Cangle.find(b);

  if(tmp < it->second) {
	  return(1);
  }

  it = Cmag.find(b);
  tmp = it->second + Ctol;
  it = Cmag.find(a);

  int ret_val = (it->second > tmp);
  return ret_val;
}

int Clt(int a, int b)
// analogous to Cgt
{
  complex ca,cb;
  if(a==b) return(0);
  ca=Cvalue(a);
  cb=Cvalue(b);

  std::map<unsigned int, mpreal>::iterator it = Cmag.find(b);

  tmp = it->second + Ctol;

  it = Cmag.find(a);
  if(it->second < tmp) {
	  return(1);
  }
  tmp = it->second + Ctol;
  it = Cmag.find(b);
  if(it->second < tmp) {
	  return(0);
  }
  it = Cangle.find(a);
  tmp = it->second + Ctol;
  it = Cangle.find(b);
  int ret_val = (tmp > it->second);
  return ret_val;
}

complex Cmake(mpreal r,mpreal i)
// make a complex value
{
  complex c;
  c.r = mpreal(r);
  c.i = mpreal(i);

  return(c);
}

complex CmakeOne(void)
{
  complex c;
  c.r = mpreal(1);

  c.i = mpreal(0);
  return c;
}

complex CmakeZero(void)
{
	  complex c;
	  c.r = mpreal(0);
	  c.i = mpreal(0);
	  return c;
}

complex CmakeMOne(void)
{
	  complex c;
	  c.r = mpreal(-1);
	  c.i = mpreal(0);
	  return c;
}

mpreal Qmake(int a, int b,int c)
// returns the complex number equal to (a+b*sqrt(2))/c
// required to be compatible with quadratic irrational-based 
// complex number package
{
	mpreal res;
	res = (a+b*sqrt(2))/c;
	return res;
}

void QMDDinitCtable(void)
// initialize the complex value table and complex operation tables to empty
{
  int i,j;
  
  Ctentries=0;

  if(VERBOSE) printf("\nDouble complex number package initialized\n\n");
}

void QMDDcomplexInit(void)
// initialization
{

	mpreal::set_default_prec(PREC);
	tmp = mpreal(0);
	tmp2 = mpreal(0);
	tmp_c.i = mpreal(0);
	tmp_c.r = mpreal(0);

	Ctol = mpreal(1.0e-12);

	mpreal mag1,mag2;
	mag1 = mpreal(0);
	Cmag.insert(std::pair<unsigned int, mpreal>(0, mag1));
	mag2 = mpreal(1);
	Cmag.insert(std::pair<unsigned int, mpreal>(1, mag2));

	QMDDinitCtable();
}

void QMDDcvalue_table_list(void)
// print the complex value table entries
{
/*  int i;
  
  printf("\nComplex value table: %d entries\n",Ctentries);
  std::cout << "index value Magnitude Angle 1) radian 2) degree" << std::endl;
  for(i=0;i<Ctentries;i++)
  {
    std::cout << i << " ";
    //printf("%d ",i);
    Cprint(Ctable[i]);
    //printf("  || %f %f\n",Cmag[i],angle(i));
    if(i!=0&&i%100==0) QMDDpause();
    std::cout << " || " << Cmag[i]<< " " << angle(i);
    std::cout << " " <<(angle(i)*180/3.141592654) << std::endl;
 }*/
}


int Clookup(complex c)
// lookup a complex value in the complex value table
// if not found add it
// this routine uses linear searching
{
  int i;


  std::map<complex, unsigned int, complex_cmp>::iterator it;
  it = Ctable2.find(c);
  if(it != Ctable2.end()) {
	  return it->second;
  }


  if(Ctentries >= 2147483640) {
	  std::cout << "ERROR too many complex numbers" << std::endl; exit(-4);
  }

  Ctentries++;
  i=Ctentries-1;

  complex new_c;

  new_c.r = mpreal(c.r);
  new_c.i = mpreal(c.i);
  Ctable[i] = new_c;

  Ctable2[new_c]=i;

  mpreal new_mag, new_angle;
  new_mag = mpreal();
  new_angle = mpreal();

  new_mag = c.r * c.r;
  new_angle = c.i*c.i;
  new_mag = sqrt(new_mag + new_angle);
  std::pair<unsigned int, mpreal> my_pair = std::make_pair(i, new_mag);

  Cmag.insert(my_pair);

  Cangle[i] = angle(i);
  return(i);
}

complex Conj(complex c)
// return complex conjugate
{
	c.i = -c.i;
  return(c);
}


// basic operations on complex values
// meanings are self-evident from the names
// NOTE arguments are the indices to the values 
// in the complex value table not the values themselves

int Cnegative(int a)
{
  complex c;
  c=Cvalue(a);
  complex c2;

  c2.r = mpreal();
  c2.i = mpreal();

  c2.r = -c.r;
  c2.i = -c2.i;

  int ret_val = Clookup(c);

  return ret_val;
}

int Cadd(int ai,int bi)
{
  complex a,b;
  int t;
  
  if(ai==0) return(bi); // identity cases
  if(bi==0) return(ai);

  std::map<unsigned long, unsigned int>::iterator it;
  unsigned long key = ai;
  key = (key << 32) | bi;
  it = cta.find(key);
  if(it != cta.end()) {
	  return it->second;
  }
  
  a=Cvalue(ai); // if new compute result
  b=Cvalue(bi); 

  tmp_c.r = a.r + b.r;
  tmp_c.i = a.i + b.i;

  t=Clookup(tmp_c); // save result
  cta.insert(std::pair<unsigned long, unsigned int>(key, t));
  key = bi;
  key = (key << 32) | ai;
  cta.insert(std::pair<unsigned long, unsigned int>(key, t));
  return(t);
}

int Csub(int ai,int bi)
{
  complex a,b;
  int t;
  
  if(bi==0) return(ai); // identity case
  
  std::map<unsigned long, unsigned int>::iterator it;
  unsigned long key = ai;
  key = (key << 32) | bi;
  it = cts.find(key);
  if(it != cts.end()) {
	  return it->second;
  }
  
  a=Cvalue(ai);  // if new compute result
  b=Cvalue(bi);

  tmp_c.r = a.r - b.r;
  tmp_c.i = a.i - b.i;

  
  t=Clookup(tmp_c); // save result
  cts.insert(std::pair<unsigned long, unsigned int>(key, t));
  return(t);
}

int Cmul(int ai,int bi)
{
  complex a,b;
  int t;
  
  if(ai==1) return(bi); // identity cases
  if(bi==1) return(ai);
  if(ai==0||bi==0) return(0);
  
  std::map<unsigned long, unsigned int>::iterator it;
  unsigned long key = ai;
  key = (key << 32) | bi;
  it = ctm.find(key);
  if(it != ctm.end()) {
	  return it->second;
  }

  
  a=Cvalue(ai); // if new compute result
  b=Cvalue(bi);

  tmp_c.r = a.r * b.r - a.i*b.i;
  tmp_c.i = a.r * b.i + a.i * b.r;

  t=Clookup(tmp_c); // save result

  ctm.insert(std::pair<unsigned long, unsigned int>(key, t));
  key = bi;
  key = (key << 32) | ai;
  ctm.insert(std::pair<unsigned long, unsigned int>(key, t));

  return(t);
}

int CintMul(int a,int bi)
{
  complex r;
  r=Cvalue(bi);

  tmp_c.r = r.r * a;
  tmp_c.i = r.i * a;
  int t = (Clookup(r));
  return t;
}

int Cdiv(int ai,int bi)
{
  complex a,b;
  int t;
  long double d;
  
  if(ai==bi) return(1); // equal case
  if(ai==0) return(0); // identity cases
  if(bi==1) return(ai);
  
  std::map<unsigned long, unsigned int>::iterator it;
  unsigned long key = ai;
  key = (key << 32) | bi;
  it = ctd.find(key);
  if(it != ctd.end()) {
	  return it->second;
  }

  a=Cvalue(ai); // if new compute result
  b=Cvalue(bi);

  if(b.i == 0)
  {
	  tmp_c.r = a.r/b.r;
	  tmp_c.i = a.i / b.r;
  } else {
	  tmp2 = mpreal();

	  tmp = b.r * b.r + b.i * b.i;

	  tmp_c.r = (a.r*b.r + a.i*b.i)/tmp;

	  tmp_c.i = (a.i * b.r - a.r * b.i)/tmp;
  }
  t=Clookup(tmp_c); // save result

  ctd.insert(std::pair<unsigned long, unsigned int>(key, t));

  return(t);
}

void QMDDmakeRootsOfUnity(void)
{
  int i;
  CTa[0]=1;

  mpreal r,s;
  r = mpreal();

  r = Pi * 2 / Radix;

  s = mpreal(r);
  r = cos(r);
  s = sin(s);

  complex c = Cmake(r,s);
  CTa[1] = Clookup(c);

  for(i=2;i<Radix;i++)
    CTa[i]=Cmul(CTa[i-1],CTa[1]);
}

/// by PN: returns the absolut value of a complex number
int CAbs(int a)
{
  int b;
  complex s;
  
  if (a<2) return a; // trivial cases 0/1
  
    s=Cvalue(a);
  //printf("CAbs: "); Cprint(s); printf(" is ");

   std::map<unsigned int, mpreal>::iterator it = Cmag.find(a);

   tmp_c.r = it->second;
   tmp_c.i = 0;
   b = Clookup(tmp_c);
  //Cprint(r);   printf("\n");
  return(b);
}

///by PN: returns whether a complex number has norm 1
int CUnit(int a)
{
 /// BETA 121017
 
 if (a<2)
   return a;

 std::map<unsigned int, mpreal>::iterator it = Cmag.find(a);

 tmp = it->second + Ctol;

 if (tmp < 1) {
	 return 0;
 }
 else {
    return 1;
 }
}

std::set<unsigned int> complex_entries;
std::set<QMDDnodeptr> visited_nodes;

void addToComplexTable(QMDDedge edge) {

	if(!QMDDterminal(edge)) {
		unsigned int before = visited_nodes.size();
		visited_nodes.insert(edge.p);
		if(before != visited_nodes.size()) {
			if(edge.p->e[0].w > 2) {
				complex_entries.erase(edge.p->e[0].w);
			}
			if(edge.p->e[1].w > 2) {
				complex_entries.erase(edge.p->e[1].w);
			}
			if(edge.p->e[2].w > 2) {
				complex_entries.erase(edge.p->e[2].w);
			}
			if(edge.p->e[3].w > 2) {
				complex_entries.erase(edge.p->e[3].w);
			}
			addToComplexTable(edge.p->e[0]);
			addToComplexTable(edge.p->e[1]);
			addToComplexTable(edge.p->e[2]);
			addToComplexTable(edge.p->e[3]);
		}
	}
}

void cleanCtable(std::vector<QMDDedge> save_edges) {

	complex_entries.clear();
	std::vector<complex> complex_entries2;
	std::map<unsigned int, complex>::iterator it;
	for(it = Ctable.begin(); it != Ctable.end(); it++) {
		complex_entries.insert(it->first);
	}
	visited_nodes.clear();


	for(std::vector<QMDDedge>::iterator it = save_edges.begin(); it != save_edges.end(); it++) {
		QMDDedge e = *it;
		complex_entries.erase(e.w);
		addToComplexTable(e);
	}
	complex_entries.erase(0);
	complex_entries.erase(1);
	complex_entries.erase(2);

	std::set<unsigned int>::iterator it2;
	std::map<unsigned int, mpreal>::iterator it3;
	for(it2 = complex_entries.begin(); it2 != complex_entries.end(); it2++) {
		complex_entries2.push_back(Ctable[*it2]);
	}
	std::vector<complex>::iterator it4;
	for(it4 = complex_entries2.begin(); it4 != complex_entries2.end(); it4++){
		Ctable2.erase(*it4);
	}

	for(it2 = complex_entries.begin(); it2 != complex_entries.end(); it2++) {
		Ctable.erase(*it2);
		it3 = Cmag.find(*it2);
		Cmag.erase(it3);
		it3 = Cangle.find(*it2);
		Cangle.erase(it3);
	}

	QMDDinitComputeTable();
	ctm.clear();
	cta.clear();
	cts.clear();
	ctd.clear();
}
