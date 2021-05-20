#include <iostream>
#include <QMDDcore.h>
#include <QMDDpackage.h>
#include <QMDDcomplex.h>
#include <stdlib.h>
#include <set>
#include <map>
#include <boost/program_options.hpp>
#include <sys/time.h>
#include <vector>
#include "mpreal.h"

using namespace std;
using mpfr::mpreal;

int nqubits;
int nq;
int* line;
int *measurements;
QMDDrevlibDescription circ;
int gatecount = 0;
mpreal epsilon;

std::map<QMDDnodeptr, mpreal> probs;

mpreal assignProbs(QMDDedge& e) {
	map<unsigned int, mpreal>::iterator it2;
	map<QMDDnodeptr, mpreal>::iterator it = probs.find(e.p);
	if(it != probs.end()) {
		it2 = Cmag.find(e.w);
		return it2->second * it2->second * it->second;
	}
	mpreal sum;
	if(QMDDterminal(e)) {
		sum = mpreal(1);
	} else {
		sum = assignProbs(e.p->e[0]) + assignProbs(e.p->e[1]) + assignProbs(e.p->e[2]) + assignProbs(e.p->e[3]);
	}

	probs.insert(std::pair<QMDDnodeptr, mpreal>(e.p, sum));

	it2 = Cmag.find(e.w);

	return it2->second * it2->second * sum;
}

void QMDDmeasureAll() {
	cout << "Measure all qubits: " << endl;

	QMDDedge e;

	map<QMDDnodeptr, mpreal>::iterator it;
	map<unsigned int, mpreal>::iterator it2;

	probs.clear();
	e = circ.e;

	mpreal p,p0,p1,tmp,w;
	p = assignProbs(e);

	if(abs(p -1) > epsilon) {
		cout << "Numerical error occurred during simulation: |alpha0|^2 + |alpha1|^2 = " << p<< ", but should be 1!"<<endl;
		exit(1);
	}
	QMDDedge cur = e;
	QMDDedge edges[4];
	QMDDedge prev = QMDDone;
	int prev_index = 0;
	for(int i = QMDDinvorder[e.p->v]; i >= 0;--i) {
		cout << "  -- measure qubit " << circ.line[cur.p->v].variable << ": " << flush;
		it = probs.find(cur.p->e[0].p);
		it2 = Cmag.find(cur.p->e[0].w);
		p0 = it->second * it2->second * it2->second;
		it2 = Cmag.find(cur.p->e[1].w);
		it = probs.find(cur.p->e[1].p);
		p0 += it->second * it2->second * it2->second;

		it = probs.find(cur.p->e[2].p);
		it2 = Cmag.find(cur.p->e[2].w);
		p1 = it->second * it2->second * it2->second;
		it = probs.find(cur.p->e[3].p);
		it2 = Cmag.find(cur.p->e[3].w);
		p1 += it->second * it2->second * it2->second;

		it2 = Cmag.find(cur.w);
		mpreal tmp = it2->second * it2->second;
		p0 *= tmp;
		p1 *= tmp;

		if(abs(p0+p1-1) > epsilon) {
			cout << "Numerical error occurred during simulation: |alpha0|^2 + |alpha1|^2 = " << p0+p1 << ", but should be 1!"<<endl;
			exit(1);
		}

		cout << "p0 = " << p0 << ", p1 = " << p1 << flush;
		mpreal n = mpreal(rand()) / RAND_MAX;

		if(n < p0) {
			cout << " -> measure 0" << endl;
			measurements[cur.p->v] = 0;
			edges[0]=cur.p->e[0];
			edges[1]=cur.p->e[1];
			edges[2]=edges[3]=QMDDzero;

			if(i == QMDDinvorder[e.p->v]) {
				it2 = Cmag.find(e.w);
				w = it2->second;
				e = QMDDmakeNonterminal(cur.p->v, edges);
				it2 = Cmag.find(e.w);
				w *= it2->second;
				e.w = 1;
				cur = e;
			} else {
				it2 = Cmag.find(cur.w);
				w = it2->second;
				cur = QMDDmakeNonterminal(cur.p->v, edges);
				it2 = Cmag.find(cur.w);
				w *= it2->second;
				cur.w = 1;
				prev.p->e[prev_index] = cur;
			}

			if(cur.p->e[0].w != 0) {
				prev_index = 0;
			} else {
				prev_index = 1;
			}
			prev=cur;
			cur = prev.p->e[prev_index];
			p0 = w / sqrt(p0);
			mpreal zero = mpreal(0);
			int a = Clookup(Cmake(p0,zero));

			cur.w = Cmul(cur.w, a);
		} else {
			cout << " -> measure 1" << endl;
			measurements[cur.p->v] = 1;
			edges[2]=cur.p->e[2];
			edges[3]=cur.p->e[3];
			edges[0]=edges[1]=QMDDzero;

			if(i == QMDDinvorder[e.p->v]) {
				it2 = Cmag.find(e.w);
				w = it2->second;
				e = QMDDmakeNonterminal(cur.p->v, edges);
				it2 = Cmag.find(e.w);
				w *= it2->second;
				e.w = 1;
				cur = e;
			} else {
				it2 = Cmag.find(cur.w);
				w = it2->second;
				cur = QMDDmakeNonterminal(cur.p->v, edges);
				it2 = Cmag.find(cur.w);
				w *= it2->second;
				cur.w = 1;
				prev.p->e[prev_index] = cur;
			}

			if(cur.p->e[2].w != 0) {
				prev_index = 2;
			} else {
				prev_index = 3;
			}

			prev=cur;
			cur = prev.p->e[prev_index];

			p1 = w / sqrt(p1);
			cur.w = Cmul(cur.w, Clookup(Cmake(p1, mpreal(0))));
		}
	}

	cout << endl;
}




int QMDDmeasureOne(int index) {
	cout << "Measure qubit " << circ.line[index].variable << ": " << flush;

	QMDDedge e;

	while(QMDDorder[index] != nq) {
		QMDDswap(QMDDorder[index]+1);
	}

	std::map<QMDDnodeptr, mpreal>::iterator it;
	std::map<unsigned int, mpreal>::iterator it2;

	probs.clear();
	e = circ.e;

	mpreal p, p0, p1, tmp, f;
	p = assignProbs(e);


	if(abs(p - 1) > epsilon) {
		cout << "Numerical error occurred during simulation: |alpha0|^2 + |alpha1|^2 != " << p<< ", but should be 1!"<<endl;
		exit(1);
	}
	QMDDedge edges[4];
	int meas;

	it = probs.find(e.p->e[0].p);
	it2 = Cmag.find(e.p->e[0].w);

	p0 = it->second * it2->second * it2->second;
	it2 = Cmag.find(e.p->e[1].w);
	it = probs.find(e.p->e[1].p);
	p0 += it->second * it2->second * it2->second;

	it = probs.find(e.p->e[2].p);
	it2 = Cmag.find(e.p->e[2].w);
	p1 = it->second * it2->second * it2->second;
	it2 = Cmag.find(e.p->e[3].w);
	it = probs.find(e.p->e[3].p);
	p1 += it->second * it2->second * it2->second;

	it2 = Cmag.find(e.w);
	tmp = it2->second * it2->second;
	p0 *= tmp;
	p1 *= tmp;

	cout << "p0 = " << p0 << ", p1 = " << p1 << flush;

	mpreal n = mpreal(rand()) / RAND_MAX;

	if(n < p0) {
		edges[0] = e.p->e[0];
		edges[1] = e.p->e[1];
		edges[2] = edges[3] = QMDDzero;
		meas = 0;

		f = mpreal(1) / sqrt(p0);
	} else {
		edges[0] = edges[1] = QMDDzero;
		edges[2] = e.p->e[2];
		edges[3] = e.p->e[3];
		meas = 1;

		f = mpreal(1) / sqrt(p1);
	}
	QMDDedge tmp_e = QMDDmakeNonterminal(nq, edges);

	tmp_e.w = Cmul(Clookup(Cmake(f, mpreal(0))), tmp_e.w);
	tmp_e.w = Cmul(tmp_e.w, circ.e.w);
	QMDDincref(tmp_e);
	QMDDdecref(circ.e);
	circ.e = tmp_e;

	cout << " -> measure " << meas << endl;

	return meas;
}

int max_active;



int complex_limit = 16000;

void QMDDsimulate(string fname)
{
	FILE *infile;

	int first,i;

	QMDDedge e,f,olde;

	// get name of input file, open it and attach it to file (a global)
	infile=openTextFile((char*)fname.c_str(),'r');
	if(infile == NULL) {
		cout << "Error: failed to open circuit file!"<< endl;
		exit(3);
	}

	//Read header from infile
	circ=QMDDrevlibHeader(infile);

    cout << "Simulate " << fname << " (requires " << circ.n << " qubits):" << endl;

	circ.ngate=circ.cgate=circ.tgate=circ.fgate=circ.pgate=circ.vgate=0;
	circ.qcost=circ.ngates=0;

	first=1;
	e = QMDDone;
	QMDDedge edges[4];
	edges[1]=edges[3]=QMDDzero;

	for(int p=0;p<circ.n;p++) {
		if(circ.line[p].ancillary=='0') {
			edges[0] = e;
			edges[2] = QMDDzero;
		} else if(circ.line[p].ancillary=='1') {
			edges[0] = QMDDzero;
			edges[2] = e;
		} else {
			cout << "Error: the input vector is not correctly specified!" << endl;
			exit(2);
		}
		e = QMDDmakeNonterminal(p, edges);
	}
	QMDDincref(e);
	first = 0;

	cout << "Start simulation (input vector = ";
	for(int i=circ.n-1; i >= 0; --i) {
		cout << circ.line[i].ancillary;
	}
	cout << "):" << endl;

	while(1) // read gates
	{
		f=QMDDreadGate(infile,&circ);
		if(f.p==NULL) break;

		if(first) // first gate in circuit
		{
			first = 0;
			e=f;
			QMDDincref(e);
		}
		else // second and subsequent gates
		{
			olde=e;
			e=QMDDmultiply(f,e); // multiply QMDD for gate * QMDD for circuit to date

			QMDDincref(e);
			QMDDdecref(olde);
			circ.e = e;
		}
		if(ActiveNodeCount > max_active) {
			max_active = ActiveNodeCount;
		}
		gatecount++;
		if(gatecount %1000 == 0) {
			cout << "  -- processed " << gatecount << " gates" << endl;
		}

		QMDDgarbageCollect();
		if(Ctable.size() > 16000) {
			vector<QMDDedge> v;
			v.push_back(circ.e);
			cleanCtable(v);
			v.clear();
			complex_limit = 2 * Ctable.size();
		}
	}

	for(i=0;i<circ.n;i++) circ.outperm[i]=i;

	skip2eof(infile); // skip rest of input file

	circ.e=e;

	i=0;
	if(circ.ngate) circ.kind[i++]='N';
	if(circ.cgate) circ.kind[i++]='C';
	if(circ.tgate) circ.kind[i++]='T';
	if(circ.fgate) circ.kind[i++]='F';
	if(circ.pgate) circ.kind[i++]='P';
	if(circ.vgate) circ.kind[i++]='V';
	circ.kind[i]=0;

	i=0;
	if(circ.nancillary>0) circ.dc[i++]='C';
	if(circ.ngarbage>0) circ.dc[i++]='G';
	circ.dc[i]=0;

	fclose(infile);

	measurements = new int[circ.n];
	QMDDmeasureAll();
	cout << "Measured: ";
	for(int i = circ.n-1; i >= 0; --i) {
		cout << measurements[i];
	}
	cout << " (starting from the most significant bit)" << endl;

	delete[] measurements;
}


set<QMDDnodeptr> visited;

void countNodes(QMDDnodeptr p) {
	if(p == QMDDone.p) {
		return;
	}
	if(visited.find(p) != visited.end()) {
		return;
	}
	visited.insert(p);
	countNodes(p->e[0].p);
	countNodes(p->e[1].p);
	countNodes(p->e[2].p);
	countNodes(p->e[3].p);
}

void apply_gate(QMDD_matrix& m, int count, int target) {
	gatecount++;

	QMDDedge e,f;
	f = QMDDmvlgate(m, circ.n, line);

	e = QMDDmultiply(f, circ.e);
	QMDDincref(e);
	QMDDdecref(circ.e);
	circ.e = e;

	QMDDgarbageCollect();

	if(ActiveNodeCount > max_active) {
		max_active = ActiveNodeCount;
	}

	if(Ctable.size() > complex_limit) {
		vector<QMDDedge> v;
		v.push_back(circ.e);

		cleanCtable(v);
		v.clear();
		if(complex_limit < 2*Ctable.size()) {
			complex_limit *= 2;
		}
	}
}

int inverse_mod(int a, int n) {
    int t = 0;
	int newt = 1;
    int r = n;
    int newr = a;
    int h;
    while(newr != 0) {
        int quotient = r / newr;
        h = t;
        t = newt;
        newt = h - quotient * newt;
        h = r;
        r = newr;
        newr = h - quotient * newr;
    }
    if(r > 1) {
    	cout << "ERROR: a is not invertible" << endl;
    }
    if(t < 0) {
    	t = t + n;
    }
    return t;
}

mpreal q_r, q_i;

void add_phi(int n, int c1, int c2) {
	int q;
	int controls = 0;
	if(c1 >= 0) controls++;
	if(c2 >= 0) controls++;
	for(int i = nqubits; i>= 0; --i) {
		q = 1;
		int fac = 0;
		for(int j = i; j >= 0; --j) {
			if((n >> j)&1) {
				fac |= 1;
			}
			q*=2;
			fac = fac << 1;
		}
		if(c1 >= 0) {
			line[nq-c1] = 1;
		}
		if(c2 >= 0) {
			line[nq-c2] = 1;
		}
		line[nq-(1+2*nqubits-i)] = 2;
		q_r = QMDDcos(fac, q);
		q_i = QMDDsin(fac, q);

		Qm[1][1]=Cmake(q_r,q_i);

		apply_gate(Qm, controls+1, nq-(1+2*nqubits-i));

		line[nq-(1+2*nqubits-i)] = -1;
		if(c1 >= 0) {
			line[nq-c1] = -1;
		}
		if(c2 >= 0) {
			line[nq-c2] = -1;
		}
	}
}

void add_phi_inv(int n, int c1, int c2) {
	int q;
	int controls = 0;
	if(c1 >= 0) controls++;
	if(c2 >= 0) controls++;
	for(int i = nqubits; i>= 0; --i) {
		q = 1;
		int fac = 0;
		for(int j = i; j >= 0; --j) {
			if((n >> j)&1) {
				fac |= 1;
			}
			fac = fac << 1;
			q*=2;
		}
		if(c1 >= 0) {
			line[nq-c1] = 1;
		}
		if(c2 >= 0) {
			line[nq-c2] = 1;
		}

		line[nq-(1+2*nqubits-i)] = 2;

		q_r = QMDDcos(fac, -q);
		q_i = QMDDsin(fac, -q);

		Qm[1][1]=Cmake(q_r,q_i);
		apply_gate(Qm, controls+1, nq-(1+2*nqubits-i));

		line[nq-(1+2*nqubits-i)] = -1;
		if(c1 >= 0) {
			line[nq-c1] = -1;
		}
		if(c2 >= 0) {
			line[nq-c2] = -1;
		}
	}
}

void qft() {
	int q;
	for(int i = nqubits+1; i < 2*nqubits+2; i++) {
		line[nq-i] = 2;
		apply_gate(Hm, 1, nq-i);
		line[nq-i] = -1;

		q = 2;
		for(int j = i+1; j < 2*nqubits+2; j++) {
			line[nq-j] = 1;
			line[nq-i] = 2;

			q_r = QMDDcos(1, q);
			q_i = QMDDsin(1, q);
			Qm[1][1]=Cmake(q_r, q_i);
			apply_gate(Qm, 2, nq-i);

			line[nq-j] = -1;
			line[nq-i] = -1;
			q *= 2;
		}
	}
}

void qft_inv() {
	int q, qq;
	for(int i = 2*nqubits+1; i >= nqubits+1; i--) {
		q = 2;
		for(int j = i+1; j < 2*nqubits+2; j++) {
			qq = q;
			line[nq-j] = 1;
			line[nq-i] = 2;

			q_r = QMDDcos(1, -qq);
			q_i = QMDDsin(1, -qq);

			Qm[1][1]=Cmake(q_r, q_i);

			apply_gate(Qm, 2, nq-i);

			line[nq-j] = -1;
			line[nq-i] = -1;
			q *= 2;
		}
		line[nq-i] = 2;
		apply_gate(Hm, 1, nq-i);
		line[nq-i] = -1;
	}
}

void mod_add_phi(int a, int N, int c1, int c2) {
	add_phi(a, c1, c2);
	add_phi_inv(N,-1,-1);
	qft_inv();

	line[nq-(nqubits+1)] = 1;
	line[nq-(2*nqubits+2)] = 2;

	apply_gate(Nm, 2, nq-(2*nqubits+2));

	line[nq-(nqubits+1)] = -1;
	line[nq-(2*nqubits+2)] = -1;

	qft();
	add_phi(N, 2*nqubits+2, -1);
	add_phi_inv(a,c1,c2);
	qft_inv();

	line[nq-(nqubits+1)] = 0;
	line[nq-(2*nqubits+2)] = 2;
	apply_gate(Nm, 2, nq-(2*nqubits+2));
	line[nq-(nqubits+1)] = -1;
	line[nq-(2*nqubits+2)] = -1;

	qft();
	add_phi(a, c1, c2);
}

void mod_add_phi_inv(int a, int N, int c1, int c2) {
	add_phi_inv(a, c1, c2);
	qft_inv();

	line[nq-(nqubits+1)] = 0;
	line[nq-(2*nqubits+2)] = 2;

	apply_gate(Nm, 2, nq-(2*nqubits+2));

	line[nq-(nqubits+1)] = -1;
	line[nq-(2*nqubits+2)] = -1;
	qft();
	add_phi(a,c1,c2);
	add_phi_inv(N, 2*nqubits+2, -1);
	qft_inv();
	line[nq-(nqubits+1)] = 1;
	line[nq-(2*nqubits+2)] = 2;

	apply_gate(Nm,2,  nq-(2*nqubits+2));

	line[nq-(nqubits+1)] = -1;
	line[nq-(2*nqubits+2)] = -1;

	qft();
	add_phi(N,-1,-1);
	add_phi_inv(a, c1, c2);
}

void cmult(int a, int N, int c) {
	qft();

	int t = a;
	for(int i = nqubits; i >= 1; i--) {
		mod_add_phi(t, N, i, c);
		t = (2*t) % N;
	}
	qft_inv();
}

void cmult_inv(int a, int N, int c) {
	qft();
	int t = a;
	for(int i = nqubits; i >= 1; i--) {
		mod_add_phi_inv(t, N, i, c);
		t = (2*t) % N;
	}
	qft_inv();
}

void u_a(int a, int N, int c) {
	cmult(a, N, c);

	for(int i = 0; i < nqubits; i++) {
		line[nq-(nqubits+2+i)] = 1;
		line[nq-(i+1)] = 2;
		apply_gate(Nm, 2, nq-(i+1));

		line[nq-(nqubits+2+i)] = 2;
		line[nq-(i+1)] = 1;
		line[nq-c] = 1;
		apply_gate(Nm, 3, nq-(nqubits+2+i));
		line[nq-(nqubits+2+i)] = 1;
		line[nq-(i+1)] = 2;
		line[nq-c] = -1;

		apply_gate(Nm, 2, nq-(i+1));

		line[nq-(nqubits+2+i)] = -1;
		line[nq-(i+1)] = -1;
	}
	cmult_inv(inverse_mod(a, N), N, c);
}

unsigned long gcd(unsigned long a, unsigned long b) {
	unsigned long c;
    while ( a != 0 ) {
    	c = a; a = b%a;  b = c;
    }
    return b;
}

unsigned long modpow(unsigned long base, unsigned long exp, unsigned long modulus) {
  base %= modulus;
  unsigned long result = 1;
  while (exp > 0) {
    if (exp & 1) result = (result * base) % modulus;
    base = (base * base) % modulus;
    exp >>= 1;
  }
  return result;
}

void shor(unsigned long n) {
    cout << "Simulate Shor's algorithm for n=" << n;

    //Prepare initial state of the circuit
    circ.nancillary=circ.ngarbage=0;

    circ.n=ceil(log2(n))*2+3;
    cout << " (requires " << circ.n << " qubits):" << endl;

    for(int p=circ.n-1;p>=0;p--) {
    	snprintf(circ.line[p].variable, MAXSTRLEN, "x%d", p);
    	snprintf(circ.line[p].input, MAXSTRLEN, "i%d", p);
    	snprintf(circ.line[p].output, MAXSTRLEN, "o%d", p);
        circ.line[p].ancillary = (p==circ.n-1-ceil(log2(n))) ? '1' : '0';
        circ.line[p].garbage='-';
    }

    for(int i=0;i<circ.n;i++) {
    	circ.inperm[i]=i;
    }

	QMDDedge e = QMDDone;
	QMDDedge edges[4];
	edges[1]=edges[3]=QMDDzero;

	for(int p=0;p<circ.n;p++) {
		if(circ.line[p].ancillary=='0') {
			edges[0] = e;
			edges[2] = QMDDzero;
		} else if(circ.line[p].ancillary=='1') {
			edges[0] = QMDDzero;
			edges[2] = e;
		} else {
			edges[0]=edges[2]=e;
		}
		e = QMDDmakeNonterminal(p, edges);
	}
	QMDDincref(e);
	circ.e = e;

	nqubits = ceil(log2(n));
	nq = nqubits*2+2;

	line = new int[2*nqubits+3];
	for(int j=0;j<2*nqubits+3;j++) {
		line[j] = -1;
	}


	int a;
	do {
		a = rand() % n;
	} while (gcd(a, n) != 1 || a == 1);

	cout << "Find a coprime to N: " << a << endl;

	measurements = new int[2*nqubits];

	int* as = new int[2*nqubits];
	as[2*nqubits-1] = a;
	unsigned long new_a = a;
	for(int i = 2*nqubits-2; i >= 0; i--) {
		new_a = new_a * new_a;
		new_a = new_a % n;
		as[i] = new_a;
	}

	//Start Shor's algorithm
	for(int i = 0; i < 2*nqubits; i++) {

		line[nq] = 2;
		apply_gate(Hm, 1, nq);
		line[nq] = -1;

		u_a(as[i], n, 0);

		line[nq] = 2;
		int q = 2;

		for(int j =i-1; j >= 0; j--) {
			if(measurements[j] == 1) {
				q_r = QMDDcos(1, -q);
				q_i = QMDDsin(1, -q);

				Qm[1][1]=Cmake(q_r, q_i);
				apply_gate(Qm, 1, nq);
			}
			q = q << 1;
		}

		apply_gate(Hm, 1, nq);
		measurements[i] = QMDDmeasureOne(nq);
		QMDDgarbageCollect();
		vector<QMDDedge> v;
		v.push_back(circ.e);
		cleanCtable(v);
		v.clear();

		if(measurements[i] == 1) {

			apply_gate(Nm, 1, nq);
		}
		line[nq] = -1;
	}
	unsigned long res = 0;
	cout << "measurement: ";
	for(int i=0; i<2*nqubits; i++) {
		cout << measurements[2*nqubits-1-i];
		res = (res << 1) + measurements[2*nqubits-1-i];
	}
	cout << " = " << res << endl;
	unsigned long denom = 1 << (2*nqubits);

	//Post processing (non-quantum part)

	if(res == 0) {
		cout << "Factorization failed (measured 0)!" << endl;
	} else {
	cout << "Continued fraction expansion of " << res << "/" << denom << " = " << flush;
	vector<int> cf;

	unsigned long old_res = res;
	unsigned long old_denom = denom;
	while(res != 0) {
		cf.push_back(denom/res);
		unsigned long tmp = denom%res;
		denom = res;
		res = tmp;
	}

	for(int i = 0; i < cf.size(); i++) {
		cout << cf[i] << " ";
	} cout << endl;

	int success = 0;
	for(int i=0; i < cf.size(); i++) {
		//determine candidate
		unsigned long denominator = cf[i];
		unsigned long numerator = 1;

		for(int j=i-1; j >= 0; j--) {
			unsigned long tmp = numerator + cf[j]*denominator;
			numerator = denominator;
			denominator = tmp;
		}
		cout << "  Candidate " << numerator << "/" << denominator << ": ";
		if(denominator > n) {
			cout << " denominator too large (greater than " << n << ")!" << endl;
			success = 0;
			cout << "Factorization failed!" << endl;break;
		} else {
			double delta = (double)old_res / (double)old_denom - (double)numerator / (double) denominator;
			if(fabs(delta) < 1.0/(2.0*old_denom)) {
				if(modpow(a, denominator, n) == 1) {
					cout << "found period = " << denominator << endl;
					if(denominator & 1) {
						cout << "Factorization failed (period is odd)!" << endl;
					} else {
						cout << "Factorization succeeded! Non-trivial factors are: " << endl;
						unsigned long f1, f2;
						f1 = modpow(a, denominator>>1, n);
						f2 = (f1+1)%n;
						f1 = (f1 == 0) ? n-1 : f1-1;
						f1 = gcd(f1, n);
						f2 = gcd(f2, n);
						cout << "  -- gcd(" << n << "^(" << denominator << "/2)-1" << "," << n << ") = " << f1 << endl;
						cout << "  -- gcd(" << n << "^(" << denominator << "/2)+1" << "," << n << ") = " << f2 << endl;
					}

					break;
				} else {
					cout << "failed" << endl;
				}
			} else {
				cout << "delta is too big (" << delta << ")"<<endl;
			}
		}
	}
	}
	delete[] measurements;
	delete[] line;

}

QMDDedge randomOracle(unsigned long n, unsigned long* expected) {
	cout << "generate random " << n << "-bit oracle: " << flush;
	*expected = 0;
	for(int p=0;p < n;p++) {
		if(rand() %2 == 0) {
			cout << "0";
			line[n-p] = 0;
			*expected = (*expected) << 1;
		} else {
			cout << "1";
			line[n-p] = 1;
			*expected = ((*expected) << 1) | 1;
		}
	}

	cout << " = " << *expected << endl;

	line[0] = 2;
	QMDDedge f = QMDDmvlgate(Nm, n+1, line);
	for(int i=0; i<n+1; i++) {
		line[i] = -1;
	}
	return f;
}

QMDDedge groverIteration(QMDDedge oracle) {

	QMDDedge e = circ.e;
	circ.e = oracle;
	QMDDincref(circ.e);
	for(int i =0; i < circ.n-1; i++) {
		line[circ.n-1-i] = 2;
		apply_gate(Hm, 1, circ.n-1-i);
		line[circ.n-1-i] = -1;
	}
	for(int i =0; i < circ.n-1; i++) {
		line[circ.n-1-i] = 2;
		apply_gate(Nm, 1, circ.n-1-i);
		line[circ.n-1-i] = -1;
	}
	line[1] = 2;
	apply_gate(Hm, 1, 1);
	line[1] = -1;

	for(int i=0; i < circ.n-2; i++) {
		line[circ.n-1-i] = 1;
	}
	line[1] = 2;
	apply_gate(Nm, circ.n-1, 1);
	for(int i=0; i < circ.n; i++) {
		line[i] = -1;
	}
	line[1] = 2;
	apply_gate(Hm, 1, 1);

	line[1] = -1;
	for(int i =0; i < circ.n-1; i++) {
		line[circ.n-1-i] = 2;
		apply_gate(Nm, 1, circ.n-1-i);
		line[circ.n-1-i] = -1;
	}
	for(int i =0; i < circ.n-1; i++) {
		line[circ.n-1-i] = 2;
		apply_gate(Hm, 1, circ.n-1-i);
		line[circ.n-1-i] = -1;
	}
	QMDDedge t = circ.e;
	circ.e = e;
	QMDDdecref(t);
	return t;
}

void grover(unsigned long n) {

	cout << "Simulate Grover's search algorithm for " << n << "-Bit number ";

    circ.nancillary=circ.ngarbage=0;

    circ.n=n+1;
    cout << " (requires " << circ.n << " qubits):" << endl;

    for(int p=circ.n-1;p>=0;p--) {
    	snprintf(circ.line[p].variable, MAXSTRLEN, "x%d", p);
    	snprintf(circ.line[p].input, MAXSTRLEN, "i%d", p);
    	snprintf(circ.line[p].output, MAXSTRLEN, "o%d", p);
        circ.line[p].ancillary = (p==0) ? '1' : '0';
        circ.line[p].garbage='-';
    }

    for(int i=0;i<circ.n;i++) {
    	circ.inperm[i]=i;
    }

	QMDDedge e = QMDDone;
	QMDDedge edges[4];
	edges[1]=edges[3]=QMDDzero;

	for(int p=0;p<circ.n;p++) {
		if(circ.line[p].ancillary=='0') {
			edges[0] = e;
			edges[2] = QMDDzero;
		} else if(circ.line[p].ancillary=='1') {
			edges[0] = QMDDzero;
			edges[2] = e;
		} else {
			edges[0]=edges[2]=e;
		}
		e = QMDDmakeNonterminal(p, edges);
	}
	QMDDincref(e);
	circ.e = e;

	line = new int[circ.n];
	measurements = new int[circ.n];
	for(int j=0;j<circ.n;j++) {
		line[j] = -1;
	}

	unsigned long iterations;
	if((circ.n-1)%2 == 0) {
		iterations = pow(2, (circ.n-1)/2);
	} else {
		iterations = floor(pow(2, (circ.n-2)/2) * sqrt(2));
	}

	unsigned long expected;
	QMDDedge t = randomOracle(n, &expected);
	QMDDincref(t);

	for(int i =0; i < circ.n; i++) {
		line[circ.n-1-i] = 2;
		apply_gate(Hm, 1, circ.n-1-i);
		line[circ.n-1-i] = -1;
	}

	QMDDedge gi = groverIteration(t);
	QMDDincref(gi);
	QMDDdecref(t);

	cout << "perform " << iterations << " Grover iterations" << endl;

	for(int i = 0; i < iterations; i++) {
		t = QMDDmultiply(gi, circ.e);
		QMDDincref(t);
		QMDDdecref(circ.e);
		if(ActiveNodeCount > max_active) {
			max_active = ActiveNodeCount;
		}
		circ.e = t;
		QMDDgarbageCollect();
		if((i+1) %1000 == 0) {
			cout << "  -- performed " << (i+1) << " of " << iterations << " grover iterations" << endl;
			vector<QMDDedge> v;

			v.push_back(circ.e);
			v.push_back(gi);
			cleanCtable(v);
			v.clear();
		}
	}

	QMDDmeasureAll();
	cout << "Measured: ";
    unsigned long measured = 0;
	for(int i = circ.n-1; i > 0; --i) {
		cout << measurements[i];
		measured = (measured << 1) | measurements[i];
	}
	cout << " = " << measured << " (expected " << expected << ")" << endl;

	gatecount = (1+(circ.n-1)*4+3)*iterations + circ.n;
	delete[] measurements;
}


int main(int argc, char** argv) {

	namespace po = boost::program_options;
	po::options_description description("Allowed options");
	description.add_options()
	    ("help", "produce help message")
		("set_seed", po::value<unsigned long>(), "set seed for random number generator")
	    ("simulate_circuit", po::value<string>(), "simulate a quantum circuit given in .real format")
		("simulate_shor", po::value<unsigned long>(), "simulate Shor's algorithm for <arg>")
		("simulate_grover", po::value<unsigned long>(), "simulate Groover's algorithm with <arg> qubits")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, description), vm);
	po::notify(vm);

	if (vm.count("help")) {
	    cout << description << "\n";
	    return 1;
	}

	max_active = 0;
	gatecount = 0;
	epsilon = mpreal(0.01);

	unsigned long seed = time(NULL);
	if (vm.count("set_seed")) {
		seed = vm["set_seed"].as<unsigned long>();
	}

	cout << "Seed for random number generator is " << seed << endl;
	srand(seed);

	QMDDinit(0);

    struct timeval t1, t2;
    double elapsedTime;

	    // start timer
	gettimeofday(&t1, NULL);

	if (vm.count("simulate_circuit")) {
		QMDDsimulate(vm["simulate_circuit"].as<string>());
	} else if (vm.count("simulate_shor")) {
		shor(vm["simulate_shor"].as<unsigned long>());
	} else if (vm.count("simulate_grover")) {
	    grover(vm["simulate_grover"].as<unsigned long>());
	} else {
		cout << description << "\n";
	    return 1;
	}

	gettimeofday(&t2, NULL);

	// compute and print the elapsed time in millisec
    elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;

    cout << endl << "SIMULATION STATS: " << endl;
    cout << "  Number of applied gates: " << gatecount << endl;
    cout << "  Simulation time: " << elapsedTime/1000 << " seconds" << endl;
    cout << "  Maximal size of QMDD (number of nodes) during simulation: " << max_active << endl;

	return 0;
}
