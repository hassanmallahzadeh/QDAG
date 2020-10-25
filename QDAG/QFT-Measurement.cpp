//
//  QFT-Measurement.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-07-07.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//
#include "QFT-Measurement.hpp"
#include "util.h"
   
Measurement::Measurement(dd::Package *dd) {
    this->dd = dd;
}
/// The 'gate' to Measurment class.
///
/// @return measurement outcome (0 to radix)
/// @param e_root root of bdd for measurement
/// @param n number of qubits
/// @param n_target targeted qubit index
/// @param unrg random generator engine
int Measurement::Measure(dd::Edge &e_root, int n/*for StateCollapseMatMul*/, int n_target, engine& unrg){
    this->e_root = e_root;
    PopulateUpProbDiagram(e_root);
    PopulateDownProbDiagram(e_root);
    array<fp,dd::RADIX> probs = QubitMeasurementProbs(n_target);
    int res = QubitMeasurementOutcome(probs, unrg);
    Mqinfo mqinfo = {n_target, res, probs[res]};
    StateCollapseRestrict(mqinfo);
    return res;
}
/// Upstream probabilities for all dd nodes. recureively each node usp is its childs usp weighted by their edge weights (val squared) each.
/// @param edge root edge
fp Measurement::PopulateUpProbDiagram(const dd::Edge& edge){
    
    if(upmap.find(edge.p) != upmap.end()){// already calculated so just return from umap
        return upmap[edge.p];
    }
    if(edge.p == dd::Package::DDone.p){// base case
        upmap.insert(make_pair(edge.p, 1));
        return 1;
    }
    if (edge.p == dd::Package::DDzero.p){// base case
        upmap.insert(make_pair(edge.p, 0));
        return 0;
    }
    {// else: not already calculated
        fp val = 0;
        for (int k = 0; k < dd::NEDGE; ++k) {
            if (k % dd::RADIX != 0) continue;// not a vector edge.
            val += CN::mag2(edge.p->e[k].w) * PopulateUpProbDiagram(edge.p->e[k]);
        }
        upmap.insert(make_pair(edge.p,val));
    }
    return upmap[edge.p];
}
/// Canlculate downstream probabilities. For each node dsp is sum of parent's dsp weighted by corresponding edge weights (val squared) each.
/// @param edge root edge
void Measurement::PopulateDownProbDiagram(const dd::Edge &edge){
    list<dd::NodePtr> l;
    assert(edge.p);// make sure not null
    dd::NodePtr curnode = edge.p;
    l.push_back(curnode);
    downmap.insert(make_pair(curnode,CN::mag2(edge.w)));
    // DDone.p, DDone.pnot are not neccessary to add to downmap. added for sake of tradition.
    downmap.insert(make_pair(dd::Package::DDone.p, 1));
    downmap.insert(make_pair(dd::Package::DDzero.p, 0));
    while(!l.empty()){// breath first search
        curnode = l.front();
        l.pop_front();
        
        for(int i = 0; i < dd::NEDGE; ++i){
            if (i % dd::RADIX != 0) continue;// not a vector edge
            if(curnode->e[i].p == dd::Package::DDone.p){// base case
                continue;
            }
            if (curnode->e[i].p == dd::Package::DDzero.p){// base case
                continue;
            }
            if(downmap.find(curnode->e[i].p) != downmap.end()){// if already in downmap modify value
                downmap[curnode->e[i].p] += downmap[curnode] * CN::mag2(curnode->e[i].w);
            }
            else{// else insert
                downmap.insert(make_pair(curnode->e[i].p, downmap[curnode] * CN::mag2(curnode->e[i].w)));
                l.push_back(curnode->e[i].p);
            }
        }
        
    }
}

void Measurement::ExtractedQubitMeasProbsOnLayer(array<fp, dd::RADIX> &a, dd::NodePtr curnode, int i) {
    a[i / dd::RADIX] +=  downmap[curnode] * CN::mag2(curnode->e[i].w)* upmap[curnode->e[i].p];
}

void Measurement::ExtractedQubitMeasProbsOnLayer(dd::NodePtr curnode, int i, list<dd::NodePtr> &l) {
    l.push_back(curnode->e[i].p);
}

void Measurement::ExtractedQubitMeasProbsLayer(array<fp, dd::RADIX> &a, dd::NodePtr curnode, int i, dd::NodePtr lastnode) {
    assert(0);//temp
    assert(curnode->v >= 0);//make sure we dont get here because of terminal nodes. they shall be blocked by if(traverseset.find(curnode) != traverseset.cend())
    a[i / dd::RADIX] +=  downmap[lastnode] * 0.5;
}

/// Measure probability of each outcome to measure qubit v.
/// @param v qubit index
array<fp,dd::RADIX> Measurement::QubitMeasurementProbs(int v) {
    //TODO: last prob in array: 1-(sum of others) do after test.
    list<dd::NodePtr> l;
    array<fp,dd::RADIX> a = {};// holds probabilities.
    dd::NodePtr curnode = e_root.p;
    dd::NodePtr lastnode = nullptr;// track parent.
    traverseset.insert(dd::Package::terminalNode);//base case
    l.push_back(curnode);
    while(!l.empty()){// breath first search
        
        curnode = l.front();
        l.pop_front();
        if(traverseset.find(curnode) != traverseset.cend()) continue;
        
        traverseset.insert(curnode);
        
        for(int i = 0; i < dd::NEDGE; ++i){
            if (i % dd::RADIX != 0) continue;// not a vector edge
            if(curnode->v > v){//above targer.
                ExtractedQubitMeasProbsOnLayer(curnode, i, l);
            }
            else if(curnode->v == v){// in target layer
                ExtractedQubitMeasProbsOnLayer(a, curnode, i);
            }
            else{////  below target level:  if node missed due to reduced graph. NOTE: won't happen with current IIC-JKU package.
                ExtractedQubitMeasProbsLayer(a, curnode, i, lastnode);
            }
        }
    }
    return a;
}
int Measurement::QubitMeasurementOutcome(array<fp, dd::RADIX> a, engine& unrg) {
    std::uniform_real_distribution<fp> dis(0.0, 1.0);
    array<fp, dd::RADIX> aincsum;// array incremental summations :)
    aincsum[0] = a[0];

    for(int i = 1; i < dd::RADIX; ++i){
        aincsum[i] = aincsum[i-1] + a[i];
    }
    assert(aincsum[dd::RADIX-1] > 0.99 && aincsum[dd::RADIX-1] < 1.01);//if state not normalized to 1.
    
    fp randnum = dis(unrg);
    for (int i = 0; i < dd::RADIX; ++i){
        if(randnum <= aincsum[i]){
            return i;
        }
    }
    assert(0);//should not get here.
    return -1;//error
}

/// Collapsed state after bits, use matrix multiplication. Wont work for radix != 2
/// @param mqinfo measured qubits info
/// @param n total number of quibits
void Measurement::StateCollapseMatMul(Mqinfo mqinfo, int n) {
    assert(0);//function does not work. Likely since the matrix is not unitary.(returns identity)
    assert(dd::RADIX == 2);
    short* line = new short[n]{};
    for(int i = 0; i < n; ++i){
        line[i] = -1;
    }
    dd::Edge e_gate;
    line[mqinfo.ix] = 2;
    if(mqinfo.val == 0)
        e_gate = dd->makeGateDD(N0mat, n, line);
    else
        e_gate = dd->makeGateDD(N1mat, n, line);
    
    dd::ComplexValue c{sqrt(1/mqinfo.pr), 0.0 };//used to normalize the state.
    dd::Complex cx = dd->cn.getCachedComplex(c.r,c.i);
    e_gate.w = dd->cn.mulCached(e_gate.w, cx);
    e_root = dd->multiply(e_gate, e_root);
    line[mqinfo.ix] = -1;//reset
    delete[] line;
}

void Measurement::ExtractedBelowLayerOnCollapse(dd::NodePtr curnode, const dd::Complex &cx, int j, Mqinfo &mqinfo) {
    assert(0);//temp
    assert(curnode->v >= 0);//make sure we dont get here because of terminal
    dd::Edge e[dd::NEDGE] = {};
    for (int k = 0; k < dd::NEDGE; ++k){
        if(k == mqinfo.val * dd::RADIX){// equal to measured value.
            e[k].w = cx;
            e[k].p = curnode;
        }
        else{// else set to zero.
            e[k].w = CN::ZERO;
            e[k].p = dd::Package::DDzero.p;
        }
    }
    dd::Edge ge  = dd->makeNonterminal(mqinfo.ix, e);// ge: generated edge.
    curnode->e[j] = ge;
}

void Measurement::ExtractedOnLayerStateCollapse(dd::NodePtr curnode, const dd::Complex &cx, int j, Mqinfo &mqifo) {
    if(j == mqifo.val * dd::RADIX){// equal to measured value.
        if(!CN::equalsZero(curnode->e[j].w)){
            if(!CN::equalsZero(curnode->e[j].w)){
                curnode->e[j].w = dd->cn.mulCached(curnode->e[j].w, cx);
            }
            //                static bool applied = false;//TODO: This way causes root weight to become 1 when Measure() returns. Ask Alan...
            //                if(applied == false){
            //                    e_root.w = dd->cn.mulCached(e_root.w, cx);
            //                    applied = true;
            //                }
        }
    }
    else{// else set to zero.
        curnode->e[j].p = dd::Package::DDzero.p;
        curnode->e[j].w = CN::ZERO;
    }
}

void Measurement::ExtractedAboveLayerOnCollapse(dd::NodePtr curnode, int j, list<dd::NodePtr> &l) {
    if(curnode->e[j].p->v >= 0){
        l.push_back(curnode->e[j].p);
    }
}

/// Collapsed state after bits, use direct RESTRICT algorithm. "Bryant's Restrict algortithm".
/// @param mqinfo measured qubits info
void Measurement::StateCollapseRestrict(Mqinfo mqinfo) {
    traverseset.clear();//clear act of QubitMeasurementProbs;
        
        dd::ComplexValue c{sqrt(1/mqinfo.pr), 0.0 };//used to normalize the state.
        dd::Complex cx = dd->cn.getCachedComplex(c.r,c.i);
        
        list<dd::NodePtr> l;
        dd::NodePtr curnode = e_root.p;
        l.push_back(curnode);
        traverseset.insert(dd::Package::DDone.p);//base case,likely not needed.
        traverseset.insert(dd::Package::DDzero.p);//base case
        while(!l.empty()){// breath first search
            curnode = l.front();
            l.pop_front();
            if(traverseset.find(curnode) != traverseset.cend())
                continue;
            traverseset.insert(curnode);
            for(int j = 0; j < dd::NEDGE; ++j){
                if(j % dd::RADIX != 0) continue;//not a vector edge.
                if(curnode->v > mqinfo.ix){//above target.
                    ExtractedAboveLayerOnCollapse(curnode, j, l);
                }
                else if(curnode->v == mqinfo.ix){//in target layer
                    
                    ExtractedOnLayerStateCollapse(curnode, cx, j, mqinfo);
                }
                else{////  below target level: if node missed due to reduced graph. NOTE: won't happen with current IIC-JKU package.
                    ExtractedBelowLayerOnCollapse(curnode, cx, j, mqinfo);
                }
            }
        }
  //  dd->garbageCollect();
}

