//
//  QFT-Measurement.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-07-07.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//
#include "commonheaders.h"
#include "QFT-Measurement.hpp"
//FIXME: replace all at() with [] after test.

Measurement::Measurement(dd::Package *dd, const dd::Edge &e_root) {
    this->dd = dd;
    this->e_root = e_root;
}
fp Measurement::PopulateUpProbDiagram(const dd::Edge& edge){
    
    if(upmap.find(edge.p) != upmap.end()){//HM: already calculated so just return from umap
        return upmap.at(edge.p);
    }
    if(edge.p == dd::Package::DDone.p){//HM: base case
        upmap.insert(make_pair(edge.p, 1));
        return 1;
    }
    if (edge.p == dd::Package::DDzero.p){//HM: base case
        upmap.insert(make_pair(edge.p, 0));
        return 0;
    }
    {//HM: else: not already calculated
        fp val = 0;
        for (int k = 0; k < dd::NEDGE; ++k) {
            if (k % dd::RADIX != 0) continue;//HM: not a vector edge.
            val += CN::mag2(edge.p->e[k].w) * PopulateUpProbDiagram(edge.p->e[k]);
        }
        upmap.insert(make_pair(edge.p,val));
    }
    return upmap.at(edge.p);
}
void Measurement::PopulateDownProbDiagram(const dd::Edge &edge){
    list<dd::NodePtr> l;
    assert(edge.p);//HM: make sure not null
    dd::NodePtr curnode = edge.p;
    l.push_back(curnode);
    downmap.insert(make_pair(curnode,CN::mag2(edge.w)));
    while(!l.empty()){//HM: breath first search
        curnode = l.front();
        l.pop_front();
        for(int i = 0; i < dd::NEDGE; ++i){
            if (i % dd::RADIX != 0) continue;//HM: not a vector edge
            if(downmap.find(curnode->e[i].p) != downmap.end()){//HM: if in downmap modify value
                downmap.at(curnode->e[i].p) += downmap.at(curnode) * CN::mag2(edge.w);
            }
            else{
                downmap.insert(make_pair(curnode->e[i].p, downmap.at(curnode) * CN::mag2(edge.w)));
            }
            if(!curnode->e[i].p) continue;//HM: base case. will be null for terminal childs.
            l.push_back(curnode->e[i].p);
        }
        
    }
}

array<fp,dd::RADIX> Measurement::QubitMeasurementProbs(int v) {
    //TODO: last prob in array: 1-(sum of others) do after test.
    list<dd::NodePtr> l;
    array<fp,dd::RADIX> a = {};//HM: holds probabilities.
    dd::NodePtr curnode = e_root.p;
    l.push_back(curnode);
    while(!l.empty()){//HM: breath first search
        curnode = l.front();
        l.pop_front();
        
        for(int i = 0; i < dd::NEDGE; ++i){
            if (i % dd::RADIX != 0) continue;//HM: not a vector edge
            if(curnode->v == v){
                if (i % dd::RADIX != 0) continue;//HM: not a vector edge
                a[i] +=  downmap.at(curnode) * CN::mag2(curnode->e[i].w)* upmap.at(curnode->e[i].p)/upmap.at(curnode);
            }
            else if(curnode->v > v){//HM: if node missed due to reduced graph.
                //TODO: this should not happen as there is a node in each layer each path. Should be included anyway as we stive for reduced bdd.
                //                for (int i = 0; i < dd::RADIX; ++i){
                //                    a[i] += 0.5 * downmap.at(curnode->e[i].p);
                //                }
            }
            else{
                l.push_back(curnode->e[i].p);//HM: nodes on upper layers.
            }
        }
    }
    return a;
}
//TODO: replace with C++11 random generator.
int Measurement::QubitMeasurementOutcome(array<fp, dd::RADIX> a) {
    array<fp, dd::RADIX> aincsum;//HM: array incremental summations :)
    aincsum[0] = a[0];
    for(int i = 1; i < dd::RADIX; ++i){
        aincsum[i] = aincsum[i-1] + a[i];
    }
    srand (static_cast<unsigned int>(time(nullptr)));
    fp randnum = (double) rand() / (RAND_MAX);
    for (int i = 0; dd::RADIX; ++i){
        if(randnum <= aincsum[i])
            return i;
    }
}


