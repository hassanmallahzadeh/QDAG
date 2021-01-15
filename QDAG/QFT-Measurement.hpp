//
//  QFT-Measurement.hpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-07-07.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#ifndef QFT_Measurement_hpp
#define QFT_Measurement_hpp
#include "IIC-JKU/DDpackage.h"
#include "IIC-JKU/DDcomplex.h"
#include "commonheaders.h"
#include <unordered_map>
#include <unordered_set>
#include <stdio.h>
using std::unordered_map;
using std::unordered_set;
using std::make_pair;
class Measurement{
private:
    dd::Edge e_root;
    dd::Package* dd = nullptr;
    unordered_map<dd::NodePtr, fp> upmap;
    unordered_map<dd::NodePtr, fp> downmap;
    unordered_set<dd::NodePtr> traverseset;
    fp PopulateUpProbDiagram(const dd::Edge &edge);
    void PopulateDownProbDiagram(const dd::Edge &edge);

    struct Mqinfo{//measured qubit info
        int ix = -1;// index
        int val = -1;// measured value
        fp pr = -1;// probability of measured outcome(before conductiong measurement).
    };
    
    void ExtractedBelowLayerOnCollapse(dd::NodePtr curnode, const dd::Complex &cx, int j, Mqinfo &mqi);
    void ExtractedOnLayerStateCollapse(dd::NodePtr curnode, const dd::Complex &cx, int j, Mqinfo &mqi);
    void ExtractedAboveLayerOnCollapse(dd::NodePtr curnode, int j, list<dd::NodePtr> &l);
    
    void ExtractedQubitMeasProbsOnLayer(array<fp, dd::RADIX> &a, dd::NodePtr curnode, int i);
    
    void ExtractedQubitMeasProbsOnLayer(dd::NodePtr curnode, int i, list<dd::NodePtr> &l);
    
    void ExtractedQubitMeasProbsLayer(array<fp, dd::RADIX> &a, dd::NodePtr curnode, int i, dd::NodePtr lastnode);
    
    array<fp,dd::RADIX> QubitMeasurementProbs(int v);
    int QubitMeasurementOutcome(array<fp, dd::RADIX>, engine&);
    void StateCollapseRestrict(Mqinfo mqinfo);
    void StateCollapseMatMul(Mqinfo mqinfo, int n);
    std::uniform_real_distribution<fp> dis;
public:
    Measurement(dd::Package* dd);
    int Measure(dd::Edge &edge, int n, int n_target, engine& unrg);
};

#endif /* QFT_Measurement_hpp */
