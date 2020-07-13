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
#include <unordered_map>
#include <stdio.h>
using std::unordered_map;
using std::make_pair;
class Measurement{
private:
    dd::Edge e_root;
    dd::Package* dd = nullptr;
    unordered_map<dd::NodePtr, fp> upmap;
    unordered_map<dd::NodePtr, fp> downmap;
    fp PopulateUpProbDiagram(const dd::Edge &edge);
    void PopulateDownProbDiagram(const dd::Edge &edge);
    array<fp,dd::RADIX> QubitMeasurementProbs(int v);
    int QubitMeasurementOutcome(array<fp,dd::RADIX>);
public:
    Measurement(dd::Package* dd, const dd::Edge &e_root);
    ~Measurement();
};

#endif /* QFT_Measurement_hpp */
